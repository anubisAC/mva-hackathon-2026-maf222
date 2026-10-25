"""Selected-gene gnomAD v4 audit with local matching and independent ablations.

Only public gene IDs are sent to gnomAD; patient alleles and genotypes stay local.
The gene selection is analysis-informed. This is not genome-wide reannotation.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import time
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .genomewide import PhenotypeIndex, build_hypotheses, rank_genes, remove_clinvar_evidence

GROUPS = {'afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'sas'}
CONFIG = Path('config/gnomad_frequency_review.json')
# Frozen public gene query; schema inspected at gnomad-browser commit
# 6836f67f92df3184e00e99f3cf97ad6283647873. No patient alleles are queried.
PUBLIC_QUERY = 'query PublicGeneFrequency($gene: String!) { gene(gene_id:$gene, reference_genome:GRCh38) { gene_id symbol chrom start stop variants(dataset:gnomad_r4) { variant_id reference_genome chrom pos ref alt flags exome { ac an homozygote_count filters flags populations { id ac an homozygote_count } fafmax { faf95_max faf95_max_gen_anc faf99_max faf99_max_gen_anc } } genome { ac an homozygote_count filters flags populations { id ac an homozygote_count } fafmax { faf95_max faf95_max_gen_anc faf99_max faf99_max_gen_anc } } joint { ac an homozygote_count filters populations { id ac an homozygote_count } fafmax { faf95_max faf95_max_gen_anc faf99_max faf99_max_gen_anc } } transcript_consequence { gene_id transcript_id hgvsc hgvsp consequence_terms } } coverage(dataset:gnomad_r4) { exome { pos mean median over_20 } genome { pos mean median over_20 } } } }'

def read(path):
    path=Path(path)
    with (gzip.open(path,'rt',encoding='utf-8') if path.suffix=='.gz' else path.open(encoding='utf-8')) as f:
        return json.load(f)

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def count_af(ac,an):
    if isinstance(ac,bool) or isinstance(an,bool): return None
    if not isinstance(ac,int) or not isinstance(an,int) or an<=0 or ac<0 or ac>an: return None
    return ac/an

def passed(data):
    return bool(data) and isinstance(data.get('filters'),list) and set(data['filters']) <= {'PASS'} and count_af(data.get('ac'),data.get('an')) is not None

def normalize_variant(chrom,pos,ref,alt):
    """Minimal representation only; no repeat-aware left alignment without FASTA."""
    chrom=str(chrom).removeprefix('chr')
    pos=int(pos);ref=ref.upper();alt=alt.upper()
    if pos < 1 or not ref or not alt or ref == alt or not set(ref+alt)<=set('ACGT'): return None
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]:ref,alt=ref[:-1],alt[:-1]
    while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:ref,alt,pos=ref[1:],alt[1:],pos+1
    return f'{chrom}:{pos}:{ref}:{alt}'

def summarize_record(record):
    if record is None:
        return {'status':'not_returned_by_gene_query','af':None,'frequency_source':None,'groupmax_af':None,'groupmax_population':None,'groupmax_an':None,'caution':'Not an assertion of absence from all gnomAD records; no zero frequency or rarity bonus assigned.'}
    result={'status':'matched_frequency_unusable','af':None,'frequency_source':None,'groupmax_af':None,'groupmax_population':None,'groupmax_an':None,'gnomad_variant_id':record['variant_id'],'flags':record.get('flags'), 'transcript_consequence':record.get('transcript_consequence'), 'datasets':{}}
    for typ in ['exome','genome','joint']:
        d=record.get(typ)
        if d is not None:
            summary=copy.deepcopy(d)
            summary['af_calculated']=count_af(d.get('ac'),d.get('an'))
            summary['pass_for_frequency']=passed(d)
            result['datasets'][typ]=summary
    usable=[(typ,record[typ]) for typ in ['exome','genome'] if passed(record.get(typ))]
    if passed(record.get('joint')):
        source,d='joint',record['joint']
    elif usable:
        source,d=max(usable,key=lambda pair:count_af(pair[1]['ac'],pair[1]['an']))
    else:return result
    result.update(status='matched_usable_frequency',af=count_af(d['ac'],d['an']),frequency_source=source)
    # Group maxima use only the selected callset, avoiding duplicated aggregate rows.
    groups=[(count_af(p.get('ac'),p.get('an')),p) for p in d.get('populations') or [] if p.get('id') in GROUPS]
    groups=[(af,p) for af,p in groups if af is not None]
    if groups:
        af,p=max(groups,key=lambda pair:(pair[0],pair[1]['an']))
        result.update(groupmax_af=af,groupmax_population=p['id'],groupmax_an=p['an'])
    return result

def match_candidates(candidates,gene):
    exact={}
    minimal={}
    for record in gene['variants']:
        if record.get('reference_genome')!='GRCh38':raise ValueError('Wrong reference assembly in response')
        key=f"{record['chrom']}:{record['pos']}:{record['ref']}:{record['alt']}"
        if key in exact:raise ValueError('Duplicate exact public allele')
        exact[key]=record
        normalized=normalize_variant(record['chrom'],record['pos'],record['ref'],record['alt'])
        if normalized is not None:minimal.setdefault(normalized,[]).append(record)
    results=[]
    for candidate in candidates:
        key=candidate['variant_id'];record=exact.get(key);kind='exact' if record else 'unmatched'
        if record is None:
            records=minimal.get(normalize_variant(candidate['chrom'],candidate['pos'],candidate['ref'],candidate['alt']),[])
            if len(records)==1:record=records[0];kind='unique_minimal_representation'
            elif len(records)>1:kind='ambiguous_minimal_representation'
        row={'variant_id':key,'gene':candidate['gene'],'protein_change':candidate['protein_change'],'consequence':candidate['consequence'],'genotype':candidate['genotype'],'legacy_population_af':candidate['population_af'],'legacy_population_source':candidate['population_af_source'],'match_type':kind,**summarize_record(record)}
        # Gene-level coverage is binned, not exact per-site AN or a callable absence test.
        row['nearby_binned_coverage']={}
        for typ,coverage in gene.get('coverage',{}).items():
            if coverage:
                nearby=min(coverage,key=lambda p:abs(p['pos']-candidate['pos']))
                row['nearby_binned_coverage'][typ]={**nearby,'distance_bp':abs(nearby['pos']-candidate['pos']),'interpretation':'Nearest displayed coverage bin; not exact site coverage or proof that an unreturned allele is absent.'}
        results.append(row)
    return results

def apply_frequencies(candidates,evidence,clinical=True,groupmax=False):
    lookup={(r['gene'],r['variant_id']):r for r in evidence}
    rows=copy.deepcopy(candidates)
    for r in rows:
        e=lookup[(r['gene'],r['variant_id'])]
        value=e['groupmax_af'] if groupmax else e['af']
        # Never silently fall back to ClinVar frequency when an independent audit is unknown.
        r['population_af']=value
        r['population_af_source']=f"gnomAD v4 {e['frequency_source']} {'groupmax AF sensitivity' if groupmax else 'AC/AN'}" if value is not None else None
        part=1 if value is not None and value<.001 else -3 if value is not None and value>.01 else 0
        r['score_components']['external_population_frequency']=part
        r['triage_score']=sum(r['score_components'].values())
        r['supplied_cis_partner_ids']=[]
        if not clinical:r['clinvar'],r['clinical_score']=[],0
    return rows

def csv_write(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})

def acquire(config=CONFIG):
    policy=read(config);cache=Path(policy['reference_cache']);cache.mkdir(parents=True,exist_ok=True)
    query=PUBLIC_QUERY
    (cache/'public_gene_query.graphql').write_text(query,encoding='utf-8')
    models=policy['public_gene_ids']
    if set(models)!=set(policy['genes']):raise ValueError('Gene ID list must match authorized scope')
    for symbol in policy['genes']:
        dest=cache/f'{symbol}.json';receipt=cache/f'{symbol}.source.json'
        if dest.exists() and receipt.exists():
            prior=read(receipt)
            if prior.get('response_sha256')==sha(dest) and prior.get('query_sha256')==hashlib.sha256(query.encode()).hexdigest():
                data=read(dest)
                if not data.get('errors') and (data.get('data') or {}).get('gene'):
                    print(f'Using verified public cache: {symbol}',flush=True);continue
        variables={'gene':models[symbol]}
        body=json.dumps({'query':query,'variables':variables}).encode()
        for attempt in range(3):
            try:
                request=urllib.request.Request(policy['endpoint'],data=body,headers={'Content-Type':'application/json','User-Agent':'MVA-local-research/0.3'})
                with urllib.request.urlopen(request,timeout=60) as response:raw=response.read();status=response.status
                data=json.loads(raw)
                if data.get('errors') or not (data.get('data') or {}).get('gene'):raise ValueError(str(data.get('errors')))
                if data['data']['gene']['gene_id']!=variables['gene']:raise ValueError('Gene response mismatch')
                dest.write_bytes(raw)
                write(receipt,{'endpoint':policy['endpoint'],'dataset':policy['dataset'],'reference_genome':'GRCh38','variables':variables,'query_sha256':hashlib.sha256(query.encode()).hexdigest(),'response_sha256':hashlib.sha256(raw).hexdigest(),'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'http_status':status,'public_variant_records':len(data['data']['gene']['variants'])})
                print(f"Cached {symbol}: {len(data['data']['gene']['variants']):,} public variants",flush=True)
                break
            except Exception:
                if attempt==2:raise
                time.sleep(2*(attempt+1))
        time.sleep(1)

def run(output,config=CONFIG):
    output=Path(output)
    if output.exists():raise ValueError('Choose a new output directory')
    policy=read(config);cache=Path(policy['reference_cache'])
    candidates=read(policy['source_candidates'])
    candidates=[r for r in candidates if r['gene'] in policy['genes']]
    phenotype=PhenotypeIndex(read(policy['phenotype_input']))
    evidence=[];inputs={str(config):sha(config),policy['source_candidates']:sha(policy['source_candidates']),policy['phenotype_input']:sha(policy['phenotype_input'])}
    for symbol in policy['genes']:
        response=cache/f'{symbol}.json';receipt=cache/f'{symbol}.source.json'
        metadata=read(receipt)
        if metadata['response_sha256']!=sha(response):raise ValueError('Cache hash mismatch')
        if metadata['query_sha256']!=hashlib.sha256(PUBLIC_QUERY.encode()).hexdigest():raise ValueError('Query hash mismatch')
        data=read(response)
        if data.get('errors'):raise ValueError('Partial API response')
        if data['data']['gene']['gene_id']!=policy['public_gene_ids'][symbol]:raise ValueError('Wrong gene in response')
        evidence.extend(match_candidates([r for r in candidates if r['gene']==symbol],data['data']['gene']))
        inputs[str(response)]=sha(response);inputs[str(receipt)]=sha(receipt)
    output.mkdir(parents=True)
    write(output/'policy.json',policy)
    write(output/'allele_evidence.json',evidence)
    flat=[]
    for r in evidence:
        row={k:r.get(k) for k in ['gene','variant_id','protein_change','consequence','genotype','match_type','status','legacy_population_af','af','frequency_source','groupmax_af','groupmax_population','groupmax_an']}
        for typ in ['exome','genome','joint']:
            d=r.get('datasets',{}).get(typ,{})
            for k in ['ac','an','homozygote_count','af_calculated','filters','flags','fafmax']:row[f'{typ}_{k}']=d.get(k)
        flat.append(row)
    csv_write(output/'allele_frequencies.csv',flat)
    variants={
        'legacy_full':copy.deepcopy(candidates),
        'legacy_remove_clinvar_and_embedded_frequency':remove_clinvar_evidence(copy.deepcopy(candidates)),
        'gnomad_full':apply_frequencies(candidates,evidence),
        'gnomad_remove_clinvar_classification':apply_frequencies(candidates,evidence,clinical=False),
        'gnomad_groupmax_sensitivity':apply_frequencies(candidates,evidence,groupmax=True),
        'gnomad_groupmax_remove_clinvar_classification':apply_frequencies(candidates,evidence,clinical=False,groupmax=True),
    }
    gene_rows=[];pair_rows=[];qc={};eligibility_rows=[]
    for scenario,rows in variants.items():
        hypotheses,quality=build_hypotheses(rows,phenotype)
        ranked=rank_genes(hypotheses,phenotype)
        pair_ranked=rank_genes([h for h in hypotheses if h['type']!='single_heterozygous_exploratory'],phenotype)
        gene_rows.extend({'audit_scenario':scenario,**r} for r in ranked)
        pair_rows.extend({'audit_scenario':scenario,**r} for r in pair_ranked)
        eligibility_rows.extend({'audit_scenario':scenario,'gene':r['gene'],'variant_id':r['variant_id'],'population_af':r['population_af'],'exclusion':r['hypothesis_exclusion_reason']} for r in rows)
        qc[scenario]=quality
    for name,rows in [('gene_ranks',gene_rows),('biallelic_ranks',pair_rows),('eligibility',eligibility_rows)]:
        write(output/(name+'.json'),rows);csv_write(output/(name+'.csv'),rows)
    summary={'status':'complete_selected_gene_audit','genes_requested':len(policy['genes']),'candidate_annotations':len(candidates),'match_status_counts':dict(Counter(r['status'] for r in evidence)),'match_type_counts':dict(Counter(r['match_type'] for r in evidence)),'bub1b_alleles':[r for r in flat if r['gene']=='BUB1B'],'bub1b_gene_ranks':[{k:r[k] for k in ['audit_scenario','rank','score','variant_ids']} for r in gene_rows if r['gene']=='BUB1B'],'bub1b_biallelic_ranks':[{k:r[k] for k in ['audit_scenario','rank','score','variant_ids']} for r in pair_rows if r['gene']=='BUB1B'],'scope':'Ranks are within the fixed selected gene set; not a newly annotated genome-wide rank.','missingness':'Unreturned or filtered records remain unknown. Minimal representation is attempted; repeat-aware left alignment is not implemented. Coverage bins do not establish callable absence.','clinical_interpretation':'Population frequency informs compatibility; rarity, zero observed homozygotes and a high rank do not establish pathogenicity or phase.'}
    write(output/'summary.json',summary);write(output/'ranking_qc.json',qc)
    write(output/'run_manifest.json',{'completed_at_utc':datetime.now(timezone.utc).isoformat(),'status':'complete','inputs_sha256':inputs,'code_sha256':{str(Path(__file__)):sha(__file__),'mva/genomewide.py':sha('mva/genomewide.py'),'mva/track1.py':sha('mva/track1.py')},'outputs_sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file()}})
    print(json.dumps({k:summary[k] for k in ['status','genes_requested','candidate_annotations','match_status_counts','bub1b_gene_ranks','bub1b_biallelic_ranks']},indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CONFIG)
    sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('download');r=sub.add_parser('run');r.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.command=='download':acquire(a.config)
    else:run(a.output,a.config)

if __name__=='__main__':main()
