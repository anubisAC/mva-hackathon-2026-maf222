"""Summarize cached predictor records with allele, protein and transcript checks.

No requests are made by this module. Missense-only models are not applicable to
stop-gained variants. Predictions remain separate from clinical classification.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path

def base_transcript(value):return value.split('.')[0]
def parse_spliceai(line,variant):
    fields=line.rstrip().split('\t')
    if len(fields)<8:raise ValueError('Invalid VCF record')
    key=(fields[0].removeprefix('chr'),int(fields[1]),fields[3],fields[4])
    if key!=(variant['chrom'].removeprefix('chr'),variant['pos'],variant['ref'],variant['alt']):return []
    info=dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
    result=[]
    for value in info.get('SpliceAI','').split(','):
        x=value.split('|')
        if len(x)!=10 or x[0]!=variant['alt'] or x[1]!=variant['gene']:continue
        scores=dict(zip(['acceptor_gain','acceptor_loss','donor_gain','donor_loss'],map(float,x[2:6])))
        if any(not 0<=s<=1 for s in scores.values()):raise ValueError('Invalid delta score')
        result.append({'gene':x[1],'scores':scores,'delta_positions':dict(zip(scores,map(int,x[6:10]))),'max_delta':max(scores.values())})
    return result

def select_missense(records,variant,model):
    if variant['consequence']!='missense_variant':return {'status':'not_applicable_missense_only','score':None}
    selected=[]
    for r in records:
        if model=='REVEL':
            key=(r['chr'],int(r['grch38_pos']),r['ref'],r['alt'])
            protein=variant['protein_change'].removeprefix('p.')
            same_protein=(r['aaref']==protein[0] and r['aaalt']==protein[-1])
            transcripts=r['Ensembl_transcriptid'].split(';');score=float(r['REVEL'])
        else:
            if r['genome']!='hg38':continue
            key=(r['CHROM'].removeprefix('chr'),int(r['POS']),r['REF'],r['ALT'])
            same_protein=r['protein_variant']==variant['protein_change'].removeprefix('p.')
            transcripts=[r['transcript_id']];score=float(r['am_pathogenicity'])
        if key!=(variant['chrom'],variant['pos'],variant['ref'],variant['alt']):continue
        if same_protein and base_transcript(variant['transcript_id']) in {base_transcript(t) for t in transcripts}:
            if not 0<=score<=1:raise ValueError('Invalid missense score')
            selected.append({'score':score,'source_record':r})
    if not selected:return {'status':'no_verified_matching_record','score':None}
    scores={r['score'] for r in selected}
    if len(scores)!=1:raise ValueError('Conflicting matched predictor scores')
    return {'status':'matched_allele_amino_acids_and_transcript_base_id','score':selected[0]['score'],'records':selected}

def run(candidates,cache,output):
    candidates=Path(candidates);cache=Path(cache);output=Path(output)
    variants=json.loads(candidates.read_text());local=json.loads((output/'local_scores.json').read_text())
    summaries=[];inputs=[candidates,output/'local_scores.json']
    sources=['CADD_BUB1B.tsv','SpliceAI_masked_BUB1B.tsv','SpliceAI_raw_BUB1B.tsv']
    records={name:(cache/name).read_text().splitlines() for name in sources};inputs.extend(cache/n for n in sources)
    if not any('CADD GRCh38-v1.7' in line for line in records['CADD_BUB1B.tsv'][:10]):raise ValueError('CADD build/version mismatch')
    for mode in ['masked','raw']:
        header=[line for line in records[f'SpliceAI_{mode}_BUB1B.tsv'] if line.startswith('#')]
        if '##reference=GRCh38/hg38' not in header or not any('SpliceAIv1.3.1' in line for line in header):raise ValueError('SpliceAI build/version mismatch')
    for v in variants:
        row={'variant_id':v['variant_id'],'protein_change':v['protein_change'],'transcript_id':v['transcript_id'],'consequence':v['consequence'],
             'REVEL':select_missense(local[v['variant_id']].get('REVEL',[]),v,'REVEL'),
             'AlphaMissense':select_missense(local[v['variant_id']].get('AlphaMissense',[]),v,'AlphaMissense')}
        cadd=[]
        for line in records['CADD_BUB1B.tsv']:
            if line.startswith('#'):continue
            f=line.split('\t')
            if (f[0],int(f[1]),f[2],f[3])==(v['chrom'],v['pos'],v['ref'],v['alt']):cadd.append({'raw':float(f[4]),'phred':float(f[5]),'record':line})
        if len(cadd)!=1:raise ValueError('Expected one exact CADD allele record')
        row['CADD']={'version':'1.7','assembly':'GRCh38',**cadd[0]}
        for mode in ['masked','raw']:
            matches=[]
            for line in records[f'SpliceAI_{mode}_BUB1B.tsv']:
                if not line.startswith('#'):matches.extend(parse_spliceai(line,v))
            unique={json.dumps(r,sort_keys=True):r for r in matches}
            if len(unique)!=1:raise ValueError('Missing or discordant BUB1B SpliceAI annotation')
            row['SpliceAI_'+mode]={'annotation':'Ensembl MANE v1.4','header_tool_version':'SpliceAIv1.3.1','distance_window_bp':50,'exact_source_matches':len(matches),'distinct_annotations':len(unique),**next(iter(unique.values()))}
        summaries.append(row)
    result={'status':'complete_cached_predictor_review','assembly':'GRCh38','variants':summaries,
            'ranking_policy':'Supporting evidence only; frozen gene/pair scores and subjective EPCR are unchanged. Correlated predictors are not counted as independent votes.',
            'transcript_version_note':'AlphaMissense uses ENST00000287598.10; local model is .11. Exact genomic allele, amino-acid change and transcript base ID match; transcript-version identity is not claimed.',
            'interpretation':'Model labels and scores are not a clinical pathogenicity classification, measured function, binding prediction or drug response. Low SpliceAI in the 50-bp window does not exclude all splice effects.'}
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    flat=[]
    for r in summaries:
        am=r['AlphaMissense'];flat.append({'variant_id':r['variant_id'],'protein_change':r['protein_change'],'REVEL':r['REVEL']['score'],'REVEL_status':r['REVEL']['status'],
          'CADD_raw':r['CADD']['raw'],'CADD_PHRED':r['CADD']['phred'],'AlphaMissense':am['score'],'AlphaMissense_label':am['records'][0]['source_record']['am_class'] if 'records' in am else 'not_applicable',
          'SpliceAI_masked_max':r['SpliceAI_masked']['max_delta'],'SpliceAI_raw_max':r['SpliceAI_raw']['max_delta'],'SpliceAI_window_bp':50})
    with (output/'scores.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    manifest={'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'analysis_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'outputs_sha256':{n:hashlib.sha256((output/n).read_bytes()).hexdigest() for n in ['summary.json','scores.csv']}}
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(flat,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--candidates',type=Path,default=Path('outputs/strategy_b/context/candidates.json'));p.add_argument('--cache',type=Path,default=Path('references/predictors'));p.add_argument('--output',type=Path,default=Path('outputs/predictor_review_2026-10-20'));a=p.parse_args();run(a.candidates,a.cache,a.output)
