"""Local site-specific read evidence from a duplicate-marked competitive alignment.

Requires pysam only for the BAM reader. Pure fragment/QC functions are testable
without patient inputs or alignment software. No network requests are made.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,statistics
from collections import Counter,defaultdict
from pathlib import Path

def exclusion_reasons(row,baseq=20,deduplicate=True,strict_k21=False):
    reasons=[]
    for field in ['unmapped','secondary','supplementary','qcfail']:
        if row[field]:reasons.append(field)
    if not row['proper_pair']:reasons.append('not_proper_pair')
    if row['mapq']<30:reasons.append('mapq_below_30')
    if row['base'] is None:reasons.append('deletion_or_refskip_at_site')
    elif row['baseq']<baseq:reasons.append(f'baseq_below_{baseq}')
    if row['distance_to_read_end'] is None or row['distance_to_read_end']<5:reasons.append('within_five_bases_of_read_end')
    if deduplicate and row['duplicate']:reasons.append('duplicate')
    if strict_k21 and not row['k21_pair']:reasons.append('not_k21_recruited')
    return reasons

def collapse_fragments(rows,baseq=20,deduplicate=True,strict_k21=False):
    groups=defaultdict(list);reasons=Counter()
    for row in rows:
        excluded=exclusion_reasons(row,baseq,deduplicate,strict_k21)
        reasons.update(excluded)
        if not excluded:groups[(row['lane'],row['read_name'])].append(row)
    fragments=[];discordant=0;overlap=0
    for records in groups.values():
        if len({r['base'] for r in records})!=1:discordant+=1;continue
        if len(records)>1:overlap+=1
        chosen=max(records,key=lambda r:(r['baseq'],r['mapq'],r['read1']))
        fragments.append({**chosen,'qualified_reads_in_fragment':len(records)})
    return fragments,{'excluded_read_reason_counts':dict(reasons),'discordant_fragments':discordant,'concordant_overlapping_fragments':overlap,'qualifying_read_alignments':sum(map(len,groups.values()))}

def aggregate(rows,ref,alt):
    counts=Counter(r['base'] for r in rows);n=sum(counts.values());denom=counts[ref]+counts[alt]
    by_lane={lane:dict(Counter(r['base'] for r in rows if r['lane']==lane)) for lane in sorted({r['lane'] for r in rows})}
    quality={}
    for allele in sorted(counts):
        group=[r for r in rows if r['base']==allele]
        quality[allele]={'n':len(group),'baseq_min':min(r['baseq'] for r in group),'baseq_median':statistics.median(r['baseq'] for r in group),
                         'mapq_min':min(r['mapq'] for r in group),'mapq_median':statistics.median(r['mapq'] for r in group),
                         'distance_to_read_end_min':min(r['distance_to_read_end'] for r in group),
                         'fragment_forward':sum(not r['fragment_reverse'] for r in group),'fragment_reverse':sum(r['fragment_reverse'] for r in group),
                         'representative_read_forward':sum(not r['reverse'] for r in group),'representative_read_reverse':sum(r['reverse'] for r in group),
                         'softclip_max':max(r['softclip_bases'] for r in group),'mismatch_edit_distance_max':max((r['edit_distance'] or 0) for r in group)}
    return {'fragments':n,'reference':counts[ref],'alternate':counts[alt],'other':n-denom,'alternate_fraction':counts[alt]/denom if denom else None,'by_lane':by_lane,'allele_quality':quality}

def csv_write(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()})

def analyze(bam_path,candidates,directory):
    import pysam
    directory=Path(directory);variants=json.loads(Path(candidates).read_text());recruited={}
    for file in directory.glob('L*_recruitment.tsv'):
        lane=file.name.split('_')[0]
        with file.open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                name=r['read_name'].removesuffix('/1').removesuffix('/2')
                recruited[(lane,name)]=bool(int(r['k21_r1']) or int(r['k21_r2']))
    bam=pysam.AlignmentFile(bam_path,'rb');all_rows=[];summaries=[];pileups=[];calls=defaultdict(set)
    for v in variants:
        chrom='chr'+v['chrom'].removeprefix('chr');pos=v['pos']-1;rows=[]
        for r in bam.fetch(chrom,pos,pos+1):
            pairs=r.get_aligned_pairs();positions={rp:qp for qp,rp in pairs if rp is not None};qp=positions.get(pos)
            lane=r.get_tag('RG') if r.has_tag('RG') else 'unknown';name=r.query_name
            if (lane,name) not in recruited:raise ValueError('Missing recruitment provenance for read')
            row={'variant_id':v['variant_id'],'protein_change':v['protein_change'],'lane':lane,'read_name':name,'read1':r.is_read1,
                 'unmapped':r.is_unmapped,'secondary':r.is_secondary,'supplementary':r.is_supplementary,'qcfail':r.is_qcfail,
                 'proper_pair':r.is_proper_pair,'duplicate':r.is_duplicate,'reverse':r.is_reverse,
                 'fragment_reverse':r.is_reverse if r.is_read1 else not r.is_reverse,
                 'mapq':r.mapping_quality,'base':r.query_sequence[qp] if qp is not None else None,
                 'baseq':r.query_qualities[qp] if qp is not None else None,
                 'distance_to_read_end':min(qp,r.query_length-1-qp) if qp is not None else None,
                 'original_read_cycle':(r.query_length-qp if r.is_reverse else qp+1) if qp is not None else None,
                 'cigar':r.cigarstring,'reference_start_0based':r.reference_start,'reference_end_exclusive':r.reference_end,
                 'template_length':r.template_length,'mate_reference_start_0based':r.next_reference_start,
                 'edit_distance':r.get_tag('NM') if r.has_tag('NM') else None,
                 'softclip_bases':sum(n for op,n in r.cigartuples or [] if op==4),'k21_pair':recruited[(lane,name)]}
            row['primary_exclusions']=exclusion_reasons(row)
            rows.append(row)
            if not row['primary_exclusions']:
                bases=[{'pos':rp+1,'base':r.query_sequence[q],'baseq':r.query_qualities[q]} for q,rp in pairs if q is not None and rp is not None and abs(rp-pos)<=45]
                pileups.append({**row,'bases':bases})
        all_rows.extend(rows)
        comparisons={}
        for label,options in [('before_duplicate_removal',{'deduplicate':False}),('primary',{}),('baseq30',{'baseq':30}),('k21_recruitment',{'strict_k21':True})]:
            fragments,qc=collapse_fragments(rows,**options);comparisons[label]={**aggregate(fragments,v['ref'],v['alt']),**qc}
            if label=='primary':
                for r in fragments:calls[(r['lane'],r['read_name'])].add(v['variant_id'])
        summaries.append({'variant_id':v['variant_id'],'protein_change':v['protein_change'],'ref':v['ref'],'alt':v['alt'],
                          'original_vcf':{k:v[k] for k in ['depth','reference_depth','alternate_depth','allele_fraction']},
                          'retrieved_alignments_over_site':len(rows),'comparisons':comparisons})
    bam.close()
    scans={file.stem:json.loads(file.read_text()) for file in directory.glob('L*_summary.json')}
    assert len(scans)==4 and all(r['complete_scan'] and r['gzip_streams_valid_to_end'] for r in scans.values())
    result={'status':'complete_targeted_read_review','variants':summaries,'lane_scans':scans,
            'total_read_pairs_scanned':sum(r['pairs_scanned'] for r in scans.values()),
            'recruited_pairs_k17':sum(r['pairs_recruited_k17'] for r in scans.values()),
            'recruited_pairs_k21':sum(r['pairs_recruited_k21'] for r in scans.values()),
            'primary_fragments_covering_both_sites':sum(len(s)==len(variants) for s in calls.values()),
            'scope':'Full FASTQ scan with targeted unique-kmer retrieval; whole-reference alignment of recruited pairs. Same sequencing experiment, not an independent specimen or full-genome re-call.',
            'phase':'Separate site support does not establish trans phase. Short fragments cannot ordinarily bridge this inter-site distance.'}
    (directory/'site_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (directory/'pileup_data.json').write_text(json.dumps(pileups,indent=2)+'\n')
    csv_write(directory/'read_observations.csv',all_rows)
    flat=[{'variant_id':v['variant_id'],'protein_change':v['protein_change'],'comparison':label,**{k:d[k] for k in ['fragments','reference','alternate','other','alternate_fraction','discordant_fragments','concordant_overlapping_fragments']}} for v in summaries for label,d in v['comparisons'].items()]
    csv_write(directory/'site_counts.csv',flat)
    print(json.dumps(result,indent=2))
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bam',type=Path,required=True);p.add_argument('--candidates',type=Path,required=True);p.add_argument('--directory',type=Path,required=True);a=p.parse_args();analyze(a.bam,a.candidates,a.directory)
if __name__=='__main__':main()
