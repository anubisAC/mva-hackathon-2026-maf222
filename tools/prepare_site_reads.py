"""WSL-only preparation; all inputs and outputs remain in the shared workspace."""
import concurrent.futures,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
import pysam

ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT)
import argparse,tempfile,shutil,re
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference',type=Path,required=True);p.add_argument('--candidates',type=Path,required=True)
p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
p.add_argument('--mapper',type=Path,required=True);p.add_argument('--scanner',type=Path,default=Path('tools/read_validation_runtime/recruit_reads'))
a=p.parse_args()
REF=a.reference.resolve()
OUT=a.output.resolve()
if OUT.exists():raise SystemExit('Choose a fresh output directory')
OUT.mkdir()
pysam.faidx(str(REF))
fa=pysam.FastaFile(str(REF))
variants=json.loads(a.candidates.read_text())
windows=[]
with (OUT/'masked_baits.fa').open('w') as f:
 for v in variants:
  chrom='chr'+v['chrom'].removeprefix('chr');start=max(0,v['pos']-1-750);end=v['pos']+750
  seq=fa.fetch(chrom,start,end).upper();assert seq[v['pos']-1-start]==v['ref']
  masked=list(seq)
  for other in variants:
   if other['chrom']==v['chrom'] and start<=other['pos']-1<end:masked[other['pos']-1-start]='N'
  f.write(f'>{chrom}:{start+1}-{end}\n'+''.join(masked)+'\n')
  windows.append({'variant_id':v['variant_id'],'chrom':chrom,'start_0based':start,'end_exclusive':end,'reference':seq})
(OUT/'reference_windows.json').write_text(json.dumps(windows,indent=2)+'\n')
policy={'reference':str(REF),'reference_basis':'NCBI GRCh38 no-alt plus hs38d1; same assembly, not the exact original maskedGRC v2 file',
 'bait_flank_bases':750,'target_base_in_baits':'masked N; no allele is required for recruitment',
 'recruitment':'any genome-unique canonical 17-mer in either mate; retain both mates; record stricter 21-mer membership',
 'scope':'Every record in all four paired FASTQ lanes; targeted retrieval followed by whole-reference competitive mapping; not all-read alignment',
 'primary_quality':{'mapq_min':30,'baseq_min':20,'read_end_distance_min':5,'proper_pair_required':True,'duplicates':'removed by samtools markdup','secondary_supplementary_qcfail':'excluded','overlapping_mates':'collapse to one fragment; discordant base calls excluded'},
 'sensitivity':'repeat counts at BQ30 and with k21-supported pairs; report raw and nonduplicate support',
 'phase':'No phase inference across the two sites unless a valid same-fragment observation spans both.',
 'started_at_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
(OUT/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
scanner=str(a.scanner.resolve())
mapper=str(a.mapper.resolve())
def run(args,log):
 with (OUT/log).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 a=pool.submit(run,[scanner,'unique',str(OUT/'masked_baits.fa'),str(REF),str(OUT/'unique_kmers.tsv')],'unique_kmers.log')
 index_tmp=tempfile.TemporaryDirectory(prefix='mva-index-')
 b=pool.submit(run,[mapper,'-x','sr','-t','4','-d',str(Path(index_tmp.name)/'reference.mmi'),str(REF)],'reference_index.log')
 a.result();print((OUT/'unique_kmers.log').read_text(),flush=True)
 # Check every 150-bp read start overlapping each site, for both target alleles.
 sets={17:set(),21:set()}
 for line in (OUT/'unique_kmers.tsv').read_text().splitlines():k,val=map(int,line.split());sets[k].add(val)
 def kmers(seq,k):
  for i in range(len(seq)-k+1):
   sub=seq[i:i+k]
   if set(sub)<=set('ACGT'):
    f=0;r=0
    for j,c in enumerate(sub):v='ACGT'.index(c);f=(f<<2)|v;r|=(3-v)<<(2*j)
    yield min(f,r)
 verification=[]
 for v,w in zip(variants,windows):
  for alt in [v['ref'],v['alt']]:
   seq=list(w['reference']);center=v['pos']-1-w['start_0based'];seq[center]=alt;seq=''.join(seq)
   for k in [17,21]:
    hits=sum(bool(set(kmers(seq[start:start+150],k))&sets[k]) for start in range(center-149,center+1))
    verification.append({'variant_id':v['variant_id'],'allele':alt,'k':k,'ideal_reads_checked':150,'recruited':hits})
    assert hits==150,'Ideal allele capture failed'
 (OUT/'ideal_read_recruitment_check.json').write_text(json.dumps(verification,indent=2)+'\n')
 print('All ideal 150-base REF and ALT reads captured equally by k17 and k21.',flush=True)
 # Scan lanes while the reference index finishes.
 def scan(r1):
  lane=re.search(r'_(L[0-9]{3})_R1_',r1.name).group(1);r2=Path(str(r1).replace('_R1_','_R2_'))
  prefix=OUT/lane
  run([scanner,'scan',str(OUT/'unique_kmers.tsv'),str(r1),str(r2),str(prefix)],lane+'_scan.log')
  print(lane+' '+Path(str(prefix)+'_summary.json').read_text().strip(),flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as scans:
  jobs=[scans.submit(scan,p) for p in sorted(a.data.glob('*_R1_*.fastq.gz'))]
  for job in jobs:job.result()
 b.result()
 shutil.copyfile(Path(index_tmp.name)/'reference.mmi',str(REF)+'.sr.mmi')
 index_tmp.cleanup()
 print('All read scans and reference indexing complete.',flush=True)
