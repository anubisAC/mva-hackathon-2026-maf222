"""Original-dataset predictor lookup. Public gene range requests; allele join local."""
import csv,gzip,hashlib,io,json,os,sys,time,urllib.request,zipfile
from pathlib import Path
import pysam

ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT)
REF=ROOT/'references/predictors';OUT=ROOT/'outputs/predictor_review_2026-10-20';OUT.mkdir(exist_ok=True)
variants=json.loads((ROOT/'outputs/strategy_b/context/candidates.json').read_text())
keys={(str(v['chrom']),v['pos'],v['ref'],v['alt']) for v in variants}
found={v['variant_id']:{} for v in variants}
def vid(key):return ':'.join(map(str,key))
def save(name,obj):
 (OUT/name).write_text(json.dumps(obj,indent=2)+'\n')
def extract_local():
 archive=REF/'revel-v1.3_segments_chrom_15.zip'
 receipt=json.loads((REF/'revel_zenodo.json').read_text())
 record=next(r for r in receipt['files'] if r['key']==archive.name)
 assert hashlib.md5(archive.read_bytes()).hexdigest()==record['checksum'].split(':')[1]
 with zipfile.ZipFile(archive) as z:
  for name in z.namelist():
   if name.endswith('.csv'):
    with z.open(name) as f:
     for r in csv.DictReader(io.TextIOWrapper(f)):
      pos=r.get('grch38_pos');key=(r['chr'],int(pos) if pos and pos!='.' else -1,r['ref'],r['alt'])
      if key in keys:found[vid(key)].setdefault('REVEL',[]).append(r)
 save('local_scores.json',found);print('REVEL',json.dumps(found),flush=True)
 with gzip.open(REF/'AlphaMissense_hg38.tsv.gz','rt') as f:
  header=None;count=0
  for line in f:
   if line.startswith('#CHROM'):header=line.lstrip('#').rstrip().split('\t');continue
   if line.startswith('#'):continue
   if header is None:raise ValueError('AlphaMissense column header not found')
   fields=line.rstrip().split('\t');chrom=fields[0].removeprefix('chr');key=(chrom,int(fields[1]),fields[2],fields[3]);count+=1
   if key in keys:found[vid(key)].setdefault('AlphaMissense',[]).append(dict(zip(header,fields)))
 save('local_scores.json',found);print('AlphaMissense records scanned',count,flush=True)

def tabix_gene(label,url,index_url=None,index_local=None):
 index=REF/(index_local or label+'.tbi')
 if not index.exists():
  with urllib.request.urlopen(index_url or url+'.tbi',timeout=120) as src:index.write_bytes(src.read())
 # The complete public BUB1B gene interval, not the patient variant coordinates.
 start,end=40161022,40221123
 t=pysam.TabixFile(url,index=str(index))
 header=list(t.header);chrom='15' if '15' in t.contigs else 'chr15'
 dest=REF/(label+'_BUB1B.tsv');n=0;matches=[]
 with dest.open('w') as out:
  out.write('\n'.join(header)+'\n')
  for line in t.fetch(chrom,start,end):
   out.write(line+'\n');n+=1;fields=line.split('\t')
   key=(fields[0].removeprefix('chr'),int(fields[1]),fields[2],fields[3]) if label=='CADD' else (fields[0].removeprefix('chr'),int(fields[1]),fields[3],fields[4])
   if key in keys:
    value={'record':line,'header':header};found[vid(key)][label]=value;matches.append(value)
 t.close()
 receipt={'url':url,'index_sha256':hashlib.sha256(index.read_bytes()).hexdigest(),'public_gene':'BUB1B','interval_0based':[chrom,start,end],'data_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'records':n,'retrieved_at_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'request_scope':'public full-gene reference; case matching local'}
 (dest.with_name(dest.name+'.source.json')).write_text(json.dumps(receipt,indent=2)+'\n')
 save('remote_scores.json',found);print(label,'public gene rows',n,'case allele matches',len(matches),flush=True)

if '--remote-only' not in sys.argv:extract_local()
if '--local-only' in sys.argv:raise SystemExit(0)
tabix_gene('CADD','https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz',index_local='CADD_v1.7_GRCh38.tsv.gz.tbi')
base='https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/'
for mode in ['masked','raw']:
 tabix_gene('SpliceAI_'+mode,base+'spliceai_scores.'+mode+'.snv.ensembl_mane_v1.4.grch38.vcf.gz')
save('remote_scores.json',found)
print(json.dumps(found,indent=2),flush=True)
