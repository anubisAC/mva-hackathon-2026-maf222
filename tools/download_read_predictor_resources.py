"""Fetch only public reference/software files; no case data are sent."""
import concurrent.futures,hashlib,json,time,urllib.request
from pathlib import Path
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]
def get(name,url,folder,expected=None):
    dest=ROOT/folder/name;dest.parent.mkdir(parents=True,exist_ok=True)
    part=dest.with_name(dest.name+'.partial')
    if dest.exists() and dest.with_name(dest.name+'.source.json').exists():
        print('Cached '+name,flush=True);return
    for attempt in range(3):
        try:
            h=hashlib.sha256();md5=hashlib.md5();total=0;last=time.monotonic()
            with urllib.request.urlopen(url,timeout=120) as src,part.open('wb') as out:
                size=src.headers.get('Content-Length')
                while block:=src.read(4*1024*1024):
                    out.write(block);h.update(block);md5.update(block);total+=len(block)
                    if time.monotonic()-last>30:
                        print(f'{name}: {total/1e6:.0f} MB'+(f' / {int(size)/1e6:.0f} MB' if size else ''),flush=True);last=time.monotonic()
            if size and total!=int(size):raise ValueError('Truncated public download')
            if expected and h.hexdigest()!=expected:raise ValueError('Published SHA256 mismatch')
            part.replace(dest)
            dest.with_name(dest.name+'.source.json').write_text(json.dumps({'url':url,'sha256':h.hexdigest(),'md5':md5.hexdigest(),'bytes':total,'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'published_sha256_verified':bool(expected)},indent=2)+'\n')
            print(f'Complete {name}: {total:,} bytes',flush=True);return
        except Exception as e:
            print(f'{name} attempt {attempt+1}: {e}',flush=True)
            if attempt==2:raise
            time.sleep(3)

base='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/'
metadata=ROOT/'references/read_validation/pysam_pypi.json'
metadata.parent.mkdir(parents=True,exist_ok=True)
if not metadata.exists():
 with urllib.request.urlopen('https://pypi.org/pypi/pysam/0.24.1/json') as src:metadata.write_bytes(src.read())
pysam=json.loads(metadata.read_text())
wheel=next(r for r in pysam['urls'] if 'cp312' in r['filename'] and 'manylinux' in r['filename'] and 'x86_64' in r['filename'])
jobs=[
 ('revel_zenodo.json','https://zenodo.org/api/records/7072866','references/predictors',None),
 ('alphamissense_objects.json','https://storage.googleapis.com/storage/v1/b/dm_alphamissense/o','references/predictors',None),
 ('GRCh38_no_alt_plus_hs38d1.fna.gz',base+'GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz','references/read_validation',None),
 ('reference_md5checksums.txt',base+'md5checksums.txt','references/read_validation',None),
 ('reference_README.txt',base+'README_analysis_sets.txt','references/read_validation',None),
 ('GRC_exclusions.bed',base+'GCA_000001405.15_GRCh38_GRC_exclusions.bed','references/read_validation',None),
 ('minimap2-2.31_x64-linux.tar.bz2','https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31_x64-linux.tar.bz2','references/read_validation','300bc287f05eb890c6211fa7db043ce98320a401621fadd1cfdbeabd1a6e4ab5'),
 ('kseq.h','https://raw.githubusercontent.com/lh3/minimap2/v2.31/kseq.h','references/read_validation',None),
 (wheel['filename'],wheel['url'],'references/read_validation',wheel['digests']['sha256']),
 ('revel-v1.3_segments_chrom_15.zip','https://zenodo.org/api/records/7072866/files/revel-v1.3_segments_chrom_15.zip/content','references/predictors',None),
 ('AlphaMissense_hg38.tsv.gz','https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz','references/predictors',None),
 ('CADD_v1.7_GRCh38.tsv.gz.tbi','https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi','references/predictors',None),
]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    futures=[pool.submit(get,*job) for job in jobs]
    for f in concurrent.futures.as_completed(futures):f.result()
