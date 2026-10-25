"""Align fully retrieved pairs, mark duplicates and summarize sites locally."""
import hashlib,json,os,subprocess,time
from pathlib import Path
import pysam

ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT)
import argparse
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
p.add_argument('--mapper',type=Path,required=True);p.add_argument('--candidates',type=Path,required=True)
p.add_argument('--index',type=Path,help='Optional faster local copy of the same complete index')
a=p.parse_args();OUT=a.output.resolve();REF=a.reference.resolve();mapper=a.mapper.resolve()
lanes=['L001','L002','L003','L004'];summaries=[OUT/(lane+'_summary.json') for lane in lanes]
index=str(a.index.resolve()) if a.index else str(REF)+'.sr.mmi'
assert all(json.loads(p.read_text())['complete_scan'] for p in summaries)
commands=[]
for lane in lanes:
 sam=OUT/(lane+'.sam');bam=OUT/(lane+'.bam')
 args=[str(mapper),'-ax','sr','--secondary=yes','-N','10','-t','4','-R',f'@RG\\tID:{lane}\\tSM:PROBAND01\\tLB:WGS_shared',index,str(OUT/(lane+'_R1.fastq')),str(OUT/(lane+'_R2.fastq'))]
 with sam.open('w') as out,(OUT/(lane+'_alignment.log')).open('w') as err:subprocess.run(args,stdout=out,stderr=err,check=True)
 with pysam.AlignmentFile(sam,'r') as src,pysam.AlignmentFile(bam,'wb',template=src) as dst:
  for record in src:dst.write(record)
 commands.append(args);print('Aligned '+lane,flush=True)
def run(name,*args):
 print('Running samtools '+name,flush=True);commands.append(['samtools',name,*map(str,args)]);return getattr(pysam,name)(*map(str,args))
run('merge','-f','-o',OUT/'merged.bam',*[OUT/(lane+'.bam') for lane in lanes])
run('sort','-n','-o',OUT/'name_sorted.bam',OUT/'merged.bam')
run('fixmate','-m',OUT/'name_sorted.bam',OUT/'fixmate.bam')
run('sort','-o',OUT/'coordinate_sorted.bam',OUT/'fixmate.bam')
run('markdup','-s','-f',OUT/'duplicate_metrics.txt',OUT/'coordinate_sorted.bam',OUT/'site_reads.bam')
run('index',OUT/'site_reads.bam')
run('quickcheck',OUT/'site_reads.bam')
(OUT/'alignment_commands.json').write_text(json.dumps({'minimap2':'2.31-r1302','pysam':pysam.__version__,'samtools':pysam.__samtools_version__,'commands':commands,'duplicate_policy':'coordinate duplicates across all lanes, treated as one library; no UMI evidence'},indent=2)+'\n')
import sys
sys.path.insert(0,str(ROOT))
from mva.read_support import analyze
analyze(OUT/'site_reads.bam',a.candidates,OUT)
print('Competitive mapping, duplicate marking and site review complete.',flush=True)
