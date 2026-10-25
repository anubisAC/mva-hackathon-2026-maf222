import contextlib,importlib.util,io,json,tempfile,unittest
from pathlib import Path

@unittest.skipUnless(importlib.util.find_spec('pysam'),'optional pysam dependency required')
class ReadSupportBamTests(unittest.TestCase):
 def test_cigar_duplicate_and_fragment_counting_through_real_bam(self):
  import pysam
  from mva.read_support import analyze
  with tempfile.TemporaryDirectory() as temp:
   p=Path(temp);header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}],
                      'RG':[{'ID':f'L00{i}','SM':'synthetic'} for i in range(1,5)]}
   variants=[{'variant_id':'1:101:T:G','chrom':'1','pos':101,'ref':'T','alt':'G','gene':'TEST','protein_change':'p.X',
              'depth':10,'reference_depth':5,'alternate_depth':5,'allele_fraction':.5}]
   (p/'candidates.json').write_text(json.dumps(variants))
   records=[]
   for name,base,flag,cigar in [('ref','T',99,'100M'),('alt','G',99,'100M'),('dup','G',99|1024,'100M'),('deletion','T',99,'50M1D50M')]:
    r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='A'*50+base+'A'*49;r.flag=flag;r.reference_id=0;r.reference_start=50;r.mapping_quality=60;r.cigarstring=cigar;r.next_reference_id=0;r.next_reference_start=250;r.template_length=300;r.query_qualities=pysam.qualitystring_to_array('D'*100);r.set_tag('RG','L001');r.set_tag('NM',0);records.append(r)
   with pysam.AlignmentFile(p/'reads.bam','wb',header=header) as bam:
    for r in records:bam.write(r)
   pysam.index(str(p/'reads.bam'))
   for lane in ['L001','L002','L003','L004']:
    names='read_name\tk17_r1\tk17_r2\tk21_r1\tk21_r2\n'
    if lane=='L001':names+=''.join(f'{r.query_name}\t1\t0\t1\t0\n' for r in records)
    (p/(lane+'_recruitment.tsv')).write_text(names)
    (p/(lane+'_summary.json')).write_text(json.dumps({'pairs_scanned':100,'pairs_recruited_k17':len(records) if lane=='L001' else 0,'pairs_recruited_k21':len(records) if lane=='L001' else 0,'complete_scan':True,'gzip_streams_valid_to_end':True}))
   with contextlib.redirect_stdout(io.StringIO()):result=analyze(p/'reads.bam',p/'candidates.json',p)
   v=result['variants'][0];primary=v['comparisons']['primary'];raw=v['comparisons']['before_duplicate_removal']
   self.assertEqual((primary['reference'],primary['alternate'],primary['fragments']),(1,1,2))
   self.assertEqual(raw['alternate'],2)
   self.assertEqual(primary['excluded_read_reason_counts']['deletion_or_refskip_at_site'],1)

if __name__=='__main__':unittest.main()
