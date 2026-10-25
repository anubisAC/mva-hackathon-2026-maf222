"""Run on Linux after compiling tools/recruit_reads.c; no patient inputs used."""
import gzip,json,random,subprocess,tempfile,unittest
from pathlib import Path

BIN=Path('tools/read_validation_runtime/recruit_reads').resolve()

@unittest.skipUnless(BIN.exists() and __import__('os').name=='posix','compiled Linux recruiter required')
class ReadRecruitmentTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.p=Path(self.tmp.name)
  rng=random.Random(7);self.seq=''.join(rng.choice('ACGT') for _ in range(1800))
  (self.p/'genome.fa').write_text('>unique\n'+self.seq+'\n>duplicated\n'+self.seq[:100]+'\n')
  (self.p/'bait.fa').write_text('>bait\n'+self.seq[600:750]+'N'+self.seq[751:900]+'\n')
  subprocess.run([str(BIN),'unique',str(self.p/'bait.fa'),str(self.p/'genome.fa'),str(self.p/'kmers')],check=True,capture_output=True)
 def run_reads(self,one,two,expected=0,limit=None):
  for name,reads in [('r1',one),('r2',two)]:
   with gzip.open(self.p/(name+'.fq.gz'),'wt') as f:
    for ident,seq in reads:f.write('@'+ident+'\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
  args=[str(BIN),'scan',str(self.p/'kmers'),str(self.p/'r1.fq.gz'),str(self.p/'r2.fq.gz'),str(self.p/'out')]
  if limit:args.append(str(limit))
  result=subprocess.run(args,capture_output=True)
  self.assertEqual(result.returncode,expected,result.stderr)
  if expected==0:return json.loads((self.p/'out_summary.json').read_text())
 def test_reference_alternate_reverse_and_mate_only_are_retained(self):
  ref=self.seq[680:830];alt=ref[:70]+('A' if ref[70]!='A' else 'C')+ref[71:]
  rc=alt.translate(str.maketrans('ACGT','TGCA'))[::-1];no='A'*150
  a=[('ref/1',ref),('alt/1',alt),('reverse/1',rc),('mate/1',no),('negative/1',no)]
  b=[(ident.replace('/1','/2'),ref if ident.startswith('mate') else no) for ident,_ in a]
  r=self.run_reads(a,b)
  self.assertEqual((r['pairs_scanned'],r['pairs_recruited_k17'],r['pairs_recruited_k21']),(5,4,4))
  self.assertTrue(r['complete_scan']);self.assertTrue(r['gzip_streams_valid_to_end'])
 def test_mate_name_mismatch_is_an_error(self):self.run_reads([('one',self.seq[:150])],[('two',self.seq[:150])],expected=2)
 def test_mate_record_count_mismatch_is_an_error(self):self.run_reads([('one',self.seq[:150])],[],expected=2)
 def test_limited_scan_never_claims_complete(self):
  reads=[('one',self.seq[600:750]),('two',self.seq[700:850])]
  r=self.run_reads(reads,reads,limit=1);self.assertFalse(r['complete_scan']);self.assertFalse(r['gzip_streams_valid_to_end'])
 def test_nonunique_reference_kmers_are_removed(self):
  (self.p/'bait_repeat.fa').write_text('>repeat\n'+self.seq[:100]+'\n')
  r=subprocess.run([str(BIN),'unique',str(self.p/'bait_repeat.fa'),str(self.p/'genome.fa'),str(self.p/'repeat_kmers')],capture_output=True)
  self.assertEqual(r.returncode,2);self.assertEqual((self.p/'repeat_kmers').read_text(),'')

if __name__=='__main__':unittest.main()
