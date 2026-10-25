import unittest
from mva.read_support import aggregate,collapse_fragments,exclusion_reasons

def row(**changes):
 r={'lane':'L1','read_name':'fragment','read1':True,'unmapped':False,'secondary':False,'supplementary':False,'qcfail':False,
    'proper_pair':True,'duplicate':False,'mapq':60,'base':'T','baseq':35,'distance_to_read_end':20,'k21_pair':True,
    'reverse':False,'fragment_reverse':False,'softclip_bases':0,'edit_distance':0}
 r.update(changes);return r

class ReadSupportTests(unittest.TestCase):
 def test_overlapping_mates_count_once(self):
  fragments,qc=collapse_fragments([row(),row(read1=False,baseq=30)])
  self.assertEqual(len(fragments),1);self.assertEqual(qc['concordant_overlapping_fragments'],1)
 def test_discordant_overlapping_mates_are_excluded(self):
  f,qc=collapse_fragments([row(),row(read1=False,base='G')])
  self.assertFalse(f);self.assertEqual(qc['discordant_fragments'],1)
 def test_same_name_in_different_lanes_is_not_merged(self):
  f,_=collapse_fragments([row(),row(lane='L2')]);self.assertEqual(len(f),2)
 def test_duplicate_threshold_and_kmer_sensitivities_are_separate(self):
  self.assertIn('duplicate',exclusion_reasons(row(duplicate=True)))
  self.assertFalse(exclusion_reasons(row(duplicate=True),deduplicate=False))
  self.assertIn('baseq_below_30',exclusion_reasons(row(baseq=25),baseq=30))
  self.assertIn('not_k21_recruited',exclusion_reasons(row(k21_pair=False),strict_k21=True))
 def test_deletion_low_mapq_and_read_edge_do_not_support_an_allele(self):
  for r in [row(base=None,baseq=None,distance_to_read_end=None),row(mapq=29),row(distance_to_read_end=4)]:
   self.assertTrue(exclusion_reasons(r))
 def test_read_quality_filter_precedes_mate_consensus(self):
  f,q=collapse_fragments([row(),row(read1=False,base='G',baseq=5)])
  self.assertEqual(len(f),1);self.assertEqual(f[0]['base'],'T');self.assertEqual(q['discordant_fragments'],0)
 def test_allele_fraction_uses_ref_plus_alt_and_reports_other(self):
  result=aggregate([row(),row(base='G'),row(base='A')],'T','G')
  self.assertEqual(result['alternate_fraction'],.5);self.assertEqual(result['other'],1)

if __name__=='__main__':unittest.main()
