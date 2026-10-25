import unittest
from mva.predictor_review import parse_spliceai,select_missense

V={'chrom':'1','pos':101,'ref':'T','alt':'G','gene':'TEST','protein_change':'p.N2K','transcript_id':'ENST_TEST.11','consequence':'missense_variant'}
AM={'CHROM':'chr1','POS':'101','REF':'T','ALT':'G','genome':'hg38','transcript_id':'ENST_TEST.10','protein_variant':'N2K','am_pathogenicity':'.92','am_class':'likely_pathogenic'}

class PredictorReviewTests(unittest.TestCase):
 def test_correct_alpha_allele_protein_and_transcript_base_match(self):
  self.assertEqual(select_missense([AM],V,'AlphaMissense')['score'],.92)
 def test_wrong_allele_transcript_protein_or_build_does_not_match(self):
  for field,value in [('ALT','C'),('POS','102'),('transcript_id','OTHER'),('protein_variant','N3K'),('genome','hg19')]:
   self.assertIsNone(select_missense([{**AM,field:value}],V,'AlphaMissense')['score'])
 def test_stop_gained_is_not_applicable_not_a_zero_missense_score(self):
  r=select_missense([AM],{**V,'consequence':'stop_gained'},'AlphaMissense')
  self.assertEqual(r['status'],'not_applicable_missense_only');self.assertIsNone(r['score'])
 def test_conflicting_scores_are_not_silently_maximized(self):
  with self.assertRaises(ValueError):select_missense([AM,{**AM,'am_pathogenicity':'.1'}],V,'AlphaMissense')
 def test_spliceai_requires_matching_allele_and_gene_and_keeps_positions(self):
  line='1\t101\t.\tT\tG\t.\t.\tSpliceAI=G|OTHER|.99|0|0|0|1|2|3|4,G|TEST|.02|.01|0|0|-48|20|-15|-38'
  r=parse_spliceai(line,V);self.assertEqual(len(r),1);self.assertEqual(r[0]['max_delta'],.02);self.assertEqual(r[0]['delta_positions']['acceptor_gain'],-48)
  self.assertEqual(parse_spliceai(line,{**V,'alt':'C'}),[])
 def test_invalid_splice_score_is_rejected(self):
  with self.assertRaises(ValueError):parse_spliceai('1\t101\t.\tT\tG\t.\t.\tSpliceAI=G|TEST|1.2|0|0|0|1|2|3|4',V)

if __name__=='__main__':unittest.main()
