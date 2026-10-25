import itertools
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from mva.bub1b_loss import score, empirical_rank, collapse_gene_scores, consistency, adjust, run
from mva.drug_signatures import file_hash
from mva.enrichment import hypergeometric_tail


class LossComparisonTests(unittest.TestCase):
    def test_mask_updates_lists_module_and_background_together(self):
        # A cycle-only hit disappears from all three denominators/numerators.
        pair = {"up": 0b0011, "down": 0b1100}
        full = score(pair, 0b0101, 0b111111, minimum=1)
        reduced = score(pair, 0b0101, 0b111111, excluded=0b0001, minimum=1)
        self.assertEqual((full['up_hits'], full['down_hits']), (1, 1))
        self.assertEqual((reduced['background_genes'], reduced['module_genes']), (5, 1))
        self.assertEqual((reduced['up_list_size'], reduced['down_list_size']), (1, 2))
        self.assertEqual(reduced['score_pp'], -50)
        self.assertEqual(reduced['up_hits'], 0)

    def test_target_removed_and_empty_module_unassessable(self):
        result = score({"up": 1, "down": 2}, 1, 15, target=1, minimum=1)
        self.assertEqual(result['module_genes'], 0)
        self.assertFalse(result['eligible'])
        self.assertIsNone(result['score_pp'])

    def test_empirical_ties_and_small_reference_are_not_false_significance(self):
        self.assertEqual(empirical_rank(0, [0]*30)['p'], 1)
        self.assertEqual(empirical_rank(0, [0]*30)['percentile_midrank'], .5)
        self.assertIsNone(empirical_rank(10, [0]*19)['p'])
        self.assertAlmostEqual(empirical_rank(10, list(range(30)))['p'], 24/31)
        self.assertAlmostEqual(empirical_rank(50, list(range(30)))['p'], 2/31)

    def test_duplicate_profiles_do_not_reweight_genes_or_cells(self):
        values = [('A', 1)]*100 + [('B', -1)]
        self.assertEqual(collapse_gene_scores(values), {'A': 1, 'B': -1})
        rows = [{'cell':'X','score_pp':1,'all_q':.01}, {'cell':'Y','score_pp':-1,'all_q':.9}]
        self.assertEqual(consistency(rows)['median_score_pp'], 0)
        self.assertFalse(consistency(rows)['leave_one_out_direction_consistent'])
        with self.assertRaises(ValueError): consistency(rows + [rows[0]])

    def test_fdr_keeps_ineligible_missing_and_covers_all_views(self):
        rows = [{'p': .01, 'view': 'full'}, {'p': .04, 'view': 'masked'}, {'p': None}]
        self.assertEqual(adjust(rows), 2)
        self.assertEqual([r['q'] for r in rows], [.02, .04, None])

    def test_hypergeometric_matches_exhaustive_enumeration(self):
        for total in range(2, 9):
            for module_size in range(1, total):
                module = set(range(module_size))
                for list_size in range(1, total):
                    overlaps = [len(module.intersection(x)) for x in itertools.combinations(range(total), list_size)]
                    for hits in range(min(module_size, list_size)+1):
                        expected = sum(x >= hits for x in overlaps)/len(overlaps)
                        self.assertAlmostEqual(hypergeometric_tail(hits,total,module_size,list_size), expected, places=11)

    def test_end_to_end_uses_only_same_source_cell_time_and_retains_two_views(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            policy = {'gene':'BUB1B','reactome_release':97,'cycle_mask_root':'CYCLE',
                      'minimum_module_genes_in_background':1,'minimum_benchmark_genes':1,
                      'views':['full','cell_cycle_removed'],
                      'modules':[{'id':'UPR','domain':'protein_handling','label':'UPR'}]}
            (root/'policy.json').write_text(json.dumps(policy))
            pathways = {'CYCLE':{'id':'CYCLE','ancestor_ids':[],'gene_symbols':['BUB1B','CYC'],'name':'Cycle','source_url':'x'},
                        'UPR':{'id':'UPR','ancestor_ids':[],'gene_symbols':['HSPA5','CYC'],'name':'UPR','source_url':'x'}}
            (root/'reference.json').write_text(json.dumps({'release':97,'pathways':pathways}))
            def pair(term, up, down):
                return f'{term} up\t\t'+'\t'.join(up)+'\n'+f'{term} down\t\t'+'\t'.join(down)+'\n'
            data = pair('BUB1B::A::::96H::::KO::Cas9::B1::g1',['HSPA5','CYC'],['Z'])
            data += pair('OTHER::A::::96H::::KO::Cas9::B1::g2',['Z'],['HSPA5','CYC'])
            data += pair('WRONGTIME::A::::24H::::KO::Cas9::B1::g3',['UNIQUE_BAD'],['Z'])
            data += pair('WRONGCELL::B::::96H::::KO::Cas9::B1::g4',['UNIQUE_BAD'],['Z'])
            kd = pair('BUB1B::KOLF2.1J::hiPSCs::144H::::KD::CRISPRi',['HSPA5'],['Z'])
            kd += pair('OTHER::KOLF2.1J::hiPSCs::144H::::KD::CRISPRi',['Z'],['HSPA5'])
            for filename, content in [('lincs-l1000-xpr.gmt.gz',data),('cm4ai.gmt.gz',kd)]:
                path=root/filename
                with gzip.open(path,'wt') as stream: stream.write(content)
                path.with_name(filename+'.receipt.json').write_text(json.dumps({'url':'https://example.org/public','sha256':file_hash(path)}))
            result = run(root/'policy.json',root/'reference.json',root,root/'result')
            self.assertEqual(result['query_contexts'],2)
            backgrounds=json.loads((root/'result/observed_backgrounds.json').read_text())
            self.assertTrue(all('UNIQUE_BAD' not in b['genes'] for b in backgrounds))
            rows=json.loads((root/'result/cell_scores.json').read_text())
            self.assertEqual(len(rows),4)
            self.assertTrue(all(r['all_n_reference_genes']==1 for r in rows))
            masked=json.loads((root/'result/query_contexts.json').read_text())
            self.assertTrue(all('CYC' not in r['up_genes'] for r in masked if r['view']=='cell_cycle_removed'))
            with self.assertRaises(FileExistsError): run(root/'policy.json',root/'reference.json',root,root/'result')


if __name__ == '__main__':
    unittest.main()
