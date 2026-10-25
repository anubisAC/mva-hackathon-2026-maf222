"""Population-data edge cases that can change rare-disease candidate ranking."""
import copy
import unittest

from mva.gnomad_review import (apply_frequencies, count_af, match_candidates,
                               normalize_variant, summarize_record)


def callset(ac=1, an=10000, filters=None, populations=None):
    return {'ac': ac, 'an': an, 'filters': [] if filters is None else filters,
            'homozygote_count': 0, 'populations': populations or [], 'fafmax': None}


def record(**changes):
    value = {'variant_id': '1-101-A-C', 'reference_genome': 'GRCh38', 'chrom': '1',
             'pos': 101, 'ref': 'A', 'alt': 'C', 'exome': callset(), 'genome': None,
             'joint': callset(an=20000), 'flags': [], 'transcript_consequence': None}
    value.update(changes)
    return value


def candidate(**changes):
    value = {'variant_id': '1:101:A:C', 'chrom': '1', 'pos': 101, 'ref': 'A', 'alt': 'C',
             'gene': 'TEST', 'protein_change': 'p.X', 'consequence': 'missense_variant',
             'genotype': '0/1', 'population_af': .2, 'population_af_source': 'ClinVar',
             'score_components': {'predicted_consequence': 4, 'quality': 1,
                                  'external_population_frequency': -3},
             'triage_score': 2, 'clinical_score': 3, 'clinvar': [{'classification': 'P'}],
             'supplied_cis_partner_ids': []}
    value.update(changes)
    return value


class GnomadReviewTests(unittest.TestCase):
    def test_invalid_denominators_and_counts_are_unknown(self):
        for ac, an in [(0, 0), (1, None), (2, 1), (-1, 10), (True, 10), (1, 1.0)]:
            self.assertIsNone(count_af(ac, an))
        self.assertEqual(count_af(0, 100), 0)

    def test_joint_an_is_used_without_manually_combining_callsets(self):
        r = summarize_record(record())
        self.assertEqual(r['af'], 1 / 20000)
        self.assertEqual(r['frequency_source'], 'joint')
        self.assertNotIn('genome', r['datasets'])
        self.assertIsNone(r['datasets']['joint']['fafmax'])

    def test_filtered_joint_falls_back_to_largest_valid_individual_af(self):
        r = summarize_record(record(joint=callset(filters=['AC0']),
                                    exome=callset(ac=2, an=100), genome=callset(ac=3, an=100)))
        self.assertEqual((r['frequency_source'], r['af']), ('genome', .03))

    def test_all_filtered_or_zero_an_remain_unknown(self):
        r = summarize_record(record(joint=callset(an=0), exome=callset(filters=['RF'])))
        self.assertEqual(r['status'], 'matched_frequency_unusable')
        self.assertIsNone(r['af'])

    def test_groupmax_ignores_aggregate_remaining_and_bad_denominators(self):
        groups = [{'id': '', 'ac': 1, 'an': 1}, {'id': 'remaining', 'ac': 1, 'an': 1},
                  {'id': 'ami', 'ac': 1, 'an': 0}, {'id': 'nfe', 'ac': 1, 'an': 100},
                  {'id': 'afr', 'ac': 1, 'an': 1000}]
        r = summarize_record(record(joint=callset(populations=groups)))
        self.assertEqual((r['groupmax_af'], r['groupmax_population'], r['groupmax_an']),
                         (.01, 'nfe', 100))

    def test_unknown_does_not_fall_back_to_legacy_or_receive_rarity_bonus(self):
        c = candidate(); e = {'gene': 'TEST', 'variant_id': c['variant_id'], **summarize_record(None)}
        r = apply_frequencies([c], [e])[0]
        self.assertIsNone(r['population_af'])
        self.assertEqual(r['triage_score'], 5)
        self.assertEqual(c['population_af'], .2)

    def test_classification_ablation_preserves_independent_rarity(self):
        c = candidate(); e = {'gene': 'TEST', 'variant_id': c['variant_id'], **summarize_record(record())}
        r = apply_frequencies([c], [e], clinical=False)[0]
        self.assertEqual(r['population_af'], 1 / 20000)
        self.assertEqual(r['triage_score'], 6)
        self.assertEqual(r['clinical_score'], 0)
        self.assertEqual(r['clinvar'], [])
        self.assertEqual(c['clinical_score'], 3)

    def test_overlapping_gene_annotations_do_not_share_missing_match(self):
        a = candidate(); b = candidate(gene='OTHER')
        e = [{'gene': 'TEST', 'variant_id': a['variant_id'], **summarize_record(record())},
             {'gene': 'OTHER', 'variant_id': a['variant_id'], **summarize_record(None)}]
        rows = apply_frequencies([a, b], e)
        self.assertIsNotNone(rows[0]['population_af'])
        self.assertIsNone(rows[1]['population_af'])

    def test_exact_and_unique_minimal_representation_matching(self):
        gene = {'variants': [record()], 'coverage': {'exome': [{'pos': 99, 'mean': 30}]}}
        exact = match_candidates([candidate()], gene)[0]
        self.assertEqual(exact['match_type'], 'exact')
        self.assertEqual(exact['nearby_binned_coverage']['exome']['distance_bp'], 2)
        alt = candidate(variant_id='1:100:TA:TC', pos=100, ref='TA', alt='TC')
        self.assertEqual(match_candidates([alt], gene)[0]['match_type'], 'unique_minimal_representation')

    def test_ambiguous_minimal_representation_is_not_selected(self):
        other = record(variant_id='1-100-TA-TC', pos=100, ref='TA', alt='TC')
        query = candidate(variant_id='1:99:GTA:GTC', pos=99, ref='GTA', alt='GTC')
        r = match_candidates([query], {'variants': [record(), other]})[0]
        self.assertEqual(r['match_type'], 'ambiguous_minimal_representation')
        self.assertIsNone(r['af'])

    def test_symbolic_unrelated_alleles_do_not_match_through_none(self):
        query = candidate(variant_id='1:201:A:<DUP>', pos=201, alt='<DUP>')
        r = match_candidates([query], {'variants': [record(alt='<DEL>')]})[0]
        self.assertEqual(r['match_type'], 'unmatched')
        self.assertIsNone(r['af'])
        self.assertIsNone(normalize_variant('1', 0, 'A', 'G'))
        self.assertIsNone(normalize_variant('1', 1, 'A', 'A'))

    def test_wrong_assembly_and_duplicate_allele_fail_closed(self):
        for variants in [[record(reference_genome='GRCh37')], [record(), copy.deepcopy(record())]]:
            with self.assertRaises(ValueError):
                match_candidates([candidate()], {'variants': variants})


if __name__ == '__main__':
    unittest.main()
