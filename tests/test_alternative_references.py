"""Small public-format fixtures exercise reference conversion without patient data."""
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from mva.alternative_references import attributes, build_gene_models, build_phenotype_reference
from mva.track1 import GeneModel


class AlternativeReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def gtf(self, transcript, symbol, intervals, *, strand="+", chrom="1", canonical=True,
            gene_biotype="protein_coding", transcript_biotype="protein_coding"):
        tags = 'tag "gencode_basic"; '
        if canonical:
            tags += 'tag "Ensembl_canonical"; '
        attrs = (f'gene_id "GENE_{symbol}"; gene_name "{symbol}"; '
                 f'gene_biotype "{gene_biotype}"; transcript_id "{transcript}"; '
                 f'transcript_version "2"; transcript_biotype "{transcript_biotype}"; ' + tags)
        return [f"{chrom}\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"
                for kind, start, end in intervals]

    def models(self, rows, sequences):
        with gzip.open(self.root / "annotation.gtf.gz", "wt", encoding="utf-8") as stream:
            stream.write("#!genome-build GRCh38\n" + "\n".join(rows) + "\n")
        with gzip.open(self.root / "cds.fa.gz", "wt", encoding="utf-8") as stream:
            for transcript, sequence in sequences.items():
                # Versioned, wrapped FASTA records exercise both common source conventions.
                stream.write(f">{transcript}.2 cds chromosome:GRCh38\n{sequence[:7]}\n{sequence[7:]}\n")
        summary = build_gene_models(self.root)
        with gzip.open(self.root / "genes_grch38.jsonl.gz", "rt", encoding="utf-8") as stream:
            metadata = json.loads(next(stream))
            models = [json.loads(line) for line in stream if line.strip()]
        self.assertEqual(metadata["usable_models"], len(models))
        self.assertEqual(summary["usable_models"], len(models))
        for name in ("annotation.gtf.gz", "cds.fa.gz"):
            self.assertEqual(metadata["source_hashes"][name],
                             hashlib.sha256((self.root / name).read_bytes()).hexdigest())
        return summary, {model["symbol"]: model for model in models}

    def test_repeated_gtf_tags_preserve_canonical_selection(self):
        parsed = attributes('gene_name "TEST"; tag "gencode_basic"; tag "MANE_Select"; '
                            'tag "Ensembl_canonical"; transcript_id "TX";')
        self.assertEqual(parsed["tag"], ["gencode_basic", "MANE_Select", "Ensembl_canonical"])
        self.assertEqual(parsed["transcript_id"], "TX")

    def test_cds_includes_stop_and_excludes_utr_on_both_strands(self):
        plus = self.gtf("TX_PLUS", "PLUS", [
            ("transcript", 90, 210), ("exon", 90, 105), ("CDS", 100, 105),
            ("exon", 200, 210), ("CDS", 200, 205), ("stop_codon", 206, 208)])
        minus = self.gtf("TX_MINUS", "MINUS", [
            ("transcript", 298, 410), ("exon", 298, 308), ("stop_codon", 300, 302),
            ("CDS", 303, 308), ("exon", 400, 410), ("CDS", 400, 405)], strand="-")
        summary, models = self.models(plus + minus,
                                     {"TX_PLUS": "ATGGAATTTCCCTAA", "TX_MINUS": "ATGGAATTTCCCTAA"})
        self.assertEqual(summary["usable_models"], 2)
        self.assertEqual(models["PLUS"]["translation_start"], 100)
        self.assertEqual(models["PLUS"]["translation_end"], 208)
        self.assertEqual(models["MINUS"]["translation_start"], 300)
        self.assertEqual(models["MINUS"]["translation_end"], 405)
        plus_gene, minus_gene = GeneModel(models["PLUS"]), GeneModel(models["MINUS"])
        self.assertTrue(plus_gene.sequence_valid)
        self.assertTrue(minus_gene.sequence_valid)
        self.assertEqual([plus_gene.cds_offset(p) for p in (100, 105, 200, 208)], [0, 5, 6, 14])
        self.assertEqual([minus_gene.cds_offset(p) for p in (405, 400, 308, 300)], [0, 5, 6, 14])
        self.assertIsNone(plus_gene.cds_offset(99))
        self.assertIsNone(minus_gene.cds_offset(406))

    def test_terminal_stop_split_across_exons_is_preserved(self):
        rows = self.gtf("TX_SPLIT", "SPLIT", [
            ("transcript", 500, 601), ("exon", 500, 509), ("CDS", 500, 508),
            ("stop_codon", 509, 509), ("exon", 600, 601), ("stop_codon", 600, 601)])
        _, models = self.models(rows, {"TX_SPLIT": "ATGAAACCCTAA"})
        gene = GeneModel(models["SPLIT"])
        self.assertTrue(gene.sequence_valid)
        self.assertEqual([gene.cds_offset(p) for p in (509, 600, 601)], [9, 10, 11])
        self.assertEqual(gene.cds[-3:], "TAA")

    def test_inconsistent_sequence_length_is_excluded_with_reason(self):
        intervals = [("transcript", 100, 111), ("exon", 100, 111),
                     ("CDS", 100, 108), ("stop_codon", 109, 111)]
        summary, models = self.models(self.gtf("TX_BAD", "BAD", intervals), {"TX_BAD": "ATGAAATAA"})
        self.assertEqual(models, {})
        self.assertEqual(summary["skipped_counts"].get("cds_length_mismatch"), 1)
        detail = json.loads((self.root / "gene_models_summary.json").read_text())
        self.assertEqual(detail["skipped_models"]["cds_length_mismatch"][0]["transcript"], "TX_BAD")

    def test_only_selected_nuclear_protein_coding_transcripts_are_retained(self):
        intervals = [("transcript", 100, 111), ("exon", 100, 111),
                     ("CDS", 100, 108), ("stop_codon", 109, 111)]
        rows = []
        cases = [("GOOD", {}), ("NONCANONICAL", {"canonical": False}),
                 ("MITO", {"chrom": "MT"}), ("NONCODING", {"transcript_biotype": "retained_intron"}),
                 ("PSEUDOGENE", {"gene_biotype": "processed_pseudogene"})]
        for symbol, options in cases:
            rows.extend(self.gtf("TX_" + symbol, symbol, intervals, **options))
        _, models = self.models(rows, {"TX_" + symbol: "ATGGAATTTTAA" for symbol, _ in cases})
        self.assertEqual(set(models), {"GOOD"})

    def test_internal_stop_excluded_but_out_of_frame_stop_sequence_allowed(self):
        intervals = [("transcript", 100, 111), ("exon", 100, 111),
                     ("CDS", 100, 108), ("stop_codon", 109, 111)]
        rows = self.gtf("TX_RECODED", "RECODED", intervals)
        rows += self.gtf("TX_STANDARD", "STANDARD", intervals)
        summary, models = self.models(rows, {"TX_RECODED": "ATGTGAAAATAA", "TX_STANDARD": "ATGATAACCTAA"})
        self.assertEqual(set(models), {"STANDARD"})
        self.assertEqual(summary["skipped_counts"].get("internal_stop_or_recoded_transcript"), 1)

    def test_duplicate_symbols_are_all_excluded_without_selecting_an_arbitrary_locus(self):
        intervals = [("transcript", 100, 111), ("exon", 100, 111),
                     ("CDS", 100, 108), ("stop_codon", 109, 111)]
        rows = self.gtf("TX_COPY_X", "DUPLICATE", intervals, chrom="X")
        rows += [row.replace('gene_id "GENE_DUPLICATE"', 'gene_id "OTHER_LOCUS"')
                 for row in self.gtf("TX_COPY_Y", "DUPLICATE", intervals, chrom="Y")]
        rows += self.gtf("TX_UNIQUE", "UNIQUE", intervals)
        summary, models = self.models(rows, {transcript: "ATGGAATTTTAA"
                                             for transcript in ("TX_COPY_X", "TX_COPY_Y", "TX_UNIQUE")})
        self.assertEqual(set(models), {"UNIQUE"})
        self.assertEqual(summary["skipped_counts"].get("ambiguous_duplicate_symbol"), 2)
        detail = json.loads((self.root / "gene_models_summary.json").read_text())
        self.assertEqual({row["transcript"] for row in detail["skipped_models"]["ambiguous_duplicate_symbol"]},
                         {"TX_COPY_X", "TX_COPY_Y"})

    def test_phenotype_deduplication_and_family_provenance(self):
        hpo = self.root / "phenotype_to_genes.txt"
        hpo.write_text("hpo_id\thpo_name\tncbi_gene_id\tgene_symbol\tdisease_id\n"
                       "HP:0000001\tFeature one\t1\tGENE_A\tOMIM:1\n"
                       "HP:0000001\tFeature one\t1\tGENE_A\tORPHA:2\n"
                       "HP:0000002\tFeature two\t1\tGENE_A\tOMIM:1\n"
                       "HP:0200067\tFamily feature\t2\tGENE_B\tOMIM:3\n", encoding="utf-8")
        phenotype = self.root / "phenotype.json"
        phenotype.write_text(json.dumps({"sha256": "synthetic-clinical-hash", "features": [
            {"hpo_id": "HP:0000001", "label": "Feature one", "subject": "proband"},
            {"hpo_id": "HP:0000999", "label": "Unmatched feature", "subject": "proband"},
            {"hpo_id": "HP:0200067", "label": "Family feature", "subject": "family_history"}]}))
        output = self.root / "nested" / "phenotype_reference.json"
        summary = build_phenotype_reference(hpo, phenotype, output)
        result = json.loads(output.read_text())
        self.assertEqual(result["genes"]["GENE_A"], ["HP:0000001", "HP:0000002"])
        self.assertEqual(result["genes"]["GENE_B"], ["HP:0200067"])
        self.assertNotIn("HP:0000999", result["genes"]["GENE_A"])
        self.assertEqual([row["hpo_id"] for row in result["patient_features"]],
                         ["HP:0000001", "HP:0000999"])
        self.assertEqual(result["family_features"], [
            {"hpo_id": "HP:0200067", "label": "Family feature", "source": "family_history"}])
        self.assertEqual(summary["association_rows"], 4)
        self.assertEqual(summary["gene_count"], 2)
        self.assertEqual(result["source"]["clinical_sha256"], "synthetic-clinical-hash")
        self.assertEqual(result["source"]["hpo_sha256"], hashlib.sha256(hpo.read_bytes()).hexdigest())
        self.assertEqual(result["gene_inheritance"], {})
        self.assertIsNone(result["source"]["inheritance_sha256"])

    def test_inheritance_unions_disease_modes_but_preserves_sources_and_excludes_negation(self):
        hpo = self.root / "phenotype_to_genes.txt"
        hpo.write_text("hpo_id\thpo_name\tncbi_gene_id\tgene_symbol\tdisease_id\n"
                       "HP:0001111\tFeature\t1\tDUAL\tOMIM:1\n"
                       "HP:0001111\tFeature\t1\tDUAL\tORPHA:2\n"
                       "HP:0002222\tOther feature\t1\tDUAL\tOMIM:1\n"
                       "HP:0001111\tFeature\t2\tRECESSIVE\tOMIM:1\n"
                       "HP:0001111\tFeature\t3\tNEGATED\tOMIM:3\n"
                       "HP:0001111\tFeature\t4\tUNKNOWN\tOMIM:4\n", encoding="utf-8")
        hpoa = self.root / "phenotype.hpoa"
        # The pinned public release uses comment metadata followed by an uncommented header.
        hpoa.write_text('#description: "Synthetic disease annotations"\n#version: 2026-09-02\n'
                        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\taspect\n"
                        "OMIM:1\tDisease one\t\tHP:0000007\tSYNTHETIC\tTAS\tI\n"
                        "OMIM:1\tDisease one\t\tHP:0000007\tSYNTHETIC\tTAS\tI\n"
                        "ORPHA:2\tDisease two\t\tHP:0000006\tSYNTHETIC\tTAS\tI\n"
                        "OMIM:3\tDisease three\tNOT\tHP:0000007\tSYNTHETIC\tTAS\tI\n"
                        "OMIM:4\tDisease four\t\tHP:0001111\tSYNTHETIC\tTAS\tP\n", encoding="utf-8")
        phenotype = self.root / "phenotype.json"
        phenotype.write_text('{"features": []}')
        output = self.root / "out.json"
        summary = build_phenotype_reference(hpo, phenotype, output)
        result = json.loads(output.read_text())
        self.assertEqual(result["gene_inheritance"], {
            "DUAL": ["autosomal_dominant", "autosomal_recessive"],
            "RECESSIVE": ["autosomal_recessive"]})
        self.assertEqual(result["gene_inheritance_sources"]["DUAL"], {
            "OMIM:1": ["autosomal_recessive"], "ORPHA:2": ["autosomal_dominant"]})
        self.assertEqual(summary["inheritance_genes"], 2)
        self.assertEqual(result["source"]["inheritance_sha256"], hashlib.sha256(hpoa.read_bytes()).hexdigest())

    def test_unrecognized_hpo_data_fails_instead_of_creating_empty_reference(self):
        hpo = self.root / "bad_hpo.tsv"
        hpo.write_text("wrong\theader\nnot\tan-hpo-row\n")
        phenotype = self.root / "phenotype.json"
        phenotype.write_text('{"features": []}')
        with self.assertRaisesRegex(ValueError, "No gene/HPO associations"):
            build_phenotype_reference(hpo, phenotype, self.root / "out.json")


if __name__ == "__main__":
    unittest.main()
