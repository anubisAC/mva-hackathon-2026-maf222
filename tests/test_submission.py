import csv
import tempfile
import unittest
from pathlib import Path

from mva.submission import FIELDS, validate


class SubmissionTests(unittest.TestCase):
    def check_rows(self, rows):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "submission.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            return validate(path)

    def row(self, **updates):
        return {**dict.fromkeys(FIELDS, ""), "proband_id": "PROBAND01", "chrom_1": "chr1", "pos_1": "10", "ref_1": "A", "alt_1": "C", "epcr": "0.5", "finding_type": "primary", **updates}

    def test_valid_single(self):
        self.assertTrue(self.check_rows([self.row()])["valid"])

    def test_empty_epcr_and_nan_rejected(self):
        for value in ("", "NaN", "inf", "0", "2"):
            self.assertFalse(self.check_rows([self.row(epcr=value)])["valid"])

    def test_partial_pair(self):
        self.assertFalse(self.check_rows([self.row(chrom_2="chr1")])["valid"])

    def test_duplicate_and_reverse_pair(self):
        row = self.row(chrom_2="chr1", pos_2="20", ref_2="G", alt_2="T")
        reverse = self.row(pos_1="20", ref_1="G", alt_1="T", chrom_2="chr1", pos_2="10", ref_2="A", alt_2="C")
        self.assertFalse(self.check_rows([row, reverse])["valid"])

    def test_unsorted_epcr(self):
        self.assertFalse(self.check_rows([self.row(epcr="0.1"), self.row(pos_1="11", epcr="0.9")])["valid"])

    def test_wrong_identifier_and_coordinates(self):
        self.assertFalse(self.check_rows([self.row(proband_id="sample", chrom_1="1", pos_1="0")])["valid"])


if __name__ == "__main__":
    unittest.main()
