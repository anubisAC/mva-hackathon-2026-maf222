import gzip
import json
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from mva.intake import inspect_data


class IntakeTests(unittest.TestCase):
    def prepare(self, directory):
        data = directory / "data"
        data.mkdir()
        (data / "test.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        cells = ["Parental recurrent miscarriages", "Recurrent spontaneous abortion", "HP:0200067", "Family history"]
        document = '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:tbl><w:tr>' + ''.join(f'<w:tc><w:p><w:r><w:t>{c}</w:t></w:r></w:p></w:tc>' for c in cells) + '</w:tr></w:tbl></w:body></w:document>'
        with ZipFile(data / "clinical.docx", "w") as archive:
            archive.writestr("word/document.xml", document)
        return data

    def test_truncated_mate_recorded_without_crash(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            data = self.prepare(directory)
            for mate, text in ((1, "@read1\nAC\n+\nII\n@read2\nAC\n+\nII\n"), (2, "@read1\nAC\n+\nII\n")):
                with gzip.open(data / f"sample_R{mate}_001.fastq.gz", "wt") as handle:
                    handle.write(text)
            result = inspect_data(data, directory / "out")
            self.assertFalse(result["fastq_pairs"][0]["valid_structure"])
            self.assertFalse(result["fastq_pairs"][0]["paired_identifiers_match"])
            features = json.loads((directory / "out/phenotype.json").read_text())["features"]
            self.assertEqual(features[0]["subject"], "family_history")

    def test_missing_mate_is_not_validated(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            data = self.prepare(directory)
            with gzip.open(data / "sample_R1_001.fastq.gz", "wt") as handle:
                handle.write("@read1\nAC\n+\nII\n")
            pair = inspect_data(data, directory / "out")["fastq_pairs"][0]
            self.assertFalse(pair["mate_exists"])
            self.assertFalse(pair["valid_structure"])


if __name__ == "__main__":
    unittest.main()
