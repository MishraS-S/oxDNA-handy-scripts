import unittest
from pathlib import Path
import argparse
import json
import json_to_sqs_tacoxdna


class TestJsonToSqs(unittest.TestCase):

    sample_designs = Path("tests/sample_designs").resolve()

    def test_compl(self):
        result = json_to_sqs_tacoxdna.compl("A")
        expected = "T"
        self.assertEqual(result, expected)
        result = json_to_sqs_tacoxdna.compl("T")
        expected = "A"
        self.assertEqual(result, expected)
        result = json_to_sqs_tacoxdna.compl("C")
        expected = "G"
        self.assertEqual(result, expected)
        result = json_to_sqs_tacoxdna.compl("G")
        expected = "C"
        self.assertEqual(result, expected)

    def test_compute_start(self):
        result = json_to_sqs_tacoxdna.compute_start(self.sample_designs / "48bp-lin.json")
        expected = [0, 8]
        self.assertEqual(result, expected)
        result = json_to_sqs_tacoxdna.compute_start(self.sample_designs / "32bp-lin.json")
        expected = [0, 16]
        self.assertEqual(result, expected)
        result = json_to_sqs_tacoxdna.compute_start(self.sample_designs / "32bp-lin_vhs1.json")
        expected = [1, 31]
        self.assertEqual(result, expected)
        result = json_to_sqs_tacoxdna.compute_start(self.sample_designs / "32bp-circ.json")
        expected = [-1, -1]
        self.assertEqual(result, expected)

    def test_load_data(self):
        json_path = self.sample_designs / "48bp-lin.json"
        seq_file = self.sample_designs / "48bp-seq.txt"
        result = json_to_sqs_tacoxdna.load_data(
            args=argparse.Namespace(json=json_path, sequence=seq_file, vHelix=-1, column=-1, out_name="caca.sqs"))
        with open(json_path, "r") as j:
            json_input = json.load(j)
        with open(seq_file, "r") as seq_file:
            seq = seq_file.read()
            seq = seq.upper()
            seq = list(seq)
        expected = [json_input, list("AGACTTCCGGCTTAAGCTCTGAAAGGGTTCTATATCTCCAGGTAGATC"), 0, 8]
        self.assertEqual(result, expected)

    def test_compute_output(self):
        with open(self.sample_designs / "48bp-lin.json", "r") as j:
            json_input = json.load(j)
        data = [json_input, list("AGACTTCCGGCTTAAGCTCTGAAAGGGTTCTATATCTCCAGGTAGATC"), 0, 8]
        result = json_to_sqs_tacoxdna.compute_output(data=data)
        expected = [list("RRRRRRRR" + "AGACTTCCGGCTTAAGCTCTGAAA" + "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"), list("RRRRRRRR" + "GGGTTCTATATCTCCAGGTAGATC"[::-1] + "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR")]
        self.assertEqual(result, expected)

        with open(self.sample_designs / "32bp-circ.json", "r") as j:
            json_input = json.load(j)
        data = [json_input, list("GGCTTAAGCTCTGAAAGGGTTCTATATCTCCA"), 0, 16]
        result = json_to_sqs_tacoxdna.compute_output(data=data)
        expected = [list("RRRRRRRRRRRRRRRR" + "GGCTTAAGCTCTGAAA" + "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"), list("RRRRRRRRRRRRRRRR" + "GGGTTCTATATCTCCA"[::-1] + "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR")]
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
