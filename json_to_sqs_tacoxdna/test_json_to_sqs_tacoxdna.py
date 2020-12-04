import unittest
import json_to_sqs_tacoxdna
from pathlib import Path
import argparse
import json


class TestJsonToSqs(unittest.TestCase):
    # TODO: read documentation about all the different testing possibilities
    folder = Path(
        "/Users/peiskert/Documents/TUM/Master-thesis/Simulations/free-form-DNA-origami/json_to_sqs_tacoxdna/")
    seq_file = Path("1033bp-seq.txt")
    json_file = Path("sample_designs/48bp-lin.json")
    sample_designs = Path("sample_designs").resolve()

    def test_proc_input(self):
        vHelix = 0
        column = 5
        args = argparse.Namespace(json=self.json_file, sequence = self.seq_file, vHelix = vHelix, column = column)
        result = json_to_sqs_tacoxdna.proc_input()  # TODO figure out how to do unittesting on input parameters
        expected = ""
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


if __name__ == "__main__":
    unittest.main()