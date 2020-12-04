import unittest
import json_to_sqs_tacoxdna
from pathlib import Path
import argparse


class TestJsonToSqs(unittest.TestCase):
    # TODO: read documentation about all the different testing possibilities
    folder = Path(
        "/Users/peiskert/Documents/TUM/Master-thesis/Simulations/free-form-DNA-origami/json_to_sqs_tacoxdna/")
    seq_file = Path("1033bp-seq.txt")
    json_file = Path("p-ff2.json")

    def test_proc_input(self):
        vHelix = 0
        column = 5
        args = argparse.Namespace(json=self.json_file, sequence = self.seq_file, vHelix = vHelix, column = column)
        result = json_to_sqs_tacoxdna.proc_input()  # TODO figure out how to do unittesting on input parameters
        expected = ""
        self.assertEqual(result, expected)
    def test_load_data(self):
        with open(self.seq_file) as seq:
            result = ""
            expected = ""
            self.assertEqual(result, expected)



if __name__ == "__main__":
    unittest.main()