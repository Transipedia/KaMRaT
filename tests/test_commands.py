import unittest
from os import path, mkdir, listdir
import os
from shutil import rmtree
import subprocess


kamrat = path.join(".", "bin", "kamrat")


class TestCommands(unittest.TestCase):

    def test_basic_command(self):
        output = subprocess.check_output(f"./bin/kamrat", stderr=subprocess.STDOUT, shell=True, text=True)
        self.assertTrue("================== k-mer Matrix, Really Tremendous! ==================" in output)

    def test_index(self):
        test_dir = "tmp_test"
        data = path.join("toyroom", "data")

        # Remove previous test remainings
        if path.exists(test_dir):
            rmtree(test_dir)
        mkdir(test_dir)

        # Define inputs and outputs
        intab = path.join(data, "kmer-counts.subset4toy.tsv.gz")
        outdir = path.join(test_dir, "kamrat.idx")
        mkdir(outdir)
        index_stdout = path.join(test_dir, "index.stdout")

        cmd = f"{kamrat} index -intab {intab} -outdir {outdir} -klen 31 -unstrand -nfbase 1000000000"
        # Command ok ?
        process = None
        with open(index_stdout, "w") as idx_out:
            process = subprocess.run(cmd.split(" "), stdout=idx_out, stderr=idx_out)
        self.assertEqual(0, process.returncode)
        # Number of outputs ok ?
        self.assertEqual(3, len(listdir(outdir)))

        # correct outputs ?
        idx_mat = path.join(outdir, "idx-mat.bin")
        self.assertTrue(path.exists(idx_mat))
        stream = os.popen(f"md5sum {idx_mat}")
        self.assertTrue(stream.read().startswith("e5904e8c5e9cce3be07e3d96a9a91bc0"))
        stream.close()

        idx_mat = path.join(outdir, "idx-meta.bin")
        self.assertTrue(path.exists(idx_mat))
        stream = os.popen(f"md5sum {idx_mat}")
        self.assertTrue(stream.read().startswith("c6ac53928b1cd9ce946de29496684cd0"))
        stream.close()

        idx_mat = path.join(outdir, "idx-pos.bin")
        self.assertTrue(path.exists(idx_mat))
        stream = os.popen(f"md5sum {idx_mat}")
        self.assertTrue(stream.read().startswith("86eca32eb46f0f701a8ba2c56e6393aa"))
        stream.close()

        # Cleaning
        rmtree(test_dir)



if __name__ == '__main__':
  unittest.main()