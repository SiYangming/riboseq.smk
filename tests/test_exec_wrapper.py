import unittest

from workflow.scripts import exec_wrapper as ew


class TestExecWrapper(unittest.TestCase):
    def test_container_command_building(self):
        image = "docker://quay.io/biocontainers/samtools:1.23--h96c455f_0"
        built = ew.build_container_command(
            image, "/host:/host", "/host", ["samtools", "--version"]
        )
        self.assertEqual(built[0], "apptainer")
        self.assertEqual(built[1], "exec")
        self.assertIn(image, built)
        self.assertIn("--bind", built)

    def test_native_and_conda(self):
        cfg = {
            "exec_mode": "native",
            "fastqc": {"fastqc_bin": "/opt/fastqc", "container_image": ""},
        }
        prefix, binary = ew.exec_wrapper_binary(cfg, "fastqc", "fastqc_bin", "fastqc")
        self.assertEqual(prefix, "")
        self.assertEqual(binary, "/opt/fastqc")
        cfg["exec_mode"] = "conda"
        prefix, binary = ew.exec_wrapper_binary(cfg, "fastqc", "fastqc_bin", "fastqc")
        self.assertEqual(prefix, "")
        self.assertEqual(binary, "fastqc")

    def test_container_mode(self):
        cfg = {
            "exec_mode": "container",
            "fastqc": {"container_image": "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"},
        }
        prefix, binary = ew.exec_wrapper_binary(cfg, "fastqc", "fastqc_bin", "fastqc")
        self.assertIn("apptainer exec", prefix)
        self.assertIn("docker://quay.io/biocontainers/fastqc", prefix)
        self.assertEqual(binary, "fastqc")

    def test_invalid_mode(self):
        with self.assertRaises(ValueError):
            ew.exec_wrapper_binary(
                {"exec_mode": "docker", "fastqc": {}}, "fastqc", "fastqc_bin", "fastqc"
            )


if __name__ == "__main__":
    unittest.main()
