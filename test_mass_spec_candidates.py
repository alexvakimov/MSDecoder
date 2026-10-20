import unittest
import contextlib
import io
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from mass_spec_candidates import (
    ADDUCT_MASSES,
    DEMO_PEAKS,
    DEMO_PRECURSOR_MZ,
    Peak,
    _score_fragment_spectrum,
    analyze_text_file,
    enumerate_formulas,
    formula_mass,
    main,
    rank_candidates,
    read_text_spectrum,
)


ROOT = Path(__file__).parent


class MassSpecCandidateTests(unittest.TestCase):
    def test_wide_tolerance_enumerates_all_hydrogen_counts(self):
        matches = enumerate_formulas(
            15.0, 2.0, max_counts={"C": 1}, hydrogen_limit=4,
            allowed_elements=("C", "H"),
        )
        self.assertEqual({counts["H"] for counts, _ in matches}, {2, 4})

    def test_enumeration_rejects_invalid_mass_and_bounds(self):
        for mass, tolerance in ((float("nan"), 1), (12, float("inf"))):
            with self.subTest(mass=mass, tolerance=tolerance):
                with self.assertRaises(ValueError):
                    enumerate_formulas(mass, tolerance)
        with self.assertRaisesRegex(ValueError, "bounds"):
            enumerate_formulas(12, 1, hydrogen_limit=-1)

    def test_empty_required_element_iterator(self):
        candidates = rank_candidates(
            DEMO_PRECURSOR_MZ, (), tolerance_ppm=50,
            required_elements=iter(()),
        )
        self.assertTrue(candidates)

    def test_cli_tolerance_overrides_text_metadata(self):
        with TemporaryDirectory() as temp_dir:
            args = [
                "mass_spec_candidates.py", "--input",
                str(ROOT / "spectra" / "slide02.txt"),
                "--output-dir", temp_dir, "--tolerance-ppm", "10",
            ]
            with patch.object(sys, "argv", args), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(), 0)
            report = (Path(temp_dir) / "report.txt").read_text()
            self.assertIn("Mass tolerance: 10 ppm", report)

    def test_cli_demo_accepts_element_restrictions(self):
        args = [
            "mass_spec_candidates.py", "--demo", "--elements", "C,H,N,O",
            "--require-elements", "N",
        ]
        output = io.StringIO()
        with patch.object(sys, "argv", args), contextlib.redirect_stdout(output):
            self.assertEqual(main(), 0)
        self.assertIn("C16H20N2O", output.getvalue())
        self.assertIn("Allowed elements: C, H, N, O", output.getvalue())

    def test_demo_rejects_ignored_peak_inputs(self):
        with patch.object(sys, "argv", ["mass_spec_candidates.py", "--demo", "--peak", "100:1"]):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                main()
        self.assertEqual(error.exception.code, 2)

    def test_empty_required_elements_override_metadata(self):
        with TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "input.txt"
            source.write_text(
                (ROOT / "spectra" / "slide02.txt").read_text()
                + "\n# required_elements: Br\n"
            )
            with contextlib.redirect_stdout(io.StringIO()):
                candidates = analyze_text_file(
                    source, Path(temp_dir) / "out", required_elements=(),
                )
            self.assertEqual(candidates[0].formula, "C16H20N2O")

    def test_first_spectrum_recovers_structure_formula(self):
        candidates = rank_candidates(
            DEMO_PRECURSOR_MZ, DEMO_PEAKS, tolerance_ppm=50, limit=10
        )

        self.assertTrue(candidates)
        self.assertEqual(candidates[0].formula, "C16H20N2O")
        self.assertAlmostEqual(candidates[0].exact_mass, 256.1576, places=3)
        self.assertTrue(any(item.startswith("225.14:C15H16N2") for item in candidates[0].explained_peaks))

    def test_digitized_slide_files_recover_reference_formulas(self):
        expected = {
            "slide02.txt": "C16H20N2O",
            "slide03.txt": "C10H15N2OBr",
            "slide04.txt": "C10H11Br",
            "slide05.txt": "C10H11Br",
            "slide06.txt": "C11H22N2O2",
            "slide07.txt": "C11H22N2O2",
            "slide08.txt": "C9H12",
            "slide09.txt": "C8H8N2",
            "slide10.txt": "C9H13N2Br",
            "slide11.txt": "C12H14N2O",
        }
        for filename, expected_formula in expected.items():
            with self.subTest(slide=filename):
                metadata, peaks = read_text_spectrum(ROOT / "spectra" / filename)
                candidates = rank_candidates(
                    float(metadata["precursor_mz"]),
                    peaks,
                    metadata.get("adduct", "proton"),
                    float(metadata.get("tolerance_ppm", 25)),
                    limit=10,
                )
                self.assertIn(expected_formula, {candidate.formula for candidate in candidates})

    def test_element_whitelist_and_required_bromine(self):
        metadata, peaks = read_text_spectrum(ROOT / "spectra" / "slide03.txt")
        precursor_mz = float(metadata["precursor_mz"])
        tolerance_ppm = float(metadata["tolerance_ppm"])
        allowed_elements = ("C", "H", "N", "O", "Br")

        candidates = rank_candidates(
            precursor_mz,
            peaks,
            metadata["adduct"],
            tolerance_ppm,
            allowed_elements=allowed_elements,
            required_elements=("Br",),
        )

        self.assertTrue(candidates)
        self.assertEqual(candidates[0].formula, "C10H15N2OBr")
        self.assertTrue(all("Br" in candidate.formula for candidate in candidates))
        self.assertTrue(all(not any(element in candidate.formula for element in ("S", "P", "Cl")) for candidate in candidates))

    def test_required_element_must_be_allowed(self):
        with self.assertRaisesRegex(ValueError, "must also be included"):
            rank_candidates(
                259.04,
                (),
                allowed_elements=("C", "H", "N", "O"),
                required_elements=("Br",),
            )

    def test_fragment_mse_uses_peak_intensities_as_weights(self):
        parent = {"C": 16, "H": 20, "N": 2, "O": 1}
        exact_fragment_mz = formula_mass({"C": 15, "H": 16, "N": 2}) + ADDUCT_MASSES["proton"]
        precursor_mz = formula_mass(parent) + ADDUCT_MASSES["proton"]
        unmatched_peak_mz = 300.123456

        strong_match = _score_fragment_spectrum(
            parent,
            (Peak(exact_fragment_mz, 100), Peak(unmatched_peak_mz, 1)),
            precursor_mz,
            ADDUCT_MASSES["proton"],
            50,
            ("C", "H", "N", "O"),
        )
        strong_mismatch = _score_fragment_spectrum(
            parent,
            (Peak(exact_fragment_mz, 1), Peak(unmatched_peak_mz, 100)),
            precursor_mz,
            ADDUCT_MASSES["proton"],
            50,
            ("C", "H", "N", "O"),
        )

        self.assertLess(strong_match[0], strong_mismatch[0])
        self.assertGreater(strong_match[1], strong_mismatch[1])

    def test_text_input_writes_plot_report_and_candidates(self):
        input_path = ROOT / "spectra" / "slide02.txt"
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "output"
            analyze_text_file(input_path, output_dir)

            self.assertTrue((output_dir / "spectrum.png").is_file())
            self.assertTrue((output_dir / "candidates.csv").is_file())
            candidates_csv = (output_dir / "candidates.csv").read_text(encoding="utf-8")
            self.assertIn("fragment_mse,matched_intensity_fraction", candidates_csv.splitlines()[0])
            report = (output_dir / "report.txt").read_text(encoding="utf-8")
            self.assertIn("C16H20N2O", report)
            self.assertIn("Allowed elements: C, H, N, O, S, P, Cl, Br", report)

    def test_batch_generates_summary_and_structure_images(self):
        from mass_spec_candidates import main
        import sys
        from unittest.mock import patch

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "results"
            args = [
                "mass_spec_candidates.py",
                "--batch",
                str(ROOT / "spectra"),
                "--structures",
                str(ROOT / "references" / "slide_structures.csv"),
                "--output-dir",
                str(output_dir),
            ]
            with patch.object(sys, "argv", args):
                self.assertEqual(main(), 0)

            self.assertTrue((output_dir / "summary.csv").is_file())
            self.assertTrue((output_dir / "slide02" / "structure.png").is_file())
            summary_header = (output_dir / "summary.csv").read_text(encoding="utf-8").splitlines()[0]
            self.assertIn("top_fragment_mse", summary_header)
            self.assertIn("top_matched_intensity_fraction", summary_header)


if __name__ == "__main__":
    unittest.main()
