"""Figure-suite outputs are stable artifacts, not alternate analyses."""
from pathlib import Path
import tempfile
import unittest
import warnings

import pandas as pd

from tests.test_four_disease_screen import batch_fixture
from twaslib.four_disease_screen import run_four_disease_screen
from twaslib.screen_figures import FIGURES, render_figure_suite


class ScreenFigureTests(unittest.TestCase):
    def test_all_figures_render_in_three_formats_with_demo_disclosure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            batch = run_four_disease_screen(batch_fixture(root))
            paths = render_figure_suite(batch)

            expected = {f"{stem}.{extension}" for stem in FIGURES for extension in ["png", "svg", "pdf"]}
            self.assertEqual({path.name for path in paths}, expected)
            self.assertTrue(all(path.stat().st_size > 1000 for path in paths))
            for stem in FIGURES:
                self.assertTrue((batch / "figures" / f"{stem}.data.tsv").is_file())
                svg = (batch / "figures" / f"{stem}.svg").read_text()
                self.assertIn("DEMO / synthetic", svg)
            captions = (batch / "figures" / "figure_captions.md").read_text()
            self.assertEqual(captions.count("DEMO / synthetic"), len(FIGURES))
            self.assertIn("BH-FDR", captions)
            candidate_svg = (batch / "figures" / "figure_07_cross_disease_candidates.svg").read_text()
            self.assertEqual(candidate_svg.count("Median primary-cell score"), 1)

    def test_empty_candidate_table_renders_an_explicit_no_data_panel(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            batch = run_four_disease_screen(batch_fixture(root))
            candidate_path = batch / "LSC" / "candidate_compounds.tsv"
            pd.read_csv(candidate_path, sep="\t").iloc[:0].to_csv(candidate_path, sep="\t", index=False)

            render_figure_suite(batch, formats=("svg",))

            svg = (batch / "figures" / "figure_07_cross_disease_candidates.svg").read_text()
            self.assertIn("No candidates", svg)
            self.assertIn("Lichen simplex chronicus", svg)

    def test_pdf_render_has_no_third_party_deprecation_noise(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            batch = run_four_disease_screen(batch_fixture(root))
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                render_figure_suite(batch, formats=("pdf",))
            self.assertEqual([str(item.message) for item in caught if issubclass(item.category, DeprecationWarning)], [])


if __name__ == "__main__":
    unittest.main()
