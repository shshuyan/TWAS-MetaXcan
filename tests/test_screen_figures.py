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
            consistency = pd.read_csv(batch / "figures" / "figure_06_condition_consistency.data.tsv", sep="\t")
            sensitivity = pd.read_csv(batch / "figures" / "figure_08_context_sensitivity.data.tsv", sep="\t")
            self.assertLessEqual(consistency.groupby("disease").pert_id.nunique().max(), 10)
            self.assertLessEqual(sensitivity.groupby("disease").compound_id.nunique().max(), 10)

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

    def test_dense_panels_limit_displayed_compounds_but_keep_tidy_sidecars(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            batch = run_four_disease_screen(batch_fixture(root))
            for disease in ["AD", "LSC", "atopic_dermatitis", "IPF"]:
                disease_root = batch / disease
                for filename, id_column, name_column in [
                    ("signature_evidence.tsv.gz", "pert_id", "pert_iname"),
                    ("all_primary_tests.tsv", "compound_id", "drug"),
                    ("secondary_contexts.tsv", "compound_id", "drug"),
                ]:
                    path = disease_root / filename
                    base = pd.read_csv(path, sep="\t")
                    copies = []
                    for index in range(12):
                        copy = base.copy()
                        copy[id_column] = copy[id_column].astype(str) + f"_{index:02}"
                        copy[name_column] = copy[name_column].astype(str) + f"_{index:02}"
                        copies.append(copy)
                    pd.concat(copies, ignore_index=True).to_csv(path, sep="\t", index=False)

            render_figure_suite(batch, formats=("svg",))

            consistency = pd.read_csv(batch / "figures" / "figure_06_condition_consistency.data.tsv", sep="\t")
            sensitivity = pd.read_csv(batch / "figures" / "figure_08_context_sensitivity.data.tsv", sep="\t")
            self.assertLessEqual(consistency.groupby("disease").pert_id.nunique().max(), 10)
            self.assertLessEqual(sensitivity.groupby("disease").compound_id.nunique().max(), 10)


if __name__ == "__main__":
    unittest.main()
