"""Direct statistical regression checks retained from the superseded extension."""
import unittest

import numpy as np
from scipy.stats import false_discovery_control

from twaslib.screen_stats import bh, cmap_calibration, cmap_connectivity, lincs_calibration


class ScreenStatisticsTests(unittest.TestCase):
    def test_bh_matches_known_values_and_scipy(self):
        p = [.01, .04, .03, 1]
        self.assertTrue(np.allclose(bh(p), [.04, .053333333, .053333333, 1]))
        self.assertTrue(np.allclose(bh(p), false_discovery_control(p)))
        self.assertEqual(len(bh([])), 0)

    def test_ks_rank_direction_and_same_sign_zero(self):
        ranks = np.arange(1, 101, dtype=float)[None, :]
        self.assertGreater(cmap_connectivity(ranks, list(range(10)), list(range(90, 100)))[0], .8)
        self.assertLess(cmap_connectivity(ranks, list(range(90, 100)), list(range(10)))[0], -.8)
        self.assertEqual(cmap_connectivity(ranks, list(range(10)), list(range(10, 20)))[0], 0)

    def test_calibration_preserves_three_or_more_group_dimensions(self):
        signs = np.array([[1, -1, 0, 0], [-1, 1, 0, 0], [0, 0, 1, -1]])
        groups = [[0], [1], [2]]
        _, observed, p, q = lincs_calibration(signs, [0], [1], groups, 100, 1)
        self.assertEqual((len(observed), len(p), len(q)), (3, 3, 3))
        ranks = np.array([[1, 4, 2, 3], [4, 1, 2, 3], [2, 3, 1, 4]])
        observed, p, q = cmap_calibration(ranks, [0], [1], np.arange(4), groups, 100, 1)
        self.assertEqual((len(observed), len(p), len(q)), (3, 3, 3))


if __name__ == "__main__":
    unittest.main()
