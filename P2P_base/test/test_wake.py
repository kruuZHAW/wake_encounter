import unittest
import pandas as pd
import wake.wake as wake


class test_wake(unittest.TestCase):
    # Test rotate_and_translate function
    def test_rotate_and_translate(self):
        # no rotation, no translation
        expected = (1, 2)
        returned = wake.rotate_and_translate(x=1, y=2, angle_deg=0, dx=0, dy=0)
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

        # no rotation, translation
        expected = (2, 4)
        returned = wake.rotate_and_translate(x=1, y=2, angle_deg=0, dx=1, dy=2)
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

        # rotation, no translation
        expected = (-2, 1)
        returned = wake.rotate_and_translate(
            x=1, y=2, angle_deg=90, dx=0, dy=0
        )
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

        # rotation, translation
        expected = (10.707, 5.707)
        returned = wake.rotate_and_translate(
            x=1, y=0, angle_deg=45, dx=10, dy=5
        )
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

       # negative rotation, no translation
        expected = (2, -1)
        returned = wake.rotate_and_translate(x=1, y=2, angle_deg=-90, dx=0, dy=0)
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

        # rotation with angle greater than 360 degrees, no translation
        expected = (-2, 1)
        returned = wake.rotate_and_translate(x=1, y=2, angle_deg=450, dx=0, dy=0)
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

        # no rotation, negative translation
        expected = (-1, 0)
        returned = wake.rotate_and_translate(x=1, y=2, angle_deg=0, dx=-2, dy=-2)
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

        # rotation, negative translation
        expected = (-8.707, 0.293)
        returned = wake.rotate_and_translate(x=1, y=0, angle_deg=45, dx=-10, dy=0)
        self.assertAlmostEqual(expected[0], returned[0], places=3)
        self.assertAlmostEqual(expected[1], returned[1], places=3)

    # Test change_ref_xy function
    def test_change_ref_xy(self):
        position_x = 100
        position_y = 50
        angle_deg = 90
        row_input = pd.Series({"x": 10, "yl": -10, "yr": 10})
        row_expected = pd.Series(
            {"xl_ref": 110, "yl_ref": 60, "xr_ref": 90, "yr_ref": 60}
        )
        returned = wake.change_ref_xy(
            row_input, position_x, position_y, angle_deg
        )
        self.assertAlmostEqual(row_expected["xl_ref"], returned["xl_ref"], places=3)
        self.assertAlmostEqual(row_expected["yl_ref"], returned["yl_ref"], places=3)
        self.assertAlmostEqual(row_expected["xr_ref"], returned["xr_ref"], places=3)
        self.assertAlmostEqual(row_expected["yl_ref"], returned["yl_ref"], places=3)
