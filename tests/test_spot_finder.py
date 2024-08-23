from .conftest import LOCAL_RESOURCES_DIR
import pytest

import numpy as np
from math import isclose

from cellspot.spot_finder import SpotFinder


def test_spot_finder():
    sf = SpotFinder(image_path= LOCAL_RESOURCES_DIR / "cell_image_001.tif")
    masks = sf.create_cell_mask()
    assert True


def test_ndarray_float_division():
    x = np.array([i for i in range(10, 0, -1)])
    y = np.array([i+1 for i in range(10)])
    div = np.divide(x, y, dtype=np.float64)
    assert isclose(div[0], 10.0, abs_tol=1e-3)
    assert isclose(div[1], 4.5, abs_tol=1e-3)
    assert isclose(div[2], 8.0/3.0, abs_tol=1e-4)


@pytest.mark.parametrize("img, msk, ans",
                         [  (np.array([10, 9, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 6,
                                      5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4,
                                      3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], dtype=np.int32),
                             np.array([10, 9, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 6,
                                      5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4,
                                      3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], dtype=np.int32),
                             np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], dtype=np.float32)
                             )
                            ]
                         )
def test_create_average_intensities(img, msk, ans):
    calc = SpotFinder._create_average_intensities(nuclei_image=img, cell_mask=msk)
    np.testing.assert_array_equal(calc, ans)

