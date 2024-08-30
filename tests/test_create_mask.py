import numpy as np
import pytest
from pathlib import Path

from attr import dataclass

from .conftest import LOCAL_RESOURCES_DIR
from cellspot.scripts.create_mask import main_with_args

RED_THRESH = 82  # 82
CLIP_LOW = 5
CLIP_HIGH = 90
CYCLES = 6

@dataclass
class TestArgs():
    filename: Path
    output_folder: Path
    clip_to: tuple[float, float]
    red_threshold: int
    pixel_cycles: int
    dead_cells: int

@dataclass
class ExpectedResult():
    cells: int
    cells_with_inclusions: int
    tolerance: float

    @property
    def ratio(self) -> float:
        return float(self.cells_with_inclusions) / float(self.cells)


@pytest.mark.parametrize( 'input,expected', [
    (TestArgs(filename=LOCAL_RESOURCES_DIR / "cell_image_001.tif", output_folder=LOCAL_RESOURCES_DIR,
              clip_to=(CLIP_LOW, CLIP_HIGH), red_threshold=RED_THRESH, pixel_cycles=CYCLES, dead_cells=70),
     ExpectedResult(cells=459, cells_with_inclusions=144, tolerance=0.1)
    ),  # passing
    (TestArgs(filename=LOCAL_RESOURCES_DIR / "cell_image_002.tif", output_folder=LOCAL_RESOURCES_DIR,
              clip_to=(CLIP_LOW, CLIP_HIGH), red_threshold=RED_THRESH, pixel_cycles=CYCLES, dead_cells=70),
     ExpectedResult(cells=680, cells_with_inclusions=421, tolerance=0.3)
    ),  # failing
    (TestArgs(filename=LOCAL_RESOURCES_DIR / "cell_image_003.tif", output_folder=LOCAL_RESOURCES_DIR,
              clip_to=(CLIP_LOW, CLIP_HIGH), red_threshold=RED_THRESH, pixel_cycles=CYCLES, dead_cells=70),
     ExpectedResult(cells=600, cells_with_inclusions=451, tolerance=0.1)
    ),
    (TestArgs(filename=LOCAL_RESOURCES_DIR / "cell_image_004.tif", output_folder=LOCAL_RESOURCES_DIR,
              clip_to=(CLIP_LOW, CLIP_HIGH), red_threshold=RED_THRESH, pixel_cycles=CYCLES, dead_cells=70),
     ExpectedResult(cells=502, cells_with_inclusions=231, tolerance=0.1)
    ),  # passing
    (TestArgs(filename=LOCAL_RESOURCES_DIR / "cell_image_005.tif", output_folder=LOCAL_RESOURCES_DIR,
              clip_to=(CLIP_LOW, CLIP_HIGH), red_threshold=RED_THRESH, pixel_cycles=CYCLES, dead_cells=70),
     ExpectedResult(cells=512, cells_with_inclusions=212, tolerance=0.3)
    ),
])
def test_create_mask(input, expected):
    ans = main_with_args(input)

    # assert( abs(ans["cells_after_clipping"] - expected.cells) < expected.tolerance*expected.cells )
    # assert (abs(ans["cells_with_inclusions"] - expected.cells_with_inclusions) < expected.tolerance * expected.cells_with_inclusions)
    ratio = ans["cells_with_inclusions"] / float(ans["cells_after_clipping"])
    assert(   abs(ratio - expected.ratio) < expected.tolerance*expected.ratio )
