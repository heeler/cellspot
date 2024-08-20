import tifffile
from bioio import BioImage
from pathlib import Path
from cellpose import models
import numpy as np


class SpotFinder(object):
    def __init__(self, image_path: Path) -> None:
        """
        Load the Tif file provided in image_path and set the
        class internals to use the cyto3 model for image segmentation.
        Parameters
        ----------
        image_path: Path to a Tif file with 3 channels / RGB in our case.

        Returns
        -------
        None
        """
        self.img_path = image_path
        self.image = tifffile.imread(self.img_path)
        self.model = models.Cellpose(gpu=False, model_type='cyto3')
        self.mask = None
        self.keepers = None

    def create_cell_mask(self) -> np.ndarray:
        """
        This launches CellPose and captures the outputs.
        The relevant output currently is the mask.

        Returns
        -------
        The mask of the detected cells. Each cell found is
        given an index and the pixels for that cell mask have
        said integer value.
        """
        masks, flows, styles, diams = self.model.eval(
            [self.image],
            channels=[3, 3],
            diameter=12.7,
            flow_threshold=4.0,
            do_3D=False
        )
        mask = masks[0]
        print(f"mask shape: {mask.shape}, mask.max: {np.max(mask)}")
        self.mask = mask
        return mask

    def keep_intersection(self) -> np.ndarray:
        """
        This function creates a list of the mask integers which
        overlap a red pixel, aka an inclusion body.

        Returns
        -------
        A unique list of mask integers.
        """
        keepers = []
        with np.nditer([self.image[:,:,0], self.mask]) as it:
            for i, m in it:
                if i > 0 and m > 0:
                    keepers.append(m)
        keepers = np.unique(keepers)
        self.keepers = keepers
        return keepers

    def remove_empty_cells(self) -> np.ndarray:
        """
        This function creates a new mask consisting of
        only those cell masks found by keep_intersection().

        Returns
        -------
        A cell mask consisting of cell masks of cells
        with inclusion bodies.
        """
        with np.nditer([self.mask, None]) as it:
            for x, y in it:
                ans = 0
                if x in self.keepers:
                    ans = x
                y[...] = ans
            return it.operands[1]




if __name__ == '__main__':
    sf = SpotFinder(
        image_path='/Users/heeler/Sandbox/Python/cellspot/tests/data/cell_image_001.tif'
    )
    mask = sf.create_cell_mask()
    keepers = sf.keep_intersection()
    print(f"masked: {len(np.unique(mask))}, keepers: {len(keepers)}")
    kept = sf.remove_empty_cells()
    tifffile.imwrite('/Users/heeler/Sandbox/Python/cellspot/tests/data/kept_001.tiff', kept)
