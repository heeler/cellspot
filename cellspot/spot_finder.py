import tifffile
from bioio import BioImage
from pathlib import Path
from cellpose import models
import numpy as np
import matplotlib.pyplot as plt
import cv2


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
        self._image = tifffile.imread(self.img_path)
        self.model = models.Cellpose(gpu=False, model_type='cyto3')
        self.mask = None
        self.keepers = None


    @property
    def image(self) -> np.ndarray:
        return self._image


    def run(self) -> np.ndarray:
        self.create_cell_mask()
        # self.remove_empty_cells()
        self.remove_dead_cells_by_size(threshold=150)
        self.remove_empty_cells(within_pixels=12)
        return self.mask


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
            diameter=25.0,
            flow_threshold=4.0,
            normalize={"lowhigh": [5.0, 90.0]}, # clipping the top 1.5 % to remove dead nuclei
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

    def remove_empty_cells(self, within_pixels: int = 0) -> np.ndarray:
        if within_pixels <= 0:
            return self._remove_empty_cells()
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
        dilate = cv2.dilate(self.mask, kernel, iterations=within_pixels)
        self.mask = dilate
        return self._remove_empty_cells() # cleanup where masks overlapped

    def _remove_empty_cells(self) -> np.ndarray:
        """
        This function creates a new mask consisting of
        only those cell masks found by keep_intersection().

        Returns
        -------
        A cell mask consisting of cell masks of cells
        with inclusion bodies.
        """
        if self.keepers is None:
            self.keep_intersection()

        with np.nditer([self.mask, None]) as it:
            for x, y in it:
                ans = 0
                if x in self.keepers:
                    ans = x
                y[...] = ans
            self.mask = it.operands[1]
            return self.mask

    def remove_dead_cells_by_size(self, threshold: int = 100) -> np.ndarray:
        count_intensities = self._calculate_count_intensities(nuclei_image=self.image[:,:,2], cell_mask=self.mask)
        keepers = []
        for index, x in np.ndenumerate([count_intensities[:, 0]]):
            if x > threshold and x != 0:
                keepers.append(index[1])
        with np.nditer([self.mask, None]) as it:
            for x, y in it:
                ans = 0
                if x in keepers:
                    ans = x
                y[...] = ans
            self.mask = it.operands[1]
            return self.mask



    @staticmethod
    def _calculate_count_intensities(nuclei_image: np.ndarray, cell_mask: np.ndarray) -> np.ndarray:
        max_index = np.max(cell_mask) + 1
        count_intensity = np.zeros((max_index, 2), dtype=np.int32)
        with np.nditer([nuclei_image, cell_mask]) as it:
            for x, y in it:
                count_intensity[y, 0] += 1
                count_intensity[y, 1] += x
        return count_intensity

    @staticmethod
    def _calculate_max_intensity(nuclei_image: np.ndarray, cell_mask: np.ndarray) -> np.ndarray:
        max_index = np.max(cell_mask) + 1
        count_intensity = np.zeros((max_index, 2), dtype=np.int32)
        with np.nditer([nuclei_image, cell_mask]) as it:
            for x, y in it:
                count_intensity[y, 0] += 1
                if x > count_intensity[y, 1]:
                    count_intensity[y, 1] = x
        return count_intensity

    @staticmethod
    def _create_average_intensities(nuclei_image: np.ndarray, cell_mask: np.ndarray) -> np.ndarray:
        count_intensity = SpotFinder._calculate_count_intensities(nuclei_image, cell_mask)
        return np.divide(count_intensity[:, 1], count_intensity[:, 0], dtype=float)

    @staticmethod
    def _calculate_median_intensity(nuclei_image: np.ndarray, cell_mask: np.ndarray) -> np.ndarray:
        max_index = np.max(cell_mask) + 1
        data_by_index = [[] for _ in range(max_index)]
        with np.nditer([nuclei_image, cell_mask]) as it:
            for x, y in it:
                if x > 0 and y > 0:
                    data_by_index[y].append(x)
        data_by_index[0].append(0)
        return np.array([np.percentile(x, 98) for x in data_by_index])

    def write_mask(self, filename: Path) -> None:
        tifffile.imwrite(str(filename), self.mask)
        file_npy = filename.with_suffix('.npy')
        self.mask.tofile(file_npy)






