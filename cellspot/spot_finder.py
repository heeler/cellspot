import tifffile
from bioio import BioImage
from pathlib import Path
from cellpose import models
import numpy as np
import matplotlib.pyplot as plt
import cv2
import json


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
        self.results = {}


    @property
    def image(self) -> np.ndarray:
        return self._image


    def run(self, clip_values: tuple[float, float]=(5.0, 90), red_threshold: int=0, pixel_cycles: int=12) -> np.ndarray:
        self.create_cell_mask(clip_values=clip_values)
        # self.remove_empty_cells()
        self.remove_dead_cells_by_size(threshold=70) #70
        self.remove_empty_cells(within_pixels=pixel_cycles, red_threshold=red_threshold)
        return self.mask


    def create_cell_mask(self, clip_values: tuple[float, float]) -> np.ndarray:
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
            channels=[3, 3], # segment using only the blue channel as requested
            diameter=25.0,
            flow_threshold=4.0,
            normalize={"lowhigh": list(clip_values)},
            do_3D=False
        )
        mask = masks[0]
        self.mask = mask
        self.results['cells_by_cellpose'] = int(np.max(mask))
        return mask

    def keep_intersection(self, red_threshold: int=0) -> np.ndarray:
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
                if i > red_threshold and m > 0:
                    keepers.append(m)
        keepers = np.unique(keepers)
        self.keepers = keepers
        return keepers

    def remove_empty_cells(self, within_pixels: int = 0, red_threshold: int = 0) -> np.ndarray:
        if within_pixels <= 0:
            return self._remove_empty_cells(red_threshold=red_threshold)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
        dilate = cv2.dilate(self.mask, kernel, iterations=within_pixels)
        self.mask = dilate
        return self._remove_empty_cells(red_threshold=red_threshold) # cleanup where masks overlapped

    def _remove_empty_cells(self, red_threshold: int) -> np.ndarray:
        """
        This function creates a new mask consisting of
        only those cell masks found by keep_intersection().

        Returns
        -------
        A cell mask consisting of cell masks of cells
        with inclusion bodies.
        """
        if self.keepers is None:
            self.keep_intersection( red_threshold=red_threshold )

        with np.nditer([self.mask, None]) as it:
            for x, y in it:
                ans = 0
                if x in self.keepers:
                    ans = x
                y[...] = ans
            self.mask = it.operands[1]
            self.results['cells_with_inclusions'] = int(len(np.unique(self.mask)))
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
            self.results['cells_after_clipping'] = int(len(np.unique(self.mask)))
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

    def write_results(self, file_out:Path) -> dict:
        with file_out.open('w') as out:
            json.dump(self.results, out)
        return self.results

    def write_mask(self, filename: Path) -> None:
        tifffile.imwrite(str(filename), self.mask)
        file_npy = filename.with_suffix('.npy')
        self.mask.tofile(file_npy)






