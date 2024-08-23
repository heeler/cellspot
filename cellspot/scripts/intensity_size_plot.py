from cellspot.spot_finder import SpotFinder

import matplotlib.pyplot as plt
import numpy as np


def main():
    sf = SpotFinder(
        image_path='/Users/heeler/Sandbox/Python/cellspot/tests/data/cell_image_001.tif'
    )
    mask = sf.create_cell_mask()
    keepers = sf.keep_intersection()
    print(f"masked: {len(np.unique(mask))}, keepers: {len(keepers)}")
    # count_intensity = SpotFinder._calculate_count_intensities(nuclei_image=sf.image[:,:,2], cell_mask=mask)
    # count_max_intensity = SpotFinder._calculate_max_intensity(nuclei_image=sf.image[:, :, 2], cell_mask=mask)
    # spot_average_intensity = np.divide(count_intensity[:,1], count_intensity[:,1], dtype=np.float64)
    # ave_intensity = SpotFinder._create_average_intensities(nuclei_image=sf.image[:,:,2], cell_mask=mask)
    #plt.hist2d(x=count_intensity[1:, 0], y=spot_average_intensity[1:], bins=[100, 20])
    median_intensity = SpotFinder._calculate_median_intensity(nuclei_image=sf.image[:, :, 2], cell_mask=mask)
    plt.hist(median_intensity, bins=100)
    plt.show()
    # kept = sf.remove_empty_cells()
    # tifffile.imwrite('/Users/heeler/Sandbox/Python/cellspot/tests/data/kept_001.tiff', kept)


if __name__ == "__main__":
    main()