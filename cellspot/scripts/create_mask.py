from fileinput import filename

from cellspot.spot_finder import SpotFinder

from pathlib import Path
import tifffile


def main():
    sf = SpotFinder(
        image_path='/Users/heeler/Sandbox/Python/cellspot/tests/data/cell_image_001.tif'
    )
    mask = sf.run()
    sf.write_mask(filename=Path('/Users/heeler/Sandbox/Python/cellspot/tests/data/mask_001.tif'))


if __name__ == '__main__':
    main()
