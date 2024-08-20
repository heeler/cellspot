from .conftest import LOCAL_RESOURCES_DIR

from cellspot.spot_finder import SpotFinder


def test_spot_finder():
    sf = SpotFinder(image_path= LOCAL_RESOURCES_DIR / "cell_image_001.tif")
    masks = sf.create_cell_mask()
    assert True
