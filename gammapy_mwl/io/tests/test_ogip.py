import pytest
import warnings
from gammapy_mwl.io.ogip import StandardOGIPDatasetReader
from gammapy.datasets import SpectrumDatasetOnOff
from astropy.io import fits
from numpy.testing import assert_allclose

@pytest.mark.parametrize(
    "pha_path",
    [
        "data/swift/xrt/swift_xrt_00036384074_src.pha",
        "data/swift/uvot/sw00036384074uvv_sk.src.pha",
        "data/swift/uvot/sw00036384074ubb_sk.src.pha",
        "data/swift/uvot/sw00036384074uuu_sk.src.pha",
        "data/swift/uvot/sw00036384074uw2_sk.src.pha",
        "data/swift/uvot/sw00036384074uw1_sk.src.pha",
        "data/swift/uvot/sw00036384074um2_sk.src.pha",
    ]
)
def test_standard_ogip_dataset_reader_runs(pha_path):
    """Test that StandardOGIPDatasetReader returns a SpectrumDatasetOnOff for a Swift/XRT or UVOT OGIP PHA file,
    and that OGIP keywords are read correctly."""
    # Suppress Astropy warnings for test
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fits.open(pha_path) as hdulist:
            header = hdulist[1].header  # Usually the spectrum is in extension 1
            ogip_class = header.get("HDUCLASS")
            ogip_clas1 = header.get("HDUCLAS1")
            exposure = header.get("EXPOSURE")
            assert ogip_class == "OGIP"
            assert ogip_clas1 == "SPECTRUM"
            assert exposure is not None

        reader = StandardOGIPDatasetReader(pha_path)
        dataset = reader.read()
        assert isinstance(dataset, SpectrumDatasetOnOff)
        assert hasattr(dataset, "counts")
        assert hasattr(dataset, "acceptance")
        assert hasattr(dataset, "counts_off")
        assert hasattr(dataset, "acceptance_off")
        assert hasattr(dataset, "edisp")
        assert hasattr(dataset, "exposure")

        # Check that the exposure and livetime are read correctly from the dataset metadata
        meta = dataset.meta_table
        assert_allclose(meta["EXPOSURE"], exposure, rtol=1e-3)
