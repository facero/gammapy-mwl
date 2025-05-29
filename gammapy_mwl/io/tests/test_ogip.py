import pytest
from gammapy_mwl.io.ogip import StandardOGIPDatasetReader
from gammapy.datasets import SpectrumDatasetOnOff
from astropy.io import fits
from numpy.testing import assert_allclose

@pytest.mark.parametrize(
    "pha_path",
    ["data/swift/xrt/swift_xrt_00036384074_src.pha"]
)
def test_standard_ogip_dataset_reader_runs(pha_path):
    """Test that StandardOGIPDatasetReader returns a SpectrumDatasetOnOff for a Swift/XRT OGIP PHA file,
    and that OGIP keywords are read correctly."""
    # Check OGIP keywords with astropy.io.fits
    with fits.open(pha_path) as hdulist:
        header = hdulist[1].header  # Usually the spectrum is in extension 1
        ogip_class = header.get("HDUCLASS")
        ogip_clas1 = header.get("HDUCLAS1")
        exposure = header.get("EXPOSURE")
        livetime = header.get("LIVETIME")
        assert ogip_class == "OGIP"
        assert ogip_clas1 == "SPECTRUM"
        assert exposure is not None
        assert livetime is not None

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
    assert_allclose(meta["LIVETIME"], livetime, rtol=1e-3)
