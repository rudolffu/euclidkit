# Euclid DR1 MER MOC products

This directory contains FITS MOCs (Multi-Order Coverage maps) describing Euclid DR1 sky footprints. Available products are:

- **MOC version**: 2.0 (FITS, `ORDERING=RANGE`, `COORDSYS=C` / ICRS)
- **Resolution**: `o13` (`max_order = 13`)

## Files

Approximate sky areas below are computed as `moc.sky_fraction * 41252.96` (deg²).

| File | What it represents | Area (deg²) |
|---|---|---:|
| `dr1_mer_wide_bins_o13_moc.fits` | Wide **catalogue** footprint. | 2014.57 |
| `dr1_mer_deep_o13_moc.fits` | Deep **catalogue** footprint. | 77.78 |
| `dr1_mer_wide_deep_union_o13_moc.fits` | Union of the wide + deep catalogue MOCs (combined footprint). | 2058.33 |
| `cgv_map_dr1input_o13_moc.fits` | CGV input footprint map. | 2108.51 |

## Quick usage (mocpy)

```python
from importlib import resources
from mocpy import MOC

resource = resources.files("euclidkit").joinpath("data", "mocs", "dr1_mer_wide_bins_o13_moc.fits")
with resources.as_file(resource) as moc_path:
    moc = MOC.from_fits(moc_path)
print("max_order:", moc.max_order)
print("area (deg^2):", moc.sky_fraction * 41252.96)
```
