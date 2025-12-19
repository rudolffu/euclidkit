# Euclid DR1 MER MOC products

This directory contains FITS MOCs (Multi-Order Coverage maps) describing the sky footprint of Euclid DR1 MER catalogues. All three products are:

- **MOC version**: 2.0 (FITS, `ORDERING=RANGE`, `COORDSYS=C` / ICRS)
- **Resolution**: `o13` (`max_order = 13`)

## Files

Approximate sky areas below are computed as `moc.sky_fraction * 41252.96` (deg²).

| File | What it represents | Area (deg²) |
|---|---|---:|
| `dr1_mer_wide_bins_o13_moc.fits` | Wide **catalogue** footprint. | 2014.57 |
| `dr1_mer_deep_o13_moc.fits` | Deep **catalogue** footprint. | 77.78 |
| `dr1_mer_wide_deep_union_o13_moc.fits` | Union of the wide + deep catalogue MOCs (combined footprint). | 2058.33 |

## Quick usage (mocpy)

```python
from mocpy import MOC

moc = MOC.from_fits("mochp/dr1_mer_wide_bins_o13_moc.fits")
print("max_order:", moc.max_order)
print("area (deg^2):", moc.sky_fraction * 41252.96)
```
