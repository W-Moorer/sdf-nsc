# Reproduce Historical `gear_error_render`

The checked-in `papers/figures/gear_error_render.png` was produced by this workflow:

1. Run the raw 1st/2nd simple-gear simulations.
2. Copy the raw outputs to:
   - `data/outputs/gear_tricubic_1st.csv`
   - `data/outputs/gear_tricubic_2nd.csv`
3. Apply the historical MAE normalization in `scripts/fix_zero.py`.
4. Plot the normalized tricubic files against the commercial reference.

Run the full chain from the repository root with:

```bash
bash scripts/run_reproduce_gear.sh
```

This regenerates:

- raw outputs:
  - `data/outputs/gear_1st_highres_0001.csv`
  - `data/outputs/gear_2nd_highres_0001.csv`
- normalized historical figure inputs:
  - `data/outputs/gear_tricubic_1st.csv`
  - `data/outputs/gear_tricubic_2nd.csv`
- figures:
  - `papers/figures/gear_error.pdf`
  - `papers/figures/gear_error_render.png`
