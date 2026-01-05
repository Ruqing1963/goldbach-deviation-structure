# Fractal-Spectral Structure of Goldbach Deviations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

This repository contains the data, code, and analysis for the paper:

**"Fractal-Spectral Structure of Goldbach Deviations: Long-Range Persistence and Amplitude Scaling"**

by Ruqing Chen (GUT Geoservice Inc., Montreal)

## Key Findings

1. **Long-Range Persistence**: Hurst exponent H = 0.84 ≫ 0.5
2. **Self-Similar Amplitude Scaling**: κ ≈ C₂ = 0.6602 (twin prime constant)
3. **Spectral Signature**: FFT peaks match L-function zeros

## Repository Structure

```
goldbach-deviation-structure/
├── README.md
├── LICENSE
├── paper/
│   ├── paper_goldbach_deviation.pdf
│   └── paper_goldbach_deviation.tex
├── data/
│   ├── goldbach_counts.csv
│   └── hurst_analysis.csv
├── figures/
│   ├── fig1_scale_separation.png
│   ├── fig2_hurst_analysis.png
│   ├── fig3_fft_spectrum.png
│   └── fig4_safety_gap.png
├── code/
│   ├── goldbach_analysis.py
│   └── requirements.txt
└── scripts/
    └── run_analysis.sh
```

## Quick Start

```bash
# Clone repository
git clone https://github.com/Ruqing1963/goldbach-deviation-structure.git
cd goldbach-deviation-structure

# Install dependencies
pip install -r code/requirements.txt

# Run analysis
python code/goldbach_analysis.py
```

## Data Description

- `goldbach_counts.csv`: Goldbach representation counts G(N) for even N from 4 to 500,000
- `hurst_analysis.csv`: R/S analysis results for Hurst exponent calculation

## Citation

If you use this code or data, please cite:

```bibtex
@article{chen2026goldbach,
  title={Fractal-Spectral Structure of Goldbach Deviations: Long-Range Persistence and Amplitude Scaling},
  author={Chen, Ruqing},
  year={2026},
  publisher={Zenodo},
  doi={10.5281/zenodo.XXXXXXX}
}
```

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

- Author: Ruqing Chen
- Email: ruqing@hotmail.com
- Affiliation: GUT Geoservice Inc., Montreal
