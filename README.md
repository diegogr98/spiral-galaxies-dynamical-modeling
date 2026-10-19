# N-Body to IFU: Synthetic Kinematics Pipeline for Spiral Galaxies

This repository contains the codebase and data processing pipeline for my astrophysics thesis. The project bridges theoretical galactic dynamics and observational astronomy by generating synthetic Integral Field Unit (IFU) kinematic maps from 3D N-body simulations and systematically comparing them with real galaxy observations from the CALIFA survey.

## Scientific Overview

The core objective of this research is to reproduce the observed kinematic properties (line-of-sight velocity $V_{\rm los}$ and velocity dispersion $\sigma_{\rm los}$) of spiral galaxies. By simulating exact analogs of specific CALIFA targets, this pipeline allows us to constrain unobservable parameters, particularly the dark matter halo mass ($M_{200}$), through dynamic stability and kinematic matching.

### Key Features
* **Automated Parameter Extraction:** Direct ingestion of observational parameters (Position Angle, Inclination, Redshift/Distance, Stellar Mass) from CALIFA FITS and CSV catalogs.
* **Initial Conditions Generation:** Integration with `MakeDisk` to generate stable multi-component (Halo, Bulge, Disk, Gas) N-body initial conditions.
* **3D-to-2D Synthetic Projection:** Geometric rotation and projection of 3D N-body particles into a 2D observer frame using small-angle approximation.
* **IFU Spectrograph Simulation:** Spatial binning of particle velocities into spaxels, mimicking PyCASSO Voronoi-tessellated kinematic maps, including Signal-to-Noise (S/N) masking based on particle count.

## Repository Structure

```text
├── data/
│   ├── raw/                 # Raw FITS catalogs (e.g., gme.fits) and simulation outputs (.dat, .giz)
│   └── processed/           # Filtered parameters and clean datasets (parametros_con_gas.csv)
├── src/
│   ├── analysis/
│   │   ├── momentos_proy.py # Core script: 3D rotation, 2D projection, and IFU map generation
│   │   └── ...              # Additional diagnostic and analysis scripts
│   └── utils/
│       └── geometria.py     # Euler rotation matrices and spatial transformations
├── outputs/
│   └── figures/             # Generated synthetic kinematic maps (PNG/PDF)
└── run_pipeline.py          # Main execution script
