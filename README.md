# APSC-NMF_and_v4r
MATLAB toolkit for adaptive pollution source-constrained nonnegative matrix factorization (APSC-NMF) source apportionment of heavy metals in PM2.5, with spatial interpolation (OK, IDW) and uncertainty analysis. Ideal for environmental science research. Includes constrained optimization and cross-validation.

The code is derived from methods used in a peer-reviewed paper and has been refactored for clarity, modularity, and reproducibility.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Testing](#testing)
- [Deployment](#deployment)
- [Contributing](#contributing)
- [License](#license)

## Requirements
### Software
- **MATLAB**: Version R2020b or later (tested on R2024a).
- **Toolboxes**:
  - Statistics and Machine Learning Toolbox (for `nnmf` and optimization functions like `fmincon`).
  - Mapping Toolbox (for shapefile reading via `shaperead` and `inpolygon`).
  - Optimization Toolbox (for constrained linear programming in NMF).
  - Optional: ooDACE Toolbox (for Kriging; download from [GitHub](https://github.com/edwin-de-jong/ooDACE) and add to path).

### Data
- Input files (e.g., Excel sheets for concentrations, shapefiles for boundaries) must be provided in the specified paths. Sample data structure:
  - `SourceApportionmentData.xls`: Columns 1-4 metadata, 5-17 heavy metal concentrations (95 rows).
  - `H_Sequence_Limit.csv`: Source contribution sequence constraints.
  - `SamplingPoint_HeavyMetalsConcentrations.xlsx`: Columns 1-2 coordinates, 3+ concentrations.
- No internet access required; all computations are local.

### Hardware
- Minimum: 8 GB RAM, multi-core CPU (for replicates).
- Recommended: 16 GB RAM for uncertainty analysis (200+ replicates).

## Installation
1. **Clone the Repository**:
