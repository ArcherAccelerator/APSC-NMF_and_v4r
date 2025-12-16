# APSC-NMF_and_v4r
MATLAB toolkit for adaptive pollution source-constrained nonnegative matrix factorization (APSC-NMF) source apportionment of heavy metals in PM2.5, with spatial interpolation (OK, IDW) and uncertainty analysis. Ideal for environmental science research. Includes constrained optimization and cross-validation.
The code is derived from methods used in a peer-reviewed paper and has been refactored for clarity, modularity, and reproducibility.

## Table of Contents
- [Scripts Overview](#scripts-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Testing](#testing)
- [Deployment](#deployment)
- [Contributing](#contributing)
- [License](#license)

## Scripts Overview

### 1. NMF_Source_Apportionment.m
**Purpose**: Implements Adaptive Pollution Source-Constrained Nonnegative Matrix Factorization (APSC-NMF) for identifying pollution sources in heavy metal concentrations, with sequential constraints on factor loadings, min-max normalization, and direct comparison to PMF results for model validation.

**Input Data**:
- `SourceApportionmentData.xls`: Excel file (95 rows/samples × 17 columns; columns 1-4: metadata, 5-17: heavy metal concentrations).
- `PMF_G_Matrix.xlsx` and `PMF_F_Matrix.xlsx`: PMF-derived factor scores and loadings for benchmark comparison.
- `H_Sequence_Limit.csv`: CSV constraints (first column: source names; remaining columns: 0/1/NaN sequence matrix for factor ordering).

**Output Results**:
- `NMF_Source_Apportionment_Results.mat`: Cell array containing source combinations, factor scores/loadings matrices, reconstruction errors (normalized and denormalized), constraint success flags, and R²/RMSE tables comparing NMF vs. PMF.
- Console output: Iteration progress, best factorization summaries, and validation diagnostics.

### 2. Spatial_Interpolation_Analysis.m
**Purpose**: Performs comparative spatial interpolation of heavy metal concentrations using Ordinary Kriging (OK), Inverse Distance Weighting (IDW), linear v4, and refined v4r methods, with random subsampling cross-validation and study area masking for accurate mapping.

**Input Data**:
- `SamplingPoint_HeavyMetalsConcentrations.xlsx`: Excel file (rows: sampling points; columns 1-2: longitude/latitude coordinates, 3+: concentrations for 12 heavy metals).
- `SamplingDataSeparatedByInteriorPoints.xlsx`: Classifications (column 1: point types, e.g., interior vs. boundary).
- `Q4.shp`: Shapefile defining the study area boundary (includes .shx and .dbf for polygon clipping).

**Output Results**:
- PNG visualization files in `outputDir` (e.g., `Cd_n5_sample1.png`): 2×2 subplot figures per metal/sample size/run, displaying interpolated grids, color-scaled maps, training/test points, min/max values, and colorbars.
- Console output: Processing logs (e.g., sample size, metal, run progress).

### 3. NMF_Uncertainty_Analysis.m
**Purpose**: Quantifies uncertainty in APSC-NMF results through Monte Carlo simulation (repeated factorizations), calculating mean and standard deviation of factor scores/loadings to assess robustness of source apportionment estimates.

**Input Data**:
- `SourceApportionmentData.xls`: Same concentration data as above (95 samples × heavy metals).
- `H_Sequence_Limit.csv`: Same sequence constraints for consistent factorization.

**Output Results**:
- `NMF_Uncertainty_Metrics_5.mat`: Matrices for mean/std deviation of factor scores and loadings, vector of denormalized errors across runs, and cell array of per-run results (including flags and errors).
- Console output: Run progress and best replicate summaries.

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
   ```
   git clone https://github.com/yourusername/APSC-NMF_and_v4r.git
   cd APSC-NMF_and_v4r
   ```

2. **Set Up MATLAB Path**:
   - Open MATLAB and add the repository root to your path:
     ```matlab
     addpath(genpath(pwd));
     ```
   - For Kriging (in Spatial_Interpolation_Analysis.m): Download and add ooDACE-1.4:
     ```matlab
     addpath('path/to/ooDACE-1.4/ooDACE');
     ```

3. **Prepare Data**:
   - Place input files (e.g., `.xls`, `.csv`, `.shp`) in the paths defined in each script (e.g., under "Configuration Parameters"). Edit these paths for your local setup.
   - Ensure shapefiles are complete (`.shp`, `.shx`, `.dbf`) and in a compatible projection (e.g., WGS84).

4. **Verify Installation**:
   - Run `clear; clc; verifySetup;` in MATLAB (create a simple `verifySetup.m` script to check paths and toolboxes, e.g., `ver` command for toolboxes).

## Usage
### Running the Scripts
1. **NMF_Source_Apportionment.m** (Base Analysis):
   - Edit parameters (e.g., `numFactors = 4`).
   - Run: `NMF_Source_Apportionment`.
   - Use outputs for downstream scripts.

2. **Spatial_Interpolation_Analysis.m** (Mapping):
   - Set `sampleSizes = [5, 10, 15, 20]` and methods.
   - Run: `Spatial_Interpolation_Analysis`.
   - Views: Generated PNGs for visual comparison.

3. **NMF_Uncertainty_Analysis.m** (Robustness Check):
   - Set `selectedSources = [1, 2, 4, 5]`.
   - Run: `NMF_Uncertainty_Analysis` (after base NMF).

### Example Workflow
```matlab
% In MATLAB Command Window
NMF_Source_Apportionment;  % Compute sources
NMF_Uncertainty_Analysis;  % Assess variability
Spatial_Interpolation_Analysis;  % Map concentrations
```
Customize via script parameters (e.g., `maxIterations = 1000`).

## Testing
### Unit Tests
- Add to `/tests/` folder (e.g., `test_nmf_helpers.m`):
  ```matlab
  % test_nmf_helpers.m
  testData = rand(10, 5);  % Synthetic matrix
  [normData, minV, maxV] = minMaxNormalization(testData(:,1));
  assert(all(normData >= 0 & normData <= 1), 'Normalization failed');
  fprintf('All tests passed.\n');
  ```
- Run: `runtests('tests/')` (requires MATLAB Test Framework).

### Integration Tests
1. **Synthetic Data Test**:
   - Generate sample inputs: `concentrationMatrix = rand(95, 13) * 10; writematrix(concentrationMatrix, 'test_data.xls');`.
   - Run scripts and verify: R² > 0.8, errors < 1e-2, figures without NaNs.

2. **Expected Outputs**:
   - NMF: Factor matrices with non-negative values; validation flags = 0.
   - Spatial: 2×2 PNGs per run with clipped masks.
   - Uncertainty: Std devs < 10% of means.

3. **Debugging**:
   - Set `'Display', 'iter'` in `statset` for verbose logs.
   - Check for warnings (e.g., constraint failures).

## Deployment
### Local Deployment
- Scripts are standalone; execute directly in MATLAB.
- For automation: Wrap in functions and use `batch` jobs.

### Server/Cloud Deployment
1. **MATLAB Production Server**:
   - Compile: `mcc -m NMF_Source_Apportionment.m`.
   - Deploy as API endpoint for input data processing.

2. **Containerization (Docker)**:
   - Sample `Dockerfile`:
     ```
     FROM mathworks/matlab:r2024a
     COPY . /app
     WORKDIR /app
     RUN matlab -batch "addpath(genpath(pwd)); run('NMF_Source_Apportionment.m');"
     CMD ["matlab", "-batch", "run('main_workflow.m');"]
     ```
   - Build/Run: `docker build -t apsc-nmf . && docker run -v /local/data:/app/data apsc-nmf`.

3. **Cloud (e.g., AWS Batch)**:
   - Upload to S3/MATLAB Drive.
   - Schedule via Lambda: Pass parameters (e.g., numFactors) and retrieve .mat outputs.

### Portability Notes
- Use `fullfile` for paths; tested on Windows/Linux.
- No external APIs; fully offline.

## Contributing
- Fork the repo and submit pull requests for enhancements (e.g., additional interpolation methods like RBF).
- Style: Descriptive variable names, `%%` section headers, vectorized code.
- Issues: Report bugs with sample data snippets.

## License
MIT License. See [LICENSE](LICENSE) for details. For academic use, cite the associated peer-reviewed paper: "Author et al. (2024). Adaptive Pollution Source-Constrained NMF for PM2.5 Heavy Metal Apportionment. *Journal of Environmental Science*."

For support: [li.air@qq.coml. Last updated: December 16, 2025.
