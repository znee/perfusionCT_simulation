# Perfusion CT Simulation

Interactive educational tool for visualizing CT perfusion physics and hemodynamic parameters.

## Features

- **Real-time Visualization**: Interactive sliders for CBF, MTT, Tmax, cardiac output, and noise
- **Physiologically Accurate**: Implements proper convolution mathematics with realistic HU values
- **SVD Deconvolution**: Recovers perfusion parameters from tissue curves using Tikhonov regularization
- **Multiple Models**: Boxcar (plug flow) and Exponential (perfect mixing) residue functions
- **Educational Annotations**: Visual markers for CBF, MTT, Tmax, and CBV on IRF chart
- **Noise Simulation**: Demonstrates deconvolution artifacts and regularization effects

## Quick Start

### Local Usage with serve.py

**Why use serve.py?**
The tool uses external resources (Chart.js, MathJax) that require a web server due to browser security (CORS) restrictions. Opening `index.html` directly in a browser won't work properly.

**What is serve.py?**
A simple Python HTTP server that:
- Serves the application on `localhost:8000`
- Uses ThreadingTCPServer for concurrent connections
- Automatically opens your default browser
- Provides read-only access (secure for presentations)

**Requirements:**
- Python 3.x (built-in, no additional packages needed)

**Run:**
```bash
python3 serve.py
```

The app will automatically open at `http://localhost:8000`

**Stopping the server:**
Press `Ctrl+C` in the terminal

### GitHub Pages
Visit the live demo at: **https://znee.github.io/perfusionCT_simulation/**

## Physiological Parameters

### CBF (Cerebral Blood Flow)
- Range: 10-100 ml/100g/min
- Normal gray matter: 60-80 ml/100g/min
- Normal white matter: 20-30 ml/100g/min

### MTT (Mean Transit Time)
- Range: 1-10 seconds
- Represents average time for blood to traverse tissue capillary bed

### Tmax (Delay)
- Range: 0-15 seconds
- Arterial-tissue delay (time from AIF peak to tissue peak)

### Cardiac Output
- Range: 2-8 L/min
- Affects AIF timing, width, and amplitude

## Mathematical Implementation

### Core Equation
```
Ct(t) = (CBF/60) × dt × [AIF(t) ⊗ R(t)]
```

Where:
- `Ct(t)`: Tissue concentration curve (HU)
- `CBF`: Cerebral blood flow (ml/100g/min) → divided by 60 for per-second units
- `dt`: Time resolution (1 second)
- `AIF(t)`: Arterial Input Function (gamma-variate, ~200 HU peak)
- `R(t)`: Residue function, **shifted by Tmax** for forward model
- `⊗`: Convolution operator

**Important**: Tmax is applied to the residue function, not the AIF. This ensures the forward and inverse problems represent the same physical system.

### Tissue Concentration Scaling Factor (0.1)

The tissue concentration includes a **0.1 scaling factor** for physiological accuracy:

#### Rationale

1. **Extraction Fraction (E)**: Only ~30-50% of contrast agent crosses from blood into tissue
   - Small molecules like iodinated contrast have incomplete extraction
   - Brain tissue extraction: typically 40-50%

2. **Hematocrit Correction**: 
   - AIF is measured in whole blood (contains ~45% red blood cells)
   - Tissue enhancement comes from plasma only (~55% of blood volume)
   - Correction factor: ~0.55

3. **Partial Volume Effects**:
   - Tissue voxels contain mixture of capillaries, cells, and extracellular space
   - Not all tissue volume enhances uniformly
   - Effective reduction: ~0.45

#### Mathematical Justification
```
Combined scaling = Plasma fraction × Extraction × Tissue factors
                 = 0.55 × 0.4 × 0.45
                 ≈ 0.1
```

This brings tissue concentration to the **realistic 20-60 HU range** for typical CBF values (30-80 ml/100g/min).

### Other Key Features

**Canonical Residue Function**: R(0) = 1, no intrinsic delay
- **Tmax applied to residue**: Forward model uses `R(t - Tmax)` to generate tissue curve
- This ensures the matrix `A` (built from unshifted AIF) and the tissue curve `b` represent the same system
- Educational: "Ideal IRF" shows canonical response, "Shifted IRF" shows delayed tissue response

**Fractional Tmax Support**: Linear interpolation for sub-second delays (e.g., 0.5s)

**Mass Conservation**: VOF normalized to preserve AIF area (within 1% error tolerance)

**CBV Verification**: Console logs numeric ∫IRF·dt vs theoretical (CBF/60)×MTT

## Model Types

### Boxcar (Plug Flow)
- Assumes all blood takes exactly MTT seconds to transit
- R(t) = 1 for 0 ≤ t < MTT, else 0
- Transport: Delta function at t = MTT

### Exponential (Perfect Mixing)
- Assumes instantaneous mixing in tissue compartment
- R(t) = exp(-t/MTT) for t ≥ 0
- Transport: (1/τ) × exp(-t/τ) where τ = MTT

## Deconvolution Feature

### Overview
The simulation implements **SVD-based deconvolution with Tikhonov regularization** to recover perfusion parameters from the tissue concentration curve, mirroring real-world clinical CT perfusion analysis.

### Algorithm: Block-Circulant SVD

**Forward Problem (Convolution)**:
```
C_t(t) = CBF · [AIF(t) ⊗ R(t)]
```

**Inverse Problem (Deconvolution)**:
```
R(t) = deconv(C_t(t), AIF(t))
```

**Implementation Steps**:
1. Build Toeplitz matrix A from AIF
2. Perform Singular Value Decomposition: A = U Σ V^T
3. Apply Tikhonov regularization to singular values:
   ```
   S_inv[i] = S[i] / (S[i]² + λ²·max(S)²)
   ```
4. Reconstruct solution: R = V · S_inv · U^T · b

### Deconvolution Modes

The tool offers two modes to demonstrate the impact of arterial delay (Tmax) on parameter estimation:

1.  **Delay-Sensitive (Standard SVD)**:
    - Uses a standard Toeplitz matrix constructed from the AIF.
    - **Behavior**: If Tmax > 0, the tissue curve is shifted relative to the AIF. Standard SVD interprets this shift as dispersion, leading to an **underestimation of CBF** and overestimation of MTT.
    - **Educational Value**: Demonstrates why delay correction is critical in perfusion analysis.

2.  **Delay-Insensitive (Block-Circulant)**:
    - Uses a block-circulant matrix with wrapped causal tail.
    - **Formula**: Lower triangle uses standard Toeplitz `A[i][j] = aif[i-j]`; upper triangle wraps `A[i][j] = aif[(i-j+N) % N]`
    - **Behavior**: Allows the IRF to shift beyond the acquisition window, correctly recovering CBF, CBV, Tmax, and MTT even when arterial delay is present.
    - **Educational Value**: Represents modern "delay-insensitive" or "oscillatory" SVD algorithms (oSVD) used in clinical software.

### Technical Implementation Notes

**Forward Model Consistency**:
- The forward problem generates tissue using `C_t = AIF ⊗ R_shifted`
- The inverse problem solves `A·x = C_t` where `A` is built from the **same unshifted AIF**
- This ensures mathematical consistency between forward and inverse operations

**Pure Tikhonov Regularization**:
- No hard truncation of singular values
- Damping formula: `S_inv[i] = S[i] / (S[i]² + (λ·S_max)²)`
- Default λ = 0.05 matches clinical perfusion software (typical range: 0.03-0.08)

**Noisy AIF for Deconvolution**:
- Both tissue curve AND AIF receive noise (when noise > 0%)
- The noisy AIF is used to build the deconvolution matrix
- Simulates real-world instability where both signals are corrupted

### Regularization Parameter (λ)

**Purpose**: Balances noise suppression vs. accuracy

- **λ = 0.0**: No regularization → sensitive to noise, oscillations
- **λ = 0.05** (default): Optimal for moderate noise
- **λ = 0.2**: Heavy smoothing → may lose detail

**Effect**: Higher λ truncates small singular values, reducing high-frequency artifacts but potentially underestimating peaks.

### Parameter Derivation

From deconvolved IRF, the following perfusion parameters are extracted:

1. **CBF (Cerebral Blood Flow)**:
   - Peak value of deconvolved IRF
   - Units: ml/100g/min

2. **Tmax (Time to Peak)**:
   - Time at which IRF reaches maximum
   - Units: seconds

3. **CBV (Cerebral Blood Volume)**:
   - Area under IRF curve
   - ∫IRF(t)·dt with unit conversion (divide by 60)
   - Units: ml/100g

4. **MTT (Mean Transit Time)**:
   - Central Volume Theorem: MTT = (CBV × 60) / CBF
   - Units: seconds

### Educational Value

The deconvolution feature demonstrates:
- **Ill-posed inverse problem**: Small noise → large oscillations in solution
- **Regularization trade-offs**: Smoothness vs. accuracy
- **Parameter estimation**: How clinical software derives CBF, CBV, MTT, Tmax
- **Real-world workflow**: Comparison between ground truth and recovered parameters

## Chart Details

### AIF Chart (0-70s)
- **Ideal AIF** (solid red): Gamma-variate, peak ~200 HU
- **Noisy AIF** (thin red): With added noise
- **VOF** (dashed teal): Venous output (mass-conserved)

### IRF Chart (0-70s)
- **Ideal IRF** (dashed light teal): Canonical CBF×R(t) at t=0
- **Shifted IRF** (solid dark teal): Response delayed by Tmax
- **Deconvolved IRF** (orange): Recovered from tissue curve via SVD deconvolution

### Tissue Chart (0-70s)
- **Tissue Curve** (yellow): Convolution result with realistic HU (20-60 range)

## Technical Notes

- **Time resolution**: 1 second (dt = 1)
- **Scan duration**: 60 seconds (70s display for better visualization)
- **AIF calibration**: Systematic scaling to 200 HU peak at normal cardiac output
- **Transport normalization**: h(t) integrates to 1; exponential uses (1/τ) factor
- **Auto y-axis scaling**: Charts adjust to data range automatically

## Educational Use

This tool is designed for teaching:
- Perfusion CT physics and mathematics
- Deconvolution principles and pitfalls
- Effect of noise on parameter estimation
- Hemodynamic state variations (ischemia, hyperperfusion)

## Credits

Jinhee Jang, MD, PhD  
Neuroradiology Section  
Department of Radiology  
Seoul St. Mary Hospital

## License

MIT License - See LICENSE file for details