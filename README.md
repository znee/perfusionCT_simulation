# Perfusion CT Simulation

Interactive educational tool for visualizing CT perfusion physics and hemodynamic parameters.

## Features

- **Real-time Visualization**: Interactive sliders for CBF, MTT, Tmax, cardiac output, and noise
- **Physiologically Accurate**: Implements proper convolution mathematics with realistic HU values
- **Multiple Models**: Boxcar (plug flow) and Exponential (perfect mixing) residue functions
- **Educational Annotations**: Visual markers for CBF, MTT, Tmax, and CBV on IRF chart
- **Noise Simulation**: Demonstrates deconvolution artifacts from noisy tissue data

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
Visit the live demo at your GitHub Pages URL.

## Physiological Parameters

### CBF (Cerebral Blood Flow)
- Range: 10-100 ml/100g/min
- Normal gray matter: 60-80 ml/100g/min
- Normal white matter: 20-30 ml/100g/min

### MTT (Mean Transit Time)
- Range: 1-10 seconds
- Represents average time for blood to traverse tissue capillary bed

### Tmax (Delay)
- Range: 0-10 seconds
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
- `R(t)`: Residue function (unitless, 0-1)
- `⊗`: Convolution operator

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
- Delay (Tmax) applied by shifting AIF in convolution
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

## Chart Details

### AIF Chart (0-70s)
- **Ideal AIF** (solid red): Gamma-variate, peak ~200 HU
- **Noisy AIF** (thin red): With added noise
- **VOF** (dashed teal): Venous output (mass-conserved)

### IRF Chart (0-70s)
- **Ideal IRF** (dashed light teal): Canonical CBF×R(t) at t=0
- **Shifted IRF** (solid dark teal): Response delayed by Tmax
- **Recovered IRF** (thin red): Noisy deconvolution result

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