# Perfusion CT Simulation Tool

An interactive web-based educational simulation tool for understanding CT perfusion principles in stroke imaging.

## Overview

This tool demonstrates the physiological relationships between:
- **Arterial Input Function (AIF)** - contrast arrival in arteries
- **Venous Output Function (VOF)** - contrast washout through veins  
- **Impulse Response Function (IRF)** - CBF × Residue Function
- **Tissue Concentration Curve** - convolution result showing tissue perfusion

## Features

### Physiological Parameters
- **CBF (Cerebral Blood Flow)**: 10-100 ml/100g/min
- **MTT (Mean Transit Time)**: 1-10 seconds  
- **Tmax (Delay)**: 0-10 seconds
- **Cardiac Output**: 2-8 L/min (affects AIF/VOF amplitude and timing)

### Educational Concepts
- **Conservation of Mass**: VOF maintains same total area as AIF
- **Bolus Dispersion**: Lower cardiac output causes wider, delayed curves
- **Convolution Principle**: Tissue curve = CBF × [AIF ⊗ R(t)]
- **Flow Models**: Boxcar (plug flow) vs Exponential (perfect mixing)

### Interactive Features
- Real-time parameter adjustment with immediate visualization
- Mobile-responsive design with touch-optimized controls
- Graphical annotations showing CBF, MTT, Tmax, and CBV relationships
- Noise simulation for realistic clinical data representation

## Clinical Relevance

This simulation helps understand:
- How **low cardiac output** creates dangerous dispersion artifacts
- Why **VOF truncation** leads to CBF/CBV overestimation
- The relationship between **hemodynamic parameters** and curve shapes
- **Stroke imaging** principles in CT perfusion analysis

## Usage

Simply open `index.html` in a web browser - no server required!

Or visit the live demo: [https://znee.github.io/perfusionCT_simulation](https://znee.github.io/perfusionCT_simulation)

🔄 *Note: If the live demo link shows 404, GitHub Pages may still be deploying (allow 5-10 minutes after enabling Pages)*

## Technical Details

- **Pure client-side**: HTML5, CSS3, JavaScript (ES6)
- **Visualization**: Chart.js for interactive plots
- **Mathematics**: MathJax for formula rendering
- **Responsive**: Works on desktop, tablet, and mobile devices

## Credits

**Jinhee Jang, MD. PhD**  
Neuroradiology Section, Department of Radiology  
Seoul St. Mary Hospital

## License

MIT License - free for educational and research use.