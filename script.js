document.addEventListener('DOMContentLoaded', () => {
    // --- Constants & State ---
    const state = {
        cbf: 50,
        mtt: 4,
        tmax: 0,
        aifWidth: 1.0,
        noise: 0,
        model: 'boxcar', // 'boxcar' or 'exponential'
        timePoints: 60, // 60 seconds
        dt: 1, // 1 second resolution
        cardiacOutput: 5.0 // L/min - affects AIF/VOF amplitude
    };

    // --- DOM Elements ---
    const cbfSlider = document.getElementById('cbfSlider');
    const cbfValue = document.getElementById('cbfValue');
    const aifWidthSlider = document.getElementById('aifWidthSlider');
    const aifWidthValue = document.getElementById('aifWidthValue');
    const mttSlider = document.getElementById('mttSlider');
    const mttValue = document.getElementById('mttValue');
    const tmaxSlider = document.getElementById('tmaxSlider');
    const tmaxValue = document.getElementById('tmaxValue');
    const modelSelect = document.getElementById('modelSelect');
    const cbvValue = document.getElementById('cbvValue');
    const cardiacOutputSlider = document.getElementById('cardiacOutputSlider');
    const cardiacOutputValue = document.getElementById('cardiacOutputValue');
    const noiseSlider = document.getElementById('noiseSlider');
    const noiseValue = document.getElementById('noiseValue');

    // --- Chart Instances ---
    let aifChart, residueChart, tissueChart;

    // --- Math Functions ---

    // Simplified gamma function to avoid recursion issues
    function gamma(z) {
        if (z <= 0) return 1;
        if (z < 1) return gamma(z + 1) / z;
        if (z === 1) return 1;
        if (z === 2) return 1;
        if (z === 3) return 2;
        if (z === 4) return 6;
        if (z === 5) return 24;

        // For other values, use a simple approximation
        // Γ(z) ≈ √(2π/z) * (z/e)^z for large z
        if (z > 5) {
            return Math.sqrt(2 * Math.PI / z) * Math.pow(z / Math.E, z);
        }

        // For intermediate values, use linear interpolation
        return Math.pow(z - 1, z - 1) * Math.exp(-(z - 1));
    }

    // Gamma-variate function for AIF
    // C(t) = K * (t - t0)^alpha * exp(-(t - t0)/beta)
    // Low cardiac output causes delay, dispersion, and lower peak
    function generateAIF(timeSteps, widthScale, cardiacOutput = 5.0) {
        // Cardiac output effects:
        // High CO (>5) = Earlier arrival, sharper peak, higher amplitude
        // Low CO (<5) = Delayed arrival, wider/flatter curve, lower peak

        const normalCO = 5.0;
        const coRatio = cardiacOutput / normalCO;

        // 1. Arrival time: Lower CO = delayed arrival
        const baseArrival = 5;
        const t0 = baseArrival + (normalCO - cardiacOutput) * 0.5; // Delay with low CO

        // 2. Shape: Lower CO = more dispersion (higher beta)
        const alpha = 3.0;
        const baseBeta = 1.5 * widthScale;
        const beta = baseBeta * (2.0 - coRatio * 0.3); // Wider with low CO

        // 3. Amplitude: Higher CO = higher peak (better delivery)
        // Realistic AIF peak should be around 100-300 HU
        const baseAmplitude = 3.5; // Reduced for realistic HU values
        const K = baseAmplitude * Math.sqrt(coRatio) / Math.pow(widthScale, alpha);

        return timeSteps.map(t => {
            if (t < t0) return 0;
            const dt = t - t0;
            return K * Math.pow(dt, alpha) * Math.exp(-dt / beta);
        });
    }

    // Transport Function h(t) - CORRECTED for conservation of mass
    // Describes distribution of transit times without normalizing by 1/MTT
    // The area under h(t) represents the fraction that exits at each time point
    function generateTransport(timeSteps, mtt, delay, model) {
        const dt = 1; // 1 second intervals

        return timeSteps.map(t => {
            if (t < delay) return 0;

            if (model === 'boxcar') {
                // For boxcar: all blood exits at exactly t = delay + mtt
                // Use a narrow Gaussian to approximate delta function
                const transitTime = delay + mtt;
                const sigma = 0.5; // Width of approximation
                const diff = t - transitTime;
                // Gaussian approximation of delta function, normalized to preserve area
                return (1 / (sigma * Math.sqrt(2 * Math.PI))) *
                    Math.exp(-0.5 * Math.pow(diff / sigma, 2)) * dt;
            } else if (model === 'exponential') {
                // Simplified exponential model to avoid mathematical issues
                const x = t - delay;
                if (x <= 0) return 0;

                // Simple exponential decay with time constant = MTT
                const tau = Math.max(1, mtt); // Ensure positive time constant
                const expCoeff = Math.exp(-x / tau);

                // Protect against invalid values
                if (!isFinite(expCoeff) || isNaN(expCoeff)) {
                    return 0;
                }

                return expCoeff;
            }
            return 0;
        });
    }

    // Residue Function R(t)
    // Boxcar: R(t) = 1 for delay <= t < delay + MTT
    // Exponential: R(t) = exp(-(t - delay) / MTT) for t >= delay
    function generateResidue(timeSteps, mtt, delay, model) {
        return timeSteps.map(t => {
            if (t < delay) return 0;

            if (model === 'boxcar') {
                if (t >= delay && t < delay + mtt) return 1;
                return 0;
            } else if (model === 'exponential') {
                const expVal = Math.exp(-(t - delay) / mtt);
                // Protect against invalid values
                if (!isFinite(expVal) || isNaN(expVal)) return 0;
                return expVal;
            }
            return 0;
        });
    }

    // Convolution for tissue concentration
    // Ct(t) = CBF * (AIF * R)(t) * dt
    // Discrete convolution: (f * g)[n] = sum(f[m] * g[n-m])
    function convolve(aif, residue, cbf, dt) {
        const n = aif.length;
        const result = new Array(n).fill(0);

        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j <= i; j++) {
                sum += aif[j] * residue[i - j];
            }
            // Scale by CBF and dt
            // Note: CBF is usually in ml/100g/min. We need to be careful with units.
            // Here we treat CBF as a relative scalar for visualization.
            // If AIF is HU, Residue is unitless, result is HU.
            // CBF scaling: The standard equation is Ct(t) = CBF * conv(AIF, R).
            // Realistic tissue concentration should be 20-100 HU
            result[i] = sum * cbf * 0.001; // Reduced scaling for realistic tissue HU
        }
        return result;
    }

    // Convolution for VOF - preserves conservation of mass
    // VOF(t) = AIF(t) ⊗ h(t) where h(t) is transport function
    // Total area of VOF should equal total area of AIF (mass conservation)
    function convolveForVOF(aif, transport, dt) {
        const n = aif.length;
        const result = new Array(n).fill(0);

        // Standard convolution
        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j <= i; j++) {
                sum += aif[j] * transport[i - j];
            }
            result[i] = sum * dt;
        }

        // Normalize to preserve total area (conservation of mass)
        const aifArea = aif.reduce((sum, val) => sum + val, 0) * dt;
        const vofArea = result.reduce((sum, val) => sum + val, 0);

        if (vofArea > 0) {
            const normalizationFactor = aifArea / vofArea;
            return result.map(val => val * normalizationFactor);
        }

        return result;
    }

    // Impulse Response Function IRF(t) = CBF * R(t)
    function generateIRF(timeSteps, cbf, mtt, delay, model) {
        const residue = generateResidue(timeSteps, mtt, delay, model);
        return residue.map(r => r * cbf * 0.5); // Scale down for better visualization
    }

    // Simulated Recovered IRF with Noise Artifacts
    // Deconvolution is ill-posed. Noise in C_t(t) leads to oscillations in R(t).
    function generateRecoveredIRF(irfData, noiseLevel) {
        if (noiseLevel === 0) return irfData;

        // Add high-frequency oscillations
        return irfData.map((val, i) => {
            // Oscillations are often larger at the beginning or edges
            const oscillation = Math.sin(i * 1.5) * (noiseLevel * 2);
            // Also add random noise
            const random = (Math.random() - 0.5) * (noiseLevel * 2);
            return val + oscillation + random;
        });
    }

    function addNoise(data, noiseLevel) {
        if (noiseLevel === 0) return data;
        const maxVal = Math.max(...data);
        return data.map(val => {
            const noise = (Math.random() - 0.5) * 2 * (noiseLevel / 100) * maxVal;
            return Math.max(0, val + noise);
        });
    }

    // --- Annotation Plugin ---
    const annotationPlugin = {
        id: 'annotationPlugin',
        afterDraw: (chart) => {
            const ctx = chart.ctx;
            const xAxis = chart.scales.x;
            const yAxis = chart.scales.y;

            // Handle residue chart annotations (CBF, MTT, Tmax, CBV)
            if (chart.canvas.id === 'residueChart') {
                // Get data parameters from state (closure)
                const { cbf, mtt, tmax, model } = state;

                ctx.save();
                ctx.font = '12px Inter';
                ctx.fillStyle = '#fff';
                ctx.strokeStyle = '#fff';
                ctx.textAlign = 'center';

                // 1. CBF Arrow (Height)
                const scaledCBF = cbf * 0.5; // Match the IRF scaling
                const xCBF = xAxis.getPixelForValue(tmax + 0.5);
                const yCBFTop = yAxis.getPixelForValue(scaledCBF);
                const yCBFBottom = yAxis.getPixelForValue(0);

                // Draw Arrow
                ctx.beginPath();
                ctx.moveTo(xCBF, yCBFBottom);
                ctx.lineTo(xCBF, yCBFTop);
                ctx.stroke();

                // Arrowhead
                ctx.beginPath();
                ctx.moveTo(xCBF - 3, yCBFTop + 5);
                ctx.lineTo(xCBF, yCBFTop);
                ctx.lineTo(xCBF + 3, yCBFTop + 5);
                ctx.stroke();

                ctx.fillText('CBF', xCBF, yCBFTop - 5);

                // 2. Tmax (Delay)
                if (tmax > 0) {
                    const xStart = xAxis.getPixelForValue(0);
                    const xEnd = xAxis.getPixelForValue(tmax);
                    const yPos = yAxis.getPixelForValue(scaledCBF / 2); // Mid-height

                    ctx.beginPath();
                    ctx.moveTo(xStart, yPos);
                    ctx.lineTo(xEnd, yPos);
                    ctx.stroke();

                    // Brackets
                    ctx.beginPath();
                    ctx.moveTo(xStart, yPos - 3);
                    ctx.lineTo(xStart, yPos + 3);
                    ctx.moveTo(xEnd, yPos - 3);
                    ctx.lineTo(xEnd, yPos + 3);
                    ctx.stroke();

                    ctx.fillText('Tmax', (xStart + xEnd) / 2, yPos - 5);
                }

                // 3. MTT (Width or Decay)
                const yMTT = yAxis.getPixelForValue(scaledCBF * 0.8);
                const xMTTStart = xAxis.getPixelForValue(tmax);
                const xMTTEnd = xAxis.getPixelForValue(tmax + mtt);

                ctx.beginPath();
                ctx.moveTo(xMTTStart, yMTT);
                ctx.lineTo(xMTTEnd, yMTT);
                ctx.stroke();

                // Brackets
                ctx.beginPath();
                ctx.moveTo(xMTTStart, yMTT - 3);
                ctx.lineTo(xMTTStart, yMTT + 3);
                ctx.moveTo(xMTTEnd, yMTT - 3);
                ctx.lineTo(xMTTEnd, yMTT + 3);
                ctx.stroke();

                ctx.fillText('MTT', (xMTTStart + xMTTEnd) / 2, yMTT - 5);

                // 4. CBV (Area)
                const xArea = xAxis.getPixelForValue(tmax + mtt / 2);
                const yArea = yAxis.getPixelForValue(scaledCBF * 0.4);
                ctx.fillStyle = 'rgba(255, 255, 255, 0.7)';
                ctx.fillText('Area = CBV', xArea, yArea);

                ctx.restore();
            }

            // Handle AIF chart scan duration line (if needed)
            if (chart.canvas.id === 'aifChart' && state.scanDuration) {
                const xPos = xAxis.getPixelForValue(state.scanDuration);

                ctx.save();
                ctx.strokeStyle = '#ff6b6b';
                ctx.lineWidth = 2;
                ctx.setLineDash([10, 5]);

                ctx.beginPath();
                ctx.moveTo(xPos, yAxis.top);
                ctx.lineTo(xPos, yAxis.bottom);
                ctx.stroke();

                // Add label
                ctx.fillStyle = '#ff6b6b';
                ctx.font = '12px Inter';
                ctx.fillText('Scan End', xPos + 5, yAxis.top + 20);

                ctx.restore();
            }
        }
    };

    // --- Chart Configuration ---
    const commonOptions = {
        responsive: true,
        maintainAspectRatio: false,
        animation: { duration: 0 }, // Disable animation for performance during sliding
        interaction: {
            mode: 'index',
            intersect: false,
        },
        plugins: {
            legend: { display: false },
            tooltip: { enabled: true }
        },
        scales: {
            x: {
                grid: { color: '#333' },
                ticks: { color: '#888' },
                title: { display: true, text: 'Time (s)', color: '#666' }
            },
            y: {
                grid: { color: '#333' },
                ticks: { color: '#888' },
                beginAtZero: true
            }
        }
    };

    function initCharts() {
        const timeLabels = Array.from({ length: state.timePoints }, (_, i) => i);

        // AIF Chart
        const ctxAIF = document.getElementById('aifChart').getContext('2d');
        aifChart = new Chart(ctxAIF, {
            type: 'line',
            data: {
                labels: timeLabels,
                datasets: [{
                    label: 'Ideal AIF',
                    data: [],
                    borderColor: '#ff6b6b',
                    backgroundColor: 'rgba(255, 107, 107, 0.1)',
                    borderWidth: 2,
                    fill: true,
                    tension: 0.4,
                    order: 2
                }, {
                    label: 'Noisy AIF',
                    data: [],
                    borderColor: '#ff6b6b', // Same color but different style? Or maybe darker red?
                    backgroundColor: 'transparent',
                    borderWidth: 1,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.1,
                    order: 1 // Draw on top
                }, {
                    label: 'VOF (Output)',
                    data: [],
                    borderColor: '#4ecdc4',
                    borderDash: [5, 5],
                    borderWidth: 2,
                    fill: false,
                    tension: 0.4,
                    order: 3
                }]
            },
            options: {
                ...commonOptions,
                plugins: {
                    ...commonOptions.plugins,
                    title: { display: false },
                    legend: { display: true } // Show legend for AIF/VOF
                },
                scales: {
                    ...commonOptions.scales,
                    x: { ...commonOptions.scales.x, max: 50 }, // Extended to 50s
                    y: {
                        ...commonOptions.scales.y,
                        min: 0,
                        max: undefined, // Auto-scale for AIF
                        title: { display: true, text: 'Concentration (HU)', color: '#666' }
                    }
                }
            }
        });

        // Residue Chart (Now IRF)
        const ctxResidue = document.getElementById('residueChart').getContext('2d');
        residueChart = new Chart(ctxResidue, {
            type: 'line',
            data: {
                labels: timeLabels,
                datasets: [{
                    label: 'Ideal IRF (CBF*R)',
                    data: [],
                    borderColor: '#4ecdc4',
                    backgroundColor: 'rgba(78, 205, 196, 0.1)',
                    borderWidth: 2,
                    fill: true,
                    stepped: true
                }, {
                    label: 'Recovered IRF (Noisy)',
                    data: [],
                    borderColor: '#e57373',
                    borderWidth: 1,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.1
                }]
            },
            options: {
                ...commonOptions,
                plugins: {
                    ...commonOptions.plugins,
                    legend: { display: true }
                },
                scales: {
                    ...commonOptions.scales,
                    y: {
                        ...commonOptions.scales.y,
                        min: 0,
                        max: undefined, // Auto-scale for IRF
                        title: { display: true, text: 'CBF × R(t)', color: '#666' }
                    },
                    x: { ...commonOptions.scales.x, max: 50 } // Extended to 50s
                }
            },
            plugins: [annotationPlugin] // Register local plugin
        });

        // Tissue Chart
        const ctxTissue = document.getElementById('tissueChart').getContext('2d');
        tissueChart = new Chart(ctxTissue, {
            type: 'line',
            data: {
                labels: timeLabels,
                datasets: [{
                    label: 'Tissue Curve',
                    data: [],
                    borderColor: '#ffe66d',
                    backgroundColor: 'rgba(255, 230, 109, 0.1)',
                    borderWidth: 3,
                    fill: true,
                    tension: 0.4
                }]
            },
            options: {
                ...commonOptions,
                scales: {
                    ...commonOptions.scales,
                    y: {
                        ...commonOptions.scales.y,
                        min: 0,
                        max: undefined, // Auto-scale for tissue
                        title: { display: true, text: 'Tissue Concentration (HU)', color: '#666' }
                    },
                    x: { ...commonOptions.scales.x, max: 50 } // Extended to 50s
                }
            }
        });
    }

    // --- Update Logic ---
    function updateSimulation() {
        // Update State
        state.cbf = parseFloat(cbfSlider.value);
        state.aifWidth = parseFloat(aifWidthSlider.value);
        state.mtt = parseFloat(mttSlider.value);
        state.tmax = parseFloat(tmaxSlider.value);
        state.cardiacOutput = parseFloat(cardiacOutputSlider.value);
        state.noise = parseFloat(noiseSlider.value);
        state.model = modelSelect.value;

        // Update UI Text
        cbfValue.textContent = state.cbf;
        aifWidthValue.textContent = state.aifWidth;
        mttValue.textContent = state.mtt;
        tmaxValue.textContent = state.tmax;
        cardiacOutputValue.textContent = state.cardiacOutput;
        noiseValue.textContent = state.noise;

        // Calculate CBV = CBF * MTT / 60 (if CBF is per min and MTT in sec)
        // Usually CBV is ml/100g. CBF is ml/100g/min. MTT is sec.
        // CBV = (CBF / 60) * MTT
        const cbv = (state.cbf / 60) * state.mtt;
        cbvValue.textContent = cbv.toFixed(2);

        // Generate Data
        const timeSteps = Array.from({ length: state.timePoints }, (_, i) => i);
        const aifData = generateAIF(timeSteps, state.aifWidth, state.cardiacOutput);
        const residueData = generateResidue(timeSteps, state.mtt, state.tmax, state.model);
        const transportData = generateTransport(timeSteps, state.mtt, state.tmax, state.model);

        // IRF Data
        const irfData = generateIRF(timeSteps, state.cbf, state.mtt, state.tmax, state.model);
        const recoveredIrfData = generateRecoveredIRF(irfData, state.noise);

        // Calculate VOF = AIF ⊗ h(t) - Conservation of Mass
        // VOF should have same total area as AIF (mass conservation)
        const vofData = convolveForVOF(aifData, transportData, state.dt);

        // Convolve for Tissue
        let tissueData = convolve(aifData, residueData, state.cbf, state.dt);

        // Debug: Log first few values to check if data is generated
        console.log('AIF first 5:', aifData.slice(0, 5));
        console.log('Residue first 5:', residueData.slice(0, 5));
        console.log('IRF first 5:', irfData.slice(0, 5));
        console.log('Tissue first 5:', tissueData.slice(0, 5));

        // Add Noise
        let noisyAIF = aifData;
        if (state.noise > 0) {
            tissueData = addNoise(tissueData, state.noise);
            noisyAIF = addNoise(aifData, state.noise);
        }

        // Update Charts
        aifChart.data.datasets[0].data = aifData;
        aifChart.data.datasets[1].data = noisyAIF;
        aifChart.data.datasets[2].data = vofData;
        aifChart.update();

        residueChart.data.datasets[0].data = irfData;
        residueChart.data.datasets[1].data = recoveredIrfData;
        // Update stepped property based on model
        residueChart.data.datasets[0].stepped = state.model === 'boxcar';
        residueChart.update();

        tissueChart.data.datasets[0].data = tissueData;
        tissueChart.update();
    }

    // --- Debounce for better mobile performance ---
    function debounce(func, wait) {
        let timeout;
        return function executedFunction(...args) {
            const later = () => {
                clearTimeout(timeout);
                func(...args);
            };
            clearTimeout(timeout);
            timeout = setTimeout(later, wait);
        };
    }

    // Use immediate update for desktop, debounced for mobile
    const isMobile = window.matchMedia('(max-width: 640px)').matches;
    const updateHandler = isMobile ? debounce(updateSimulation, 150) : updateSimulation;

    // --- Event Listeners ---
    cbfSlider.addEventListener('input', updateHandler);
    aifWidthSlider.addEventListener('input', updateHandler);
    mttSlider.addEventListener('input', updateHandler);
    tmaxSlider.addEventListener('input', updateHandler);
    cardiacOutputSlider.addEventListener('input', updateHandler);
    modelSelect.addEventListener('change', updateSimulation);
    noiseSlider.addEventListener('input', updateHandler);

    // --- Initialization ---
    initCharts();
    updateSimulation();
});
