document.addEventListener('DOMContentLoaded', () => {
    // Check if numeric.js is loaded
    if (typeof numeric === 'undefined') {
        console.error('CRITICAL ERROR: numeric.js library not loaded! Deconvolution will fail.');
        alert('Error: numeric.js library not loaded. Please refresh the page.');
    }

    // --- Constants & State ---
    const state = {
        cbf: 50,
        mtt: 6,
        tmax: 0,
        aifWidth: 1.0,
        noise: 0,
        deconvMode: 'insensitive', // 'sensitive' or 'insensitive'
        model: 'exponential', // 'boxcar' or 'exponential' - exponential is default
        timePoints: 60, // 60 seconds
        dt: 1, // 1 second resolution
        cardiacOutput: 5.0, // L/min - affects AIF/VOF amplitude
        lambda: 0.05 // Deconvolution regularization parameter
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
    const lambdaSlider = document.getElementById('lambdaSlider');
    const lambdaValue = document.getElementById('lambdaValue');
    const deconvModeSelect = document.getElementById('deconvModeSelect');

    // Comparison table elements
    const cbfOriginal = document.getElementById('cbfOriginal');
    const cbvOriginal = document.getElementById('cbvOriginal');
    const tmaxOriginal = document.getElementById('tmaxOriginal');
    const mttOriginal = document.getElementById('mttOriginal');
    const cbfDeconv = document.getElementById('cbfDeconv');
    const cbvDeconv = document.getElementById('cbvDeconv');
    const tmaxDeconv = document.getElementById('tmaxDeconv');
    const mttDeconv = document.getElementById('mttDeconv');
    const cbfError = document.getElementById('cbfError');
    const cbvError = document.getElementById('cbvError');
    const tmaxError = document.getElementById('tmaxError');
    const mttError = document.getElementById('mttError');

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

    // Gamma-variate function for AIF - CALIBRATED to 150-250 HU
    // C(t) = K * (t - t0)^alpha * exp(-(t - t0)/beta)
    function generateAIF(timeSteps, widthScale, cardiacOutput = 5.0) {
        const normalCO = 5.0;
        const coRatio = cardiacOutput / normalCO;

        // 1. Arrival time: Lower CO = delayed arrival
        const baseArrival = 5;
        const t0 = baseArrival + (normalCO - cardiacOutput) * 0.5;

        // 2. Shape parameters
        const alpha = 3.0;
        const baseBeta = 1.5 * widthScale;
        const beta = baseBeta * (2.0 - coRatio * 0.3);

        // 3. Fix 4: Systematic calibration to target peak HU
        // Generate raw AIF, find max, scale to target
        const rawAIF = timeSteps.map(t => {
            if (t < t0) return 0;
            const dt = t - t0;
            return Math.pow(dt, alpha) * Math.exp(-dt / beta);
        });

        const rawMax = Math.max(...rawAIF);
        const targetPeakHU = 200; // Clinical range: 150-300 HU
        const scaleFactor = targetPeakHU / rawMax;

        return rawAIF.map(val => val * scaleFactor * Math.sqrt(coRatio));
    }

    // Transport Function h(t) - Mass-conservative
    // Describes distribution of transit times
    // Normalized to integrate to 1 (dt applied in convolution)
    function generateTransport(timeSteps, mtt, delay, model) {
        return timeSteps.map(t => {
            if (t < delay) return 0;

            if (model === 'boxcar') {
                // Delta function approximation at t = delay + mtt
                // Gaussian normalized to area = 1 (no dt here)
                const transitTime = delay + mtt;
                const sigma = Math.max(0.7 * 1, 0.5); // sigma >= 0.7*dt for stable discrete area
                const diff = t - transitTime;
                return (1 / (sigma * Math.sqrt(2 * Math.PI))) *
                    Math.exp(-0.5 * Math.pow(diff / sigma, 2));
            } else if (model === 'exponential') {
                // Exponential with proper normalization: (1/τ) * exp(-x/τ)
                const x = t - delay;
                if (x <= 0) return 0;

                const tau = Math.max(1, mtt);
                return (1 / tau) * Math.exp(-x / tau); // Added 1/tau normalization
            }
            return 0;
        });
    }

    // Residue Function R(t) - CANONICAL (no delay)
    // R(0) = 1 for both models (correct for deconvolution teaching)
    // Delay is handled by shifting AIF in convolution
    function generateResidue(timeSteps, mtt, model) {
        return timeSteps.map(t => {
            if (model === 'boxcar') {
                return (t >= 0 && t < mtt) ? 1 : 0;
            } else if (model === 'exponential') {
                return (t >= 0) ? Math.exp(-t / mtt) : 0;
            }
            return 0;
        });
    }

    // Convolution for tissue concentration with AIF delay shift
    // Ct(t) = CBF * (AIF * R)(t) * dt
    // Tmax applied as AIF shift (arterial-tissue delay) with fractional interpolation
    // Units: CBF (ml/100g/min) needs conversion to (ml/100g/sec)
    function convolve(aif, residue, cbf, dt, tmax) {
        const n = aif.length;
        const result = new Array(n).fill(0);

        // Fractional delay support via linear interpolation
        const shifted = tmax / dt;
        const i0 = Math.floor(shifted);
        const frac = shifted - i0;

        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j <= i; j++) {
                // Interpolated AIF for fractional delay
                const k = Math.max(0, j - i0);
                const k1 = Math.max(0, k - 1);
                const interpolatedAIF = (1 - frac) * aif[k] + frac * aif[k1];
                sum += interpolatedAIF * residue[i - j];
            }
            result[i] = (cbf / 60) * dt * sum * 0.1;
        }
        return result;
    }

    // Convolution for VOF - Auto mass-conservative with normalized h(t)
    // VOF(t) = AIF(t) ⊗ h(t) where h(t) is transport function
    function convolveForVOF(aif, transport, dt) {
        const n = aif.length;
        const result = new Array(n).fill(0);

        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j <= i; j++) {
                sum += aif[j] * transport[i - j];
            }
            result[i] = sum * dt;
        }

        // Fix 6: Optional normalization (only if >1% error)
        const aifArea = aif.reduce((s, v) => s + v, 0) * dt;
        const vofArea = result.reduce((s, v) => s + v, 0);

        if (Math.abs(vofArea - aifArea) / aifArea > 0.01) {
            // Apply correction only if > 1% error
            const factor = aifArea / vofArea;
            return result.map(val => val * factor);
        }

        return result;
    }

    // Impulse Response Function IRF(t) = CBF * R(t)
    // Fix 5: No arbitrary scaling - literal CBF × R(t)
    function generateIRF(timeSteps, cbf, mtt, model) {
        const residue = generateResidue(timeSteps, mtt, model);
        return residue.map(r => r * cbf);
    }

    // Shifted IRF for visualization (shows effective IRF in tissue with Tmax delay)
    // Shifted Residue Function for Forward Model
    // Applies Tmax delay directly to the residue function R(t)
    function generateShiftedResidue(residue, tmax, dt) {
        const n = residue.length;
        const result = new Array(n).fill(0);

        // Fractional delay support via linear interpolation
        const shifted = tmax / dt;
        const i0 = Math.floor(shifted);
        const frac = shifted - i0;

        for (let i = 0; i < n; i++) {
            // R_shifted[i] = R[i - delay]
            // Interpolate between i - i0 and i - i0 - 1
            const k = i - i0;
            const k1 = k - 1;

            let val = 0;
            if (k >= 0 && k < n) val += (1 - frac) * residue[k];
            if (k1 >= 0 && k1 < n) val += frac * residue[k1];

            result[i] = val;
        }
        return result;
    }

    // Shifted IRF for visualization (shows effective IRF in tissue with Tmax delay)
    function generateShiftedIRF(irfData, tmax, dt) {
        // Reuse the shifted residue logic but for IRF
        return generateShiftedResidue(irfData, tmax, dt);
    }

    function addNoise(data, noiseLevel) {
        if (noiseLevel === 0) return data;
        const maxVal = Math.max(...data);
        return data.map(val => {
            const noise = (Math.random() - 0.5) * 2 * (noiseLevel / 100) * maxVal;
            return Math.max(0, val + noise);
        });
    }

    // --- Deconvolution Functions ---

    // Build Toeplitz matrix from AIF for convolution operation
    // Each row represents the AIF shifted by one time step
    function buildToeplitzMatrix(aif) {
        const n = aif.length;
        const matrix = [];
        for (let i = 0; i < n; i++) {
            const row = new Array(n).fill(0);
            for (let j = 0; j <= i; j++) {
                row[j] = aif[i - j];
            }
            matrix.push(row);
        }
        return matrix;
    }

    // Build Block-Circulant Matrix for Delay-Insensitive Deconvolution
    // Lower triangle: standard Toeplitz (causal)
    // Upper triangle: wraps from lower-left (allows IRF shift beyond window)
    function buildCirculantMatrix(aif) {
        const N = aif.length;
        const matrix = [];

        for (let i = 0; i < N; i++) {
            const row = new Array(N).fill(0);
            for (let j = 0; j < N; j++) {
                if (j <= i) {
                    // Lower triangle: standard Toeplitz convolution
                    row[j] = aif[i - j];
                } else {
                    // Upper triangle: wrap (allows delayed IRF to be recovered)
                    row[j] = aif[i - j + N];
                }
            }
            matrix.push(row);
        }

        return matrix;
    }

    // SVD-based deconvolution with Tikhonov regularization
    // Solves: C_t(t) = A  * R(t) where A is Toeplitz or Circulant matrix from AIF
    // Returns: Deconvolved IRF (scaled residue function)
    function deconvolveSVD(tissueCurve, aif, dt, lambda, mode) {
        try {
            // Choose matrix construction based on mode
            let A;
            if (mode === 'insensitive') {
                A = buildCirculantMatrix(aif);
            } else {
                // Default to standard Toeplitz (sensitive to delay)
                A = buildToeplitzMatrix(aif);
            }

            // Scale tissue curve and matrix by CBF conversion factor and dt
            // This reverses the scaling in convolve(): (cbf/60) * dt * 0.1
            // We want to recover CBF * R(t) where CBF is in ml/100g/min
            // Tissue = (CBF/60) * dt * 0.1 * (A * R)
            // Tissue = (CBF * R) * (dt * 0.1 / 60) * A
            // Let x = CBF * R (what we want)
            // Tissue = (Tissue * 60) / (dt * 0.1 * A)
            // So we scale tissue by 60 / (dt * 0.1)
            const scaleFactor = (dt * 0.1) / 60;
            const b = tissueCurve.map(val => val / scaleFactor);

            // Perform SVD using numeric.js: A = U * S * V^T
            const svd = numeric.svd(A);
            const U = svd.U;
            const S = svd.S; // Singular values (1D array)
            const V = svd.V;

            // Tikhonov regularization: invert singular values with pure Tikhonov damping
            // S_inv[i] = S[i] / (S[i]^2 + (lambda * maxS)^2)
            const maxS = Math.max(...S);
            const threshold = lambda * maxS; // Relative threshold

            const S_inv = S.map(s => {
                // Pure Tikhonov regularization (no hard truncation)
                return s / (s * s + threshold * threshold);
            });

            // Reconstruct solution: R = V * S_inv * U^T * b
            // Step 1: U^T * b
            const UtB = numeric.dot(numeric.transpose(U), b);

            // Step 2: S_inv * (U^T * b) - element-wise multiplication
            const S_invUtB = UtB.map((val, i) => val * S_inv[i]);

            // Step 3: V * (S_inv * U^T * b)
            const deconvolved = numeric.dot(V, S_invUtB);

            // Ensure non-negative and handle numerical artifacts
            return deconvolved.map(val => Math.max(0, val));

        } catch (error) {
            console.error('Deconvolution error:', error);
            // Return zeros if deconvolution fails
            return new Array(tissueCurve.length).fill(0);
        }
    }

    // Extract perfusion parameters from deconvolved IRF
    // Returns: { cbf, cbv, tmax, mtt }
    function deriveParametersFromDeconvolution(deconvolvedIRF, dt) {
        // CBF: Peak value of deconvolved IRF (in ml/100g/min)
        const cbf_deconv = Math.max(...deconvolvedIRF);

        // Find time to peak (Tmax) in seconds
        const peakIndex = deconvolvedIRF.indexOf(cbf_deconv);
        const tmax_deconv = peakIndex * dt;

        // CBV: Area under IRF curve (integral of IRF with unit conversion)
        // IRF units: ml/100g/min, so divide by 60 to get ml/100g/sec, then multiply by dt
        const cbv_deconv = deconvolvedIRF.reduce((sum, val) => sum + val, 0) * dt / 60;

        // MTT: Central Volume Theorem - MTT = CBV / CBF
        // CBV is in ml/100g, CBF is in ml/100g/min
        // MTT = (CBV * 60) / CBF to get seconds
        const mtt_deconv = cbf_deconv > 0 ? (cbv_deconv * 60) / cbf_deconv : 0;

        return {
            cbf: cbf_deconv,
            cbv: cbv_deconv,
            tmax: tmax_deconv,
            mtt: mtt_deconv
        };
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
                ctx.font = '16px Inter';
                ctx.fillStyle = '#fff';
                ctx.strokeStyle = '#fff';
                ctx.textAlign = 'center';

                // 1. CBF Arrow (Height)
                // IRF = CBF × R(t), so height is literally CBF
                const xCBF = xAxis.getPixelForValue(tmax + 0.5);
                const yCBFTop = yAxis.getPixelForValue(cbf); // Literal CBF height
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
                // Note: IRF is canonical (no delay), but Tmax represents
                // the arterial-tissue delay applied in convolution
                if (tmax > 0) {
                    const xStart = xAxis.getPixelForValue(0);
                    const xEnd = xAxis.getPixelForValue(tmax);
                    const yPos = yAxis.getPixelForValue(cbf / 2); // Mid-height

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
                const yMTT = yAxis.getPixelForValue(cbf * 0.8);
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
                const yArea = yAxis.getPixelForValue(cbf * 0.4);
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
                beginAtZero: true,
                grace: '10%'
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
                    x: { ...commonOptions.scales.x, max: 70 }, // Extended to 70s
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
                    borderColor: '#95e1d3',
                    borderDash: [5, 5],
                    backgroundColor: 'transparent',
                    borderWidth: 2,
                    fill: false,
                    stepped: true,
                    order: 4
                }, {
                    label: 'Shifted IRF (with Tmax)',
                    data: [],
                    borderColor: '#4ecdc4',
                    backgroundColor: 'rgba(78, 205, 196, 0.1)',
                    borderWidth: 2,
                    fill: true,
                    stepped: true,
                    order: 3
                }, {
                    label: 'Deconvolved IRF',
                    data: [],
                    borderColor: '#f59e0b',
                    backgroundColor: 'transparent',
                    borderWidth: 3,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.2,
                    order: 1
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
                    x: { ...commonOptions.scales.x, max: 70 } // Extended to 70s
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
                    x: { ...commonOptions.scales.x, max: 70 } // Extended to 70s
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
        state.lambda = parseFloat(lambdaSlider.value);
        state.model = modelSelect.value;

        // Update UI Text
        cbfValue.textContent = state.cbf;
        aifWidthValue.textContent = state.aifWidth;
        mttValue.textContent = state.mtt;
        tmaxValue.textContent = state.tmax;
        cardiacOutputValue.textContent = state.cardiacOutput;
        noiseValue.textContent = state.noise;
        lambdaValue.textContent = state.lambda.toFixed(2);

        // Calculate CBV = CBF * MTT / 60 (if CBF is per min and MTT in sec)
        // Usually CBV is ml/100g. CBF is ml/100g/min. MTT is sec.
        // CBV = (CBF / 60) * MTT
        const cbv = (state.cbf / 60) * state.mtt;
        cbvValue.textContent = cbv.toFixed(2);

        // Generate Data
        const timeSteps = Array.from({ length: state.timePoints }, (_, i) => i);
        const aifData = generateAIF(timeSteps, state.aifWidth, state.cardiacOutput);
        const residueData = generateResidue(timeSteps, state.mtt, state.model); // No delay
        const transportData = generateTransport(timeSteps, state.mtt, state.tmax, state.model);

        // IRF Data (no delay in function, tmax handled in convolution)
        const irfData = generateIRF(timeSteps, state.cbf, state.mtt, state.model);
        const shiftedIrfData = generateShiftedIRF(irfData, state.tmax, state.dt);

        // Calculate VOF = AIF ⊗ h(t) - Auto mass-conservative
        const vofData = convolveForVOF(aifData, transportData, state.dt);

        // Convolve for Tissue (with tmax as AIF shift)
        const tissueDataRaw = convolve(aifData, residueData, state.cbf, state.dt, state.tmax);
        const tissueData = addNoise(tissueDataRaw, state.noise);

        // Perform Deconvolution on (potentially noisy) tissue data
        const aifForDeconv = addNoise(aifData, state.noise);
        const deconvolvedIRF = deconvolveSVD(tissueData, aifForDeconv, state.dt, state.lambda, state.deconvMode);

        // Derive parameters from deconvolution
        const deconvParams = deriveParametersFromDeconvolution(deconvolvedIRF, state.dt);

        // Update ground truth (original) parameters in table
        cbfOriginal.textContent = state.cbf.toFixed(1);
        cbvOriginal.textContent = cbv.toFixed(2);
        tmaxOriginal.textContent = state.tmax.toFixed(1);
        mttOriginal.textContent = state.mtt.toFixed(1);

        // Update deconvolved parameters in table
        cbfDeconv.textContent = deconvParams.cbf.toFixed(1);
        cbvDeconv.textContent = deconvParams.cbv.toFixed(2);
        tmaxDeconv.textContent = deconvParams.tmax.toFixed(1);
        mttDeconv.textContent = deconvParams.mtt.toFixed(1);

        // Calculate and display errors
        const cbfErrorPct = ((Math.abs(deconvParams.cbf - state.cbf) / state.cbf) * 100).toFixed(1);
        const cbvErrorPct = ((Math.abs(deconvParams.cbv - cbv) / cbv) * 100).toFixed(1);
        const tmaxErrorAbs = Math.abs(deconvParams.tmax - state.tmax).toFixed(1);
        const mttErrorPct = ((Math.abs(deconvParams.mtt - state.mtt) / state.mtt) * 100).toFixed(1);

        cbfError.textContent = `${cbfErrorPct}%`;
        cbvError.textContent = `${cbvErrorPct}%`;
        tmaxError.textContent = `${tmaxErrorAbs}s`;
        mttError.textContent = `${mttErrorPct}%`;

        // Color-code errors (green if < 10%, yellow if < 25%, red if >= 25%)
        // For badges, we set the background color and ensure text is white/readable
        const getBadgeColor = (val, threshold1, threshold2) => {
            return val < threshold1 ? 'rgba(129, 199, 132, 0.2)' : val < threshold2 ? 'rgba(255, 183, 77, 0.2)' : 'rgba(229, 115, 115, 0.2)';
        };

        const getTextColor = (val, threshold1, threshold2) => {
            return val < threshold1 ? '#81c784' : val < threshold2 ? '#ffb74d' : '#e57373';
        };

        cbfError.style.backgroundColor = getBadgeColor(parseFloat(cbfErrorPct), 10, 25);
        cbfError.style.color = getTextColor(parseFloat(cbfErrorPct), 10, 25);

        cbvError.style.backgroundColor = getBadgeColor(parseFloat(cbvErrorPct), 10, 25);
        cbvError.style.color = getTextColor(parseFloat(cbvErrorPct), 10, 25);

        tmaxError.style.backgroundColor = getBadgeColor(parseFloat(tmaxErrorAbs), 0.5, 1.0);
        tmaxError.style.color = getTextColor(parseFloat(tmaxErrorAbs), 0.5, 1.0);

        mttError.style.backgroundColor = getBadgeColor(parseFloat(mttErrorPct), 10, 25);
        mttError.style.color = getTextColor(parseFloat(mttErrorPct), 10, 25);

        // Fix 7: CBV numeric verification (with unit conversion)
        // IRF has units of CBF (ml/100g/min), so divide by 60 for ml/100g/sec
        const computedCBV = irfData.reduce((sum, val) => sum + val, 0) * state.dt / 60;
        const theoreticalCBV = (state.cbf / 60) * state.mtt;
        // console.log(`CBV numeric: ${computedCBV.toFixed(2)}, theoretical: ${theoreticalCBV.toFixed(2)}, error: ${((Math.abs(computedCBV - theoreticalCBV) / theoreticalCBV) * 100).toFixed(1)}%`);

        // Debug: Log deconvolution results
        // console.log('Deconvolution Results:', {
        //     original: { cbf: state.cbf, cbv: cbv.toFixed(2), tmax: state.tmax, mtt: state.mtt },
        //     deconvolved: deconvParams,
        //     errors: {
        //         cbf: ((Math.abs(deconvParams.cbf - state.cbf) / state.cbf) * 100).toFixed(1) + '%',
        //         cbv: ((Math.abs(deconvParams.cbv - cbv) / cbv) * 100).toFixed(1) + '%',
        //         tmax: Math.abs(deconvParams.tmax - state.tmax).toFixed(1) + 's',
        //         mtt: ((Math.abs(deconvParams.mtt - state.mtt) / state.mtt) * 100).toFixed(1) + '%'
        //     }
        // });

        // Update Charts
        aifChart.data.datasets[0].data = aifData;
        aifChart.data.datasets[1].data = aifForDeconv;
        aifChart.data.datasets[2].data = vofData;
        aifChart.update();

        residueChart.data.datasets[0].data = irfData;
        residueChart.data.datasets[1].data = shiftedIrfData;
        residueChart.data.datasets[2].data = deconvolvedIRF; // Deconvolved IRF
        // Update stepped property based on model
        residueChart.data.datasets[0].stepped = state.model === 'boxcar';
        residueChart.data.datasets[1].stepped = state.model === 'boxcar';
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
    deconvModeSelect.addEventListener('change', (e) => {
        state.deconvMode = e.target.value;
        updateSimulation();
    });
    noiseSlider.addEventListener('input', updateHandler);
    lambdaSlider.addEventListener('input', updateHandler);

    // --- Initialization ---
    initCharts();
    // Initial simulation update
    updateSimulation();

    // Trigger MathJax typeset if available
    if (typeof MathJax !== 'undefined' && MathJax.typeset) {
        MathJax.typeset();
    }
});
