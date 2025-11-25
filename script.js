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
        timePoints: 121, // 0-120 seconds (inclusive)
        dt: 1, // 1 second resolution
        cardiacOutput: 5.0, // L/min - affects AIF/VOF amplitude
        lambda: 0.05, // Deconvolution regularization parameter
        scanDuration: 60, // seconds slider
        scanDurationEffective: 60, // actual truncated duration applied to data
        noiseModel: 'poisson'
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
    const noiseModelSelect = document.getElementById('noiseModelSelect');
    const lambdaSlider = document.getElementById('lambdaSlider');
    const lambdaValue = document.getElementById('lambdaValue');
    const deconvModeSelect = document.getElementById('deconvModeSelect');
    const scanDurationSlider = document.getElementById('scanDurationSlider');
    const scanDurationValue = document.getElementById('scanDurationValue');
    const scanDurationEffectiveDisplay = document.getElementById('scanDurationEffectiveDisplay');

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
    let flowChart, residueChart;

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
    function computeVofDelay(cardiacOutput) {
        const normalCO = 5.0;
        const baseDelay = 4.0; // seconds for normal venous transit
        const coFactor = normalCO / Math.max(1.0, cardiacOutput);
        const delay = baseDelay * coFactor;
        return Math.min(12, Math.max(2, delay));
    }

    function generateTransport(timeSteps, mtt, delay, model, dt) {
        return timeSteps.map(t => {
            if (t < delay) return 0;

            if (model === 'boxcar') {
                // Delta function approximation at t = delay + mtt
                // Gaussian normalized to area = 1 (no dt here)
                const transitTime = delay + mtt;
                const sigma = Math.max(0.7 * dt, 0.5); // sigma >= 0.7*dt for stable discrete area
                const diff = t - transitTime;
                return (1 / (sigma * Math.sqrt(2 * Math.PI))) *
                    Math.exp(-0.5 * Math.pow(diff / sigma, 2));
            } else if (model === 'exponential') {
                // Gamma-variate transport: peaks at delay + mtt
                // Shape parameter = 2 gives realistic peaked distribution
                const x = t - delay;
                if (x <= 0) return 0;

                const tau = Math.max(1, mtt);
                const shape = 2.0; // Shape parameter for gamma distribution
                // Gamma PDF: (x^(k-1) * exp(-x/θ)) / (θ^k * Γ(k))
                // For k=2: (x * exp(-x/θ)) / θ^2
                // Peak at x = (k-1)*θ = 1*θ = mtt
                return (x / (tau * tau)) * Math.exp(-x / tau);
            }
            return 0;
        });
    }

    // Residue Function R(t) - CANONICAL (no delay)
    // R(0) = 1 for both models (correct for deconvolution teaching)
    // Delay is handled by shifting the residue (forward model matches inverse)
    function generateResidue(timeSteps, mtt, model, dt) {
        return timeSteps.map(t => {
            if (model === 'boxcar') {
                const remaining = mtt - t;
                if (remaining <= 0) return 0;
                if (remaining >= dt) return 1;
                return remaining / dt;
            } else if (model === 'exponential') {
                return (t >= 0) ? Math.exp(-t / mtt) : 0;
            }
            return 0;
        });
    }

    // Convolution for tissue concentration with residue delay shift
    // Ct(t) = CBF * (AIF * R)(t) * dt
    // Tmax applied to the residue (arterial-tissue delay) with fractional interpolation
    // Units: CBF (ml/100g/min) needs conversion to (ml/100g/sec)
    function convolve(aif, residue, cbf, dt, tmax, perfusionScale = 1, circular = false) {
        const n = aif.length;
        const result = new Array(n).fill(0);
        let delayedResidue = residue;
        if (tmax > 0) {
            if (circular) {
                delayedResidue = residue.map((_, i) => {
                    const shift = Math.floor(tmax / dt);
                    const frac = (tmax / dt) - shift;
                    const idx = (i - shift + n) % n;
                    const idxPrev = (idx - 1 + n) % n;
                    return (1 - frac) * residue[idx] + frac * residue[idxPrev];
                });
            } else {
                delayedResidue = generateShiftedResidue(residue, tmax, dt);
            }
        }

        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j <= i; j++) {
                sum += aif[j] * delayedResidue[i - j];
            }
            result[i] = (cbf / 60) * dt * sum * 0.1 * perfusionScale;
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
    function generateIRF(timeSteps, cbf, mtt, model, dt) {
        const residue = generateResidue(timeSteps, mtt, model, dt);
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

    function gaussianRandom() {
        let u = 0, v = 0;
        while (u === 0) u = Math.random();
        while (v === 0) v = Math.random();
        return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    }

    function poissonRandom(mean) {
        if (!Number.isFinite(mean) || mean <= 0) return 0;
        if (mean > 50) {
            // Normal approximation for efficiency at high means
            return Math.max(0, Math.round(mean + Math.sqrt(mean) * gaussianRandom()));
        }
        const L = Math.exp(-mean);
        let p = 1;
        let k = 0;
        do {
            k++;
            p *= Math.random();
        } while (p > L);
        return k - 1;
    }

    function addNoise(data, noiseLevel, model = 'gaussian') {
        if (noiseLevel === 0) return [...data];
        const maxVal = Math.max(...data);
        const scale = noiseLevel / 100;

        return data.map(val => {
            if (model === 'poisson') {
                // CT Physics: Noise proportional to 1/√(photon_count)
                // Simulate realistic photon statistics

                const basePhotonCount = 100000; // Typical CT: 10^5 photons/pixel
                const baseHU = 50; // Reference HU for normalization

                // Higher HU (more attenuation) = fewer transmitted photons
                // Approximate relationship: photons ∝ exp(-μ*x) ∝ 1/(1 + HU/1000)
                const attenuationFactor = 1 + Math.abs(val) / 1000;
                const effectivePhotons = basePhotonCount / attenuationFactor;

                // Poisson noise in photon domain: σ = √N
                const noisyPhotons = poissonRandom(effectivePhotons);
                const relativeNoise = (noisyPhotons - effectivePhotons) / Math.sqrt(effectivePhotons);

                // Convert to HU domain
                // Noise in HU ≈ (ΔN/N) * sensitivity, where sensitivity relates to CT numbers
                const huNoise = relativeNoise * (baseHU / Math.sqrt(basePhotonCount / 100)) * scale * maxVal;

                return Math.max(0, val + huNoise);
            }

            // Gaussian (electronics noise - constant variance)
            const noise = gaussianRandom() * maxVal * scale;
            return Math.max(0, val + noise);
        });
    }

    function getAcquisitionSamples(scanDuration, dt, maxPoints) {
        const clamped = Math.max(0, scanDuration);
        const samples = Math.floor(clamped / dt) + 1;
        return Math.min(maxPoints, Math.max(1, samples));
    }

    function applyAcquisitionWindow(data, samples, totalLength) {
        const output = new Array(totalLength).fill(null);
        const limit = Math.min(totalLength, data.length);
        for (let i = 0; i < limit; i++) {
            if (i < samples) {
                output[i] = data[i];
            }
        }
        return output;
    }

    function getCardiacPerfusionScale(cardiacOutput) {
        const normalCO = 5.0;
        const ratio = Math.max(0.2, cardiacOutput / normalCO);
        return Math.min(1.4, Math.max(0.6, Math.pow(ratio, 0.9)));
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
    function deconvolveSVD(tissueCurve, aif, dt, lambda, mode, perfusionScale = 1) {
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
            const scaleFactor = (dt * 0.1 * perfusionScale) / 60;
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

            return deconvolved;

        } catch (error) {
            console.error('Deconvolution error:', error);
            // Return zeros if deconvolution fails
            return new Array(tissueCurve.length).fill(0);
        }
    }

    // Extract perfusion parameters from deconvolved IRF
    // Returns: { cbf, cbv, tmax, mtt }
    function deriveParametersFromDeconvolution(deconvolvedIRF, dt, maxTime = Infinity) {
        if (!deconvolvedIRF.length) {
            return { cbf: 0, cbv: 0, tmax: 0, mtt: 0 };
        }

        const cbfPeak = Math.max(...deconvolvedIRF);
        if (!Number.isFinite(cbfPeak)) {
            return { cbf: 0, cbv: 0, tmax: 0, mtt: 0 };
        }

        const cbf = Math.max(0, cbfPeak);
        const peakIndex = Math.max(0, deconvolvedIRF.indexOf(cbfPeak));
        const rawTmax = peakIndex * dt;
        const boundedWindow = Number.isFinite(maxTime) ? maxTime : Infinity;
        const tmax = Math.max(0, Math.min(boundedWindow, rawTmax));

        const cbvRaw = deconvolvedIRF.reduce((sum, val) => sum + (Number.isFinite(val) ? val : 0), 0) * dt / 60;
        const cbv = Math.max(0, cbvRaw);

        const mttRaw = cbf > 0 ? (cbv * 60) / cbf : 0;
        const mtt = Math.max(0, Math.min(boundedWindow, mttRaw));

        return {
            cbf,
            cbv,
            tmax,
            mtt
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
            if ((chart.canvas.id === 'flowChart' || chart.canvas.id === 'residueChart') && state.scanDurationEffective) {
                const xPos = xAxis.getPixelForValue(state.scanDurationEffective);

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
                const label = chart.canvas.id === 'residueChart' ? 'Scan Window' : 'Scan End';
                ctx.fillText(label, xPos + 5, yAxis.top + 20);

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
        const timeLabels = Array.from({ length: state.timePoints }, (_, i) => i * state.dt);

        // Combined Flow Chart (AIF, Tissue, VOF)
        const ctxFlow = document.getElementById('flowChart').getContext('2d');
        flowChart = new Chart(ctxFlow, {
            type: 'line',
            data: {
                labels: timeLabels,
                datasets: [{
                    label: 'Ideal AIF (Input)',
                    data: [],
                    borderColor: '#ff6b6b',
                    backgroundColor: 'rgba(255, 107, 107, 0.1)',
                    borderWidth: 2,
                    fill: true,
                    tension: 0.35,
                    order: 4
                }, {
                    label: 'Measured AIF',
                    data: [],
                    borderColor: '#ff6b6b',
                    borderDash: [4, 3],
                    borderWidth: 1.5,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.1,
                    order: 3
                }, {
                    label: 'Tissue Response',
                    data: [],
                    borderColor: '#ffe66d',
                    backgroundColor: 'rgba(255, 230, 109, 0.15)',
                    borderWidth: 2.5,
                    fill: true,
                    tension: 0.35,
                    order: 2
                }, {
                    label: 'VOF (Venous Output)',
                    data: [],
                    borderColor: '#1dd1a1',
                    borderDash: [6, 4],
                    borderWidth: 2,
                    backgroundColor: 'transparent',
                    fill: false,
                    tension: 0.25,
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
                    x: { ...commonOptions.scales.x, max: 120 },
                    y: {
                        ...commonOptions.scales.y,
                        min: 0,
                        max: undefined,
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
                    x: { ...commonOptions.scales.x, max: 120 } // Extended to 120s
                }
            },
            plugins: [annotationPlugin] // Register local plugin
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
        state.noiseModel = noiseModelSelect.value;
        state.lambda = parseFloat(lambdaSlider.value);
        state.scanDuration = parseFloat(scanDurationSlider.value);
        state.model = modelSelect.value;

        // Update UI Text
        cbfValue.textContent = state.cbf;
        aifWidthValue.textContent = state.aifWidth;
        mttValue.textContent = state.mtt;
        tmaxValue.textContent = state.tmax;
        cardiacOutputValue.textContent = state.cardiacOutput;
        noiseValue.textContent = state.noise;
        lambdaValue.textContent = state.lambda.toFixed(2);
        scanDurationValue.textContent = state.scanDuration;

        // Calculate CBV = CBF * MTT / 60 (if CBF is per min and MTT in sec)
        // Usually CBV is ml/100g. CBF is ml/100g/min. MTT is sec.
        // CBV = (CBF / 60) * MTT
        const cbv = (state.cbf / 60) * state.mtt;
        cbvValue.textContent = cbv.toFixed(2);

        // Generate Data
        const timeSteps = Array.from({ length: state.timePoints }, (_, i) => i * state.dt);
        const aifData = generateAIF(timeSteps, state.aifWidth, state.cardiacOutput);
        const residueData = generateResidue(timeSteps, state.mtt, state.model, state.dt); // No delay
        const vofDelay = computeVofDelay(state.cardiacOutput) + state.tmax;
        const transportData = generateTransport(timeSteps, state.mtt, vofDelay, state.model, state.dt);

        // IRF Data (no delay in function, tmax handled in convolution)
        const irfData = generateIRF(timeSteps, state.cbf, state.mtt, state.model, state.dt);
        const shiftedIrfData = generateShiftedIRF(irfData, state.tmax, state.dt);

        const acquisitionSamples = getAcquisitionSamples(state.scanDuration, state.dt, state.timePoints);
        state.scanDurationEffective = (acquisitionSamples - 1) * state.dt;
        scanDurationEffectiveDisplay.textContent = state.scanDurationEffective.toFixed(0);

        // Calculate VOF = AIF ⊗ h(t) - Auto mass-conservative
        const vofData = convolveForVOF(aifData, transportData, state.dt);

        // Convolve for Tissue (with tmax applied to residue)
        const perfusionScale = getCardiacPerfusionScale(state.cardiacOutput);
        const useCircularDelay = state.deconvMode === 'insensitive';
        const tissueDataRaw = convolve(aifData, residueData, state.cbf, state.dt, state.tmax, perfusionScale, useCircularDelay);
        const tissueMeasuredFull = addNoise(tissueDataRaw, state.noise, state.noiseModel);
        const tissueObservedDisplay = tissueDataRaw.map((val, idx) => idx < acquisitionSamples ? tissueMeasuredFull[idx] : null);

        // Perform Deconvolution on (potentially noisy) tissue data
        const aifMeasuredFull = addNoise(aifData, state.noise, state.noiseModel);
        const tissueForDeconv = tissueMeasuredFull.slice(0, acquisitionSamples);
        const aifForDeconv = aifMeasuredFull.slice(0, acquisitionSamples);
        const deconvolvedIRF = deconvolveSVD(tissueForDeconv, aifForDeconv, state.dt, state.lambda, state.deconvMode, perfusionScale);
        const deconvolvedDisplay = applyAcquisitionWindow(deconvolvedIRF, deconvolvedIRF.length, state.timePoints);

        const aifObservedDisplay = applyAcquisitionWindow(aifMeasuredFull, acquisitionSamples, state.timePoints);
        const vofObservedDisplay = applyAcquisitionWindow(vofData, acquisitionSamples, state.timePoints);

        // Derive parameters from deconvolution
        const deconvParams = deriveParametersFromDeconvolution(deconvolvedIRF, state.dt, state.scanDurationEffective);

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
        flowChart.data.datasets[0].data = aifData;
        flowChart.data.datasets[1].data = aifObservedDisplay;
        flowChart.data.datasets[2].data = tissueObservedDisplay;
        flowChart.data.datasets[3].data = vofObservedDisplay;
        flowChart.update();

        residueChart.data.datasets[0].data = irfData;
        residueChart.data.datasets[1].data = shiftedIrfData;
        residueChart.data.datasets[2].data = deconvolvedDisplay; // Deconvolved IRF (windowed)
        // Update stepped property based on model
        residueChart.data.datasets[0].stepped = state.model === 'boxcar';
        residueChart.data.datasets[1].stepped = state.model === 'boxcar';
        residueChart.update();
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
    noiseModelSelect.addEventListener('change', updateSimulation);
    scanDurationSlider.addEventListener('input', updateHandler);

    // --- Initialization ---
    initCharts();
    // Initial simulation update
    updateSimulation();

    // Trigger MathJax typeset if available
    if (typeof MathJax !== 'undefined' && MathJax.typeset) {
        MathJax.typeset();
    }
});
