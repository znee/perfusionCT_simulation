// modules/perfusion-models.js

import { gaussianRandom, poissonRandom } from './utils.js';

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
export function generateAIF(timeSteps, widthScale, cardiacOutput = 5.0) {
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
export function computeVofDelay(cardiacOutput) {
    const normalCO = 5.0;
    const baseDelay = 4.0; // seconds for normal venous transit
    const coFactor = normalCO / Math.max(1.0, cardiacOutput);
    const delay = baseDelay * coFactor;
    return Math.min(12, Math.max(2, delay));
}

export function generateTransport(timeSteps, mtt, delay, model, dt) {
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
export function generateResidue(timeSteps, mtt, model, dt) {
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
export function convolve(aif, residue, cbf, dt, tmax, perfusionScale = 1, circular = false) {
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
export function convolveForVOF(aif, transport, dt) {
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
export function generateIRF(timeSteps, cbf, mtt, model, dt) {
    const residue = generateResidue(timeSteps, mtt, model, dt);
    return residue.map(r => r * cbf);
}

// Shifted IRF for visualization (shows effective IRF in tissue with Tmax delay)
// Shifted Residue Function for Forward Model
// Applies Tmax delay directly to the residue function R(t)
export function generateShiftedResidue(residue, tmax, dt) {
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
export function generateShiftedIRF(irfData, tmax, dt) {
    // Reuse the shifted residue logic but for IRF
    return generateShiftedResidue(irfData, tmax, dt);
}

export function addNoise(data, noiseLevel, model = 'gaussian') {
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
// Solves: C_t(t) = A * R(t) where A is Toeplitz or Circulant matrix from AIF
// Returns: Deconvolved IRF (scaled residue function)
export function deconvolveSVD(tissueCurve, aif, dt, lambda, mode, perfusionScale = 1) {
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
export function deriveParametersFromDeconvolution(deconvolvedIRF, dt, maxTime = Infinity) {
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
