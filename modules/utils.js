// modules/utils.js

export function gaussianRandom() {
    let u = 0, v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

export function poissonRandom(mean) {
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

export function getAcquisitionSamples(scanDuration, dt, maxPoints) {
    const clamped = Math.max(0, scanDuration);
    const samples = Math.floor(clamped / dt) + 1;
    return Math.min(maxPoints, Math.max(1, samples));
}

export function applyAcquisitionWindow(data, samples, totalLength) {
    const output = new Array(totalLength).fill(null);
    const limit = Math.min(totalLength, data.length);
    for (let i = 0; i < limit; i++) {
        if (i < samples) {
            output[i] = data[i];
        }
    }
    return output;
}

export function getCardiacPerfusionScale(cardiacOutput) {
    const normalCO = 5.0;
    const ratio = Math.max(0.2, cardiacOutput / normalCO);
    return Math.min(1.4, Math.max(0.6, Math.pow(ratio, 0.9)));
}

export function debounce(func, wait) {
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
