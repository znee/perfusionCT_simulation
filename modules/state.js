// modules/state.js

export const state = {
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
