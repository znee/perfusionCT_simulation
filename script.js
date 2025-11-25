// script.js

import { state } from './modules/state.js';
import { initCharts, updateCharts } from './modules/charts.js';
import { updateUI, updateResults, setupEventListeners } from './modules/ui.js';
import {
    generateAIF,
    generateResidue,
    computeVofDelay,
    generateTransport,
    generateIRF,
    generateShiftedIRF,
    convolve,
    convolveForVOF,
    addNoise,
    deconvolveSVD,
    deriveParametersFromDeconvolution
} from './modules/perfusion-models.js';
import { getAcquisitionSamples, applyAcquisitionWindow, getCardiacPerfusionScale } from './modules/utils.js';

document.addEventListener('DOMContentLoaded', () => {
    // Check if numeric.js is loaded
    if (typeof numeric === 'undefined') {
        console.error('CRITICAL ERROR: numeric.js library not loaded! Deconvolution will fail.');
        alert('Error: numeric.js library not loaded. Please refresh the page.');
        return;
    }

    function updateSimulation() {
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

        // Calculate VOF = AIF ⊗ h(t) - Auto mass-conservative
        const vofData = convolveForVOF(aifData, transportData, state.dt);

        // Convolve for Tissue (with tmax applied to residue)
        const perfusionScale = getCardiacPerfusionScale(state.cardiacOutput);
        const useCircularDelay = state.deconvMode === 'insensitive';
        const tissueDataRaw = convolve(aifData, residueData, state.cbf, state.dt, state.tmax, perfusionScale, useCircularDelay);
        const tissueMeasuredFull = addNoise(tissueDataRaw, state.noise, state.noiseModel);
        const tissueObservedDisplay = tissueMeasuredFull.map((val, idx) => idx < acquisitionSamples ? val : null);

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

        // Update UI and Charts
        updateUI();
        updateResults(deconvParams);
        updateCharts({
            aifData,
            aifObservedDisplay,
            tissueObservedDisplay,
            vofObservedDisplay,
            irfData,
            shiftedIrfData,
            deconvolvedDisplay
        });
    }

    // --- Initialization ---
    initCharts();
    setupEventListeners(updateSimulation);
    updateSimulation();

    // Trigger MathJax typeset if available
    if (typeof MathJax !== 'undefined' && MathJax.typeset) {
        MathJax.typeset();
    }
});
