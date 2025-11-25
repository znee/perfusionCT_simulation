// modules/ui.js

import { state } from './state.js';
import { debounce } from './utils.js';

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

export function updateUI() {
    cbfValue.textContent = state.cbf;
    aifWidthValue.textContent = state.aifWidth;
    mttValue.textContent = state.mtt;
    tmaxValue.textContent = state.tmax;
    cardiacOutputValue.textContent = state.cardiacOutput;
    noiseValue.textContent = state.noise;
    lambdaValue.textContent = state.lambda.toFixed(2);
    scanDurationValue.textContent = state.scanDuration;
    scanDurationEffectiveDisplay.textContent = state.scanDurationEffective.toFixed(0);

    const cbv = (state.cbf / 60) * state.mtt;
    cbvValue.textContent = cbv.toFixed(2);
}

export function updateResults(deconvParams) {
    const cbv = (state.cbf / 60) * state.mtt;

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
}

export function setupEventListeners(updateSimulation) {
    const isMobile = window.matchMedia('(max-width: 640px)').matches;
    const updateHandler = isMobile ? debounce(updateSimulation, 150) : updateSimulation;

    cbfSlider.addEventListener('input', () => {
        state.cbf = parseFloat(cbfSlider.value);
        updateHandler();
    });

    aifWidthSlider.addEventListener('input', () => {
        state.aifWidth = parseFloat(aifWidthSlider.value);
        updateHandler();
    });

    mttSlider.addEventListener('input', () => {
        state.mtt = parseFloat(mttSlider.value);
        updateHandler();
    });

    tmaxSlider.addEventListener('input', () => {
        state.tmax = parseFloat(tmaxSlider.value);
        updateHandler();
    });

    cardiacOutputSlider.addEventListener('input', () => {
        state.cardiacOutput = parseFloat(cardiacOutputSlider.value);
        updateHandler();
    });

    noiseSlider.addEventListener('input', () => {
        state.noise = parseFloat(noiseSlider.value);
        updateHandler();
    });

    lambdaSlider.addEventListener('input', () => {
        state.lambda = parseFloat(lambdaSlider.value);
        updateHandler();
    });

    scanDurationSlider.addEventListener('input', () => {
        state.scanDuration = parseFloat(scanDurationSlider.value);
        updateHandler();
    });

    modelSelect.addEventListener('change', (e) => {
        state.model = e.target.value;
        updateSimulation();
    });

    deconvModeSelect.addEventListener('change', (e) => {
        state.deconvMode = e.target.value;
        updateSimulation();
    });

    noiseModelSelect.addEventListener('change', (e) => {
        state.noiseModel = e.target.value;
        updateSimulation();
    });
}
