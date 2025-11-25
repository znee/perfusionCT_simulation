// modules/charts.js

import { state } from './state.js';

let flowChart, residueChart;

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

export function initCharts() {
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
                label: 'IRF',
                data: [],
                borderColor: '#4ecdc4',
                backgroundColor: 'rgba(78, 205, 196, 0.1)',
                borderWidth: 2,
                fill: true,
                stepped: true,
                order: 3
            }, {
                label: 'Deconvolved IRF (no noise)',
                data: [],
                borderColor: '#ec4899',
                borderDash: [5, 5],
                backgroundColor: 'transparent',
                borderWidth: 2,
                pointRadius: 0,
                fill: false,
                tension: 0.2,
                order: 2
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

export function updateCharts(data) {
    flowChart.data.datasets[0].data = data.aifData;
    flowChart.data.datasets[1].data = data.aifObservedDisplay;
    flowChart.data.datasets[2].data = data.tissueObservedDisplay;
    flowChart.data.datasets[3].data = data.vofObservedDisplay;
    flowChart.update();

    residueChart.data.datasets[0].data = data.shiftedIrfData;
    residueChart.data.datasets[1].data = data.deconvolvedNoNoiseDisplay; // Deconvolved IRF without noise
    residueChart.data.datasets[2].data = data.deconvolvedDisplay; // Deconvolved IRF (with noise)
    // Update stepped property based on model
    residueChart.data.datasets[0].stepped = state.model === 'boxcar';
    residueChart.update();
}
