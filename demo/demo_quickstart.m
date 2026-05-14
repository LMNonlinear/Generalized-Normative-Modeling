%% demo_quickstart.m - GNM-ToolBox Quick Start Demo
%
% This script demonstrates the core GNM workflow on a small synthetic dataset.
% It takes ~5 seconds to run and requires no external data downloads.
%
% Steps:
%   1. Load synthetic multi-site data
%   2. Fit a GNM model with constant batch effects
%   3. Examine z-score calibration
%   4. Predict on the same data (for demonstration)
%
% See also: gnm, gnm_fit, gnm_predict

%% Setup
clear; clc; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup.m'));

fprintf('\n=== GNM-ToolBox Quick Start Demo ===\n\n');

%% Step 1: Load synthetic data
% The test dataset has columns: name, site, age, y1, y2

data_path = fullfile(fileparts(mfilename('fullpath')), 'test_synthetic_data.csv');
if ~exist(data_path, 'file')
    error('Demo data not found: %s\nRun from the demo/ directory.', data_path);
end

fprintf('Data: %s\n', data_path);
T = readtable(data_path, 'VariableNamingRule', 'preserve');
fprintf('  Rows: %d, Columns: %d\n', height(T), width(T));
fprintf('  Sites: %s\n', strjoin(unique(T.site), ', '));
fprintf('  Age range: [%.3f, %.3f]\n\n', min(T.age), max(T.age));

%% Step 2: Fit GNM model
% Formula: y1 response, smooth function of age, site as batch with constant shift
formula = 'y1 ~ s(age) | site[constant]';
fprintf('Formula: %s\n', formula);

tic;
mnhs = gnm(data_path, formula, ...
    'inform', {'name'}, ...           % informational columns (not used in model)
    'tag', 'quickstart');             % project tag
mnhs = gnm_fit(mnhs);
t_fit = toc;

fprintf('\nFit completed in %.1f seconds.\n\n', t_fit);

%% Step 3: Examine results
% Access z-scores from the trained model
z_col = mnhs.tregs{end}.opt.resp_z{1};
z_all = real(mnhs.T.(z_col));
z_all = z_all(isfinite(z_all));

fprintf('=== Z-Score Calibration ===\n');
fprintf('  Mean:     %.4f (ideal: 0)\n', mean(z_all));
fprintf('  Std:      %.4f (ideal: 1)\n', std(z_all));
fprintf('  Skewness: %.4f\n', skewness(z_all));
fprintf('  Kurtosis: %.4f (excess)\n\n', kurtosis(z_all) - 3);

% Per-site calibration
sites = mnhs.T.site;
u_sites = unique(sites);
fprintf('  Per-site z-score means:\n');
for i = 1:numel(u_sites)
    mask = strcmp(sites, u_sites{i});
    zi = real(mnhs.T.(z_col)(mask));
    zi = zi(isfinite(zi));
    fprintf('    %-15s mean=%+.4f, std=%.4f (n=%d)\n', ...
        u_sites{i}, mean(zi), std(zi), numel(zi));
end

%% Step 4: Predict on new data (using same file as demo)
fprintf('\n=== Prediction Demo ===\n');
fprintf('Predicting on the same data (for demonstration)...\n');

% In practice, you would use a different CSV with new subjects
mnhs_for_predict = mnhs;
mnhs_for_predict.T = [];  % clear training data to save memory

T_pred = gnm_predict(mnhs_for_predict, data_path);

z_pred = real(T_pred.(z_col));
z_pred = z_pred(isfinite(z_pred));
fprintf('\n  Prediction z-scores: mean=%+.4f, std=%.4f\n', mean(z_pred), std(z_pred));

%% Step 5: Visualization (optional)
try
    figure('Position', [100, 100, 900, 400], 'Color', 'w');

    % Panel A: Z-score histogram
    subplot(1, 2, 1);
    histogram(z_all, 30, 'Normalization', 'pdf', ...
        'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none');
    hold on;
    x = linspace(-4, 4, 200);
    plot(x, normpdf(x), 'r-', 'LineWidth', 1.5);
    hold off;
    xlabel('z-score'); ylabel('Density');
    title('Z-Score Distribution vs N(0,1)');
    legend({'GNM z-scores', 'N(0,1)'}, 'Location', 'northeast');
    xlim([-4, 4]);
    grid on; box on;

    % Panel B: Per-site box plot
    subplot(1, 2, 2);
    valid = isfinite(real(mnhs.T.(z_col)));
    boxplot(real(mnhs.T.(z_col)(valid)), mnhs.T.site(valid));
    yline(0, '--r');
    xlabel('Site'); ylabel('z-score');
    title('Per-Site Z-Score Calibration');
    grid on; box on;

    sgtitle('GNM Quick Start Demo Results', 'FontWeight', 'bold');
    fprintf('  Figure displayed.\n');
catch ME
    fprintf('  (Visualization skipped: %s)\n', ME.message);
end

fprintf('\n=== Demo Complete ===\n');
fprintf('Next steps:\n');
fprintf('  - Try different formulas (e.g., functional batch: site[functional(age)])\n');
fprintf('  - Use your own data (CSV with age, site, and response columns)\n');
fprintf('  - See demo_harmnqeeg_log.m for a real-world EEG example\n\n');
