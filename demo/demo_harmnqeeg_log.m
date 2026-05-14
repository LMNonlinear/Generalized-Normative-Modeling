%% demo_harmnqeeg_log.m
% End-to-end demo: train GNM on HarMNqEEG log data, then predict z-scores
% for out-of-sample Barbados1978 malnutrition subjects.
%
% Training data: T_global_log_1564.CSV (1564 subjects, 14 sites, 18 log responses)
% Predict data:  T_log_BBDSmalnu44.csv  (44 subjects from DEDAAS Barbados1978)
%
% NOTE: This demo references data that is NOT bundled with the public
%       toolbox. Update the paths below to point to your own copies of
%       the HarMNqEEG and Barbados1978 datasets before running. The
%       HarMNqEEG dataset is available from Hernandez-Gonzalez et al.
%       (2024); Barbados1978 access is governed by its data agreement.
%
% NOTE: The predict CSV must have factor columns (age, freq) on the SAME
%       scale as the training data. In this demo, training used log10(age),
%       so we transform the predict age here. For general use, users should
%       prepare their CSV with matching scales before calling gnm_predict.

close all; clc;
rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir);
run(fullfile(rootDir, 'setup.m'));

% --- USER: edit the four paths below ---
data_dir   = fullfile(rootDir, 'data');                  % directory containing the CSVs below
outdir     = fullfile(rootDir, 'result', 'HarMNqEEG_log');
path_train   = fullfile(data_dir, 'T_global_log_1564.CSV');
path_predict = fullfile(data_dir, 'T_log_BBDSmalnu44.csv');
% ----------------------------------------

if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ========================================================================
%  Step 1: Train
%  ========================================================================
result_path = fullfile(outdir, 'gnm_result.mat');

if isfile(result_path)
    fprintf('=== Step 1: Loading cached training result ===\n');
    load(result_path, 'mnhs');
    fprintf('  Loaded from: %s\n', result_path);
else
    fprintf('=== Step 1: Train GNM on HarMNqEEG log data ===\n');

    % Build response formula: 18 diagonal log power columns
    resp_names = arrayfun(@(i) sprintf('log%d_%d', i, i), 1:18, 'UniformOutput', false);
    resp_str = strjoin(resp_names, '+');
    formula = [resp_str ' ~ s(freq,age) | datasetDetail[functional(age)]'];

    mnhs = gnm(path_train, formula, ...
        'inform', {'name', 'country', 'device', 'sex', 'disease'}, ...
        'hRangeMu',    [0.4, 0.4; 0.4, 0.4], ...
        'num_hMu',     [1, 1], ...
        'hRangeSigma', [0.6, 0.6; 0.6, 0.6], ...
        'num_hSigma',  [1, 1], ...
        'tag', 'harmnqeeg_log');

    fprintf('  Starting fit...\n');
    tic;
    mnhs = gnm_fit(mnhs);
    t_elapsed = toc;
    fprintf('  Fit completed in %.1f s (%.1f min)\n', t_elapsed, t_elapsed / 60);

    % Save result
    save(result_path, 'mnhs', '-v7.3');
    fprintf('  Saved to: %s\n', result_path);
end

% Print training summary
z_col = mnhs.tregs{2}.opt.resp_z{1};
z_train = real(mnhs.T.(z_col));
z_clean = z_train(isfinite(z_train));
fprintf('  Training z-score: mean=%.4f, std=%.4f, N=%d\n', ...
    mean(z_clean), std(z_clean), numel(z_clean));

%% ========================================================================
%  Step 2: Predict z-scores for new data
%  ========================================================================
fprintf('\n=== Step 2: Predict z-scores for new data ===\n');

% --- Demo-specific preprocessing ---
% This predict CSV has age in raw years, but training used log10(age).
% Transform to match training scale. For general use, users should prepare
% their CSV so that factor columns match the training scale beforehand.
T_pred_raw = readtable(path_predict);
fprintf('  Predict age range (raw years): [%.2f, %.2f]\n', ...
    min(T_pred_raw.age), max(T_pred_raw.age));
T_pred_raw.age = log10(T_pred_raw.age);
fprintf('  Predict age range (log10):     [%.4f, %.4f]\n', ...
    min(T_pred_raw.age), max(T_pred_raw.age));

[~, pred_basename] = fileparts(path_predict);
path_predict_ready = fullfile(outdir, [pred_basename '_log10age.csv']);
writetable(T_pred_raw, path_predict_ready);

% --- General prediction (reusable for any prepared CSV) ---
mnhs.T = [];  % free training table memory
T_pred = gnm_predict(mnhs, path_predict_ready);

% Save prediction result (.mat)
pred_mat = fullfile(outdir, ['T_predict_' pred_basename '.mat']);
save(pred_mat, 'T_pred', '-v7.3');
fprintf('  Saved: %s\n', pred_mat);

% Export z-score + ystar columns to CSV (real part only for compatibility)
z_cols = mnhs.tregs{end}.opt.resp_z;
ystar_cols = append('ystar', mnhs.resp);
meta_cols = intersect({'name', 'age', 'freq', 'sex', 'disease', 'datasetDetail'}, ...
    T_pred.Properties.VariableNames, 'stable');
export_cols = [meta_cols, z_cols(:)', ystar_cols(:)'];
T_export = T_pred(:, export_cols);
for ic = 1:numel(z_cols)
    T_export.(z_cols{ic}) = real(T_export.(z_cols{ic}));
end
for ic = 1:numel(ystar_cols)
    T_export.(ystar_cols{ic}) = real(T_export.(ystar_cols{ic}));
end
csv_out = fullfile(outdir, ['T_predict_' pred_basename '_zscore.csv']);
writetable(T_export, csv_out);
fprintf('  Saved: %s\n', csv_out);

%% ========================================================================
%  Step 3: Summary statistics
%  ========================================================================
fprintf('\n=== Step 3: Prediction Summary ===\n');
n_subj = numel(unique(T_pred.name));
fprintf('  Subjects: %d\n', n_subj);
fprintf('  Rows: %d\n', height(T_pred));

% z-score statistics per response
fprintf('\n  Batch z-score summary (first 5 responses):\n');
fprintf('  %-30s %8s %8s %8s\n', 'Column', 'Mean', 'Std', 'N>2');
for ir = 1:min(5, numel(z_cols))
    zv = real(T_pred.(z_cols{ir}));
    zv = zv(isfinite(zv));
    fprintf('  %-30s %8.3f %8.3f %8d\n', z_cols{ir}, mean(zv), std(zv), sum(abs(zv) > 2));
end

fprintf('\n=== DONE ===\n');
