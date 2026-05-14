function T_new = gnm_predict(mnhs_trained, path_csv_new, varargin)
% gnm_predict  Predict batch-corrected z-scores for new data using a trained GNM model.
%
%   T_new = gnm_predict(mnhs_trained, path_csv_new)
%   T_new = gnm_predict(mnhs_trained, path_csv_new, 'Name', Value, ...)
%
%   Applies a trained GNM normative model to new (out-of-sample) data.
%   Produces z-scores, ystar (harmonized responses), and intermediate
%   mu/sigma estimates for each hierarchical level.
%
%   The input CSV must contain the same response, factor, and batch columns
%   as the training data. If the new data comes from a batch not seen during
%   training, use the 'reRefBatch' option to map it to an existing batch.
%
%   Inputs:
%     mnhs_trained - trained GNM struct (from gnm_fit or loaded .mat file).
%                    Required fields: tregs, resp, opt, factor, batch, level, inform.
%     path_csv_new - path to CSV file with new observations.
%
%   Name-Value parameters:
%     'reRefBatch' - N x 2 cell array mapping new batch names to trained
%                    batch names. Format: {'newBatch1', 'trainedBatch1'; ...}
%                    Use when new data comes from a site not in the training set.
%     'tag'        - project identifier tag (default: 'predict')
%
%   Output:
%     T_new - table with original data plus predicted columns:
%               z-scores:  z{level}{resp}       (batch-corrected z ~ N(0,1))
%               ystar:     ystar{resp}           (harmonized response)
%               mu:        mu{level}{resp}       (predicted mean)
%               musigma:   musigma{level}{resp}  (predicted variance of mean)
%
%   Examples:
%     % Train a model
%     mnhs = gnm('train.csv', 'y1+y2 ~ s(freq,age) | site[functional(age)]');
%     mnhs = gnm_fit(mnhs);
%     save('norm.mat', 'mnhs');
%
%     % Predict on new data (same site exists in training)
%     trained = load('norm.mat'); trained = trained.mnhs;
%     T_new = gnm_predict(trained, 'new_data.csv');
%
%     % Predict on new site (map to closest training site)
%     T_new = gnm_predict(trained, 'new_data.csv', ...
%         'reRefBatch', {'NewSite_BrainAmp', 'BrainProduct Germany'});
%
%   See also: gnm, gnm_fit, gnm_link, tnureg_zmap_predict

%% Parse inputs
p = inputParser;
p.addRequired('mnhs_trained', @isstruct);
p.addRequired('path_csv_new', @(x) ischar(x) || isstring(x));
p.addParameter('reRefBatch', {}, @(x) iscell(x) && (isempty(x) || size(x, 2) == 2));
p.addParameter('tag', 'predict', @(x) ischar(x) || isstring(x));
p.parse(mnhs_trained, path_csv_new, varargin{:});

reRefBatch = p.Results.reRefBatch;
tag = char(p.Results.tag);

%% Validate trained model (#5: include 'inform' in required fields)
required_fields = {'tregs', 'resp', 'opt', 'factor', 'batch', 'level', 'inform'};
for i = 1:numel(required_fields)
    assert(isfield(mnhs_trained, required_fields{i}), ...
        'GNM:gnm_predict', 'Missing field "%s" in trained model.', required_fields{i});
end

path_csv_new = char(path_csv_new);
assert(isfile(path_csv_new), 'GNM:gnm_predict', 'CSV file not found: %s', path_csv_new);

fprintf('\n=== GNM Predict ===\n');
fprintf('  Model : %d levels, %d responses\n', numel(mnhs_trained.tregs), numel(mnhs_trained.resp));
fprintf('  Data  : %s\n', path_csv_new);

%% Read new data
% Build header list for CSV import
level = mnhs_trained.level;
idglobal = find_char(level, 'global');
level_import = level;
if idglobal
    level_import(idglobal) = [];
end
header2import = cat(1, mnhs_trained.inform(:), level_import(:), ...
    mnhs_trained.factor(:), mnhs_trained.resp(:), mnhs_trained.batch(:));
header2import = unique(header2import, 'stable');

csv_opts = detectImportOptions(path_csv_new);
available_vars = csv_opts.VariableNames;
missing_cols = setdiff(header2import, available_vars);

% (#5) Fail fast on missing critical columns (resp, factor, batch)
critical_cols = [mnhs_trained.resp(:); mnhs_trained.factor(:); mnhs_trained.batch(:)];
missing_critical = setdiff(critical_cols, available_vars);
if ~isempty(missing_critical)
    error('GNM:gnm_predict:missingColumns', ...
        'Required columns not found in CSV: %s', strjoin(missing_critical, ', '));
end

% Warn about non-critical missing columns (inform)
missing_noncritical = setdiff(missing_cols, missing_critical);
if ~isempty(missing_noncritical)
    warning('GNM:gnm_predict', 'Optional columns not found in CSV (skipped): %s', ...
        strjoin(missing_noncritical, ', '));
    header2import = intersect(header2import, available_vars, 'stable');
end

csv_opts.SelectedVariableNames = header2import;
T_new = readtable(path_csv_new, csv_opts);
fprintf('  Loaded: %d rows x %d cols\n', size(T_new, 1), size(T_new, 2));

% Initialize global column (same as gnm_initialize_model)
if idglobal
    T_new = initializeTableVar(T_new, 'global', [], {'used'});
end

%% Handle reRefBatch: remap new batch names to trained ones (#3: per-row save/restore)
if ~isempty(reRefBatch)
    batch_col = mnhs_trained.batch{1};
    % Save original batch names before remapping
    T_new.originalBatch = T_new.(batch_col);
    for i = 1:size(reRefBatch, 1)
        new_name = reRefBatch{i, 1};
        trained_name = reRefBatch{i, 2};
        fprintf('  [reRefBatch] %s -> %s\n', new_name, trained_name);
        % Replace batch column value for matching rows
        mask = strcmp(T_new.(batch_col), new_name);
        T_new.(batch_col)(mask) = {trained_name};
    end
end

%% Apply link forward transform (if training used one)
has_link = isfield(mnhs_trained, 'link') && ~strcmp(mnhs_trained.link.name, 'identity');
if has_link
    link = mnhs_trained.link;

    % Link functions are undefined for complex responses
    if any(mnhs_trained.opt.sep_complex == 1)
        error('GNM:gnm_predict', ...
            'Link function "%s" is not supported for complex responses.', link.name);
    end

    for i = 1:numel(mnhs_trained.resp)
        col = mnhs_trained.resp{i};
        if ismember(col, T_new.Properties.VariableNames)
            y_col = T_new.(col);
            valid = link.domain_check(y_col);
            valid(isnan(y_col)) = true;
            assert(all(valid, 'all'), 'GNM:gnm_predict:domainViolation', ...
                'Link "%s" domain violation in column "%s".', link.name, col);
            T_new.(col) = link.forward(y_col);
        end
    end
    fprintf('  [link] Applied forward transform: %s\n', link.name);
end

%% Core prediction: tnureg_zmap_predict
opt_pred = mnhs_trained.opt;
% Enable compaction of heavy intermediate outputs to save memory
opt_pred.compact.eps = true;
opt_pred.compact.sigma = true;
opt_pred.compact.musigma = true;
opt_pred.compact.mu = true;        % compact mu — recomputed per-response for ystar
opt_pred.compact.z = false;        % keep z (final output)
opt_pred.compact.ystar = false;    % keep ystar (final output)
% Configure computation flags for prediction
for il = 1:numel(mnhs_trained.tregs)
    mnhs_trained.tregs{il}.opt.calc_eps = false;
    mnhs_trained.tregs{il}.opt.calc_sigma = false;
    % calc_musigma: needed internally for z computation,
    % but never saved to table in predict mode — recomputed for ystar
    mnhs_trained.tregs{il}.opt.calc_musigma = true;
    mnhs_trained.tregs{il}.opt.save_musigma = false;
end

T_new = tnureg_zmap_predict(T_new, mnhs_trained.resp, mnhs_trained.tregs, opt_pred);

%% Apply link inverse transform on ystar and global mu
if has_link
    link = mnhs_trained.link;
    resp_ystar = append('ystar', mnhs_trained.resp);

    % Inverse transform ystar
    for i = 1:numel(resp_ystar)
        col = resp_ystar{i};
        if ismember(col, T_new.Properties.VariableNames)
            T_new.(col) = link.inverse(T_new.(col));
        end
    end

    % Inverse transform global mu (level 1) for interpretability
    if ~opt_pred.compact.mu
        mu_cols = mnhs_trained.tregs{1}.opt.resp_mu;
        for i = 1:numel(mu_cols)
            col = mu_cols{i};
            if ismember(col, T_new.Properties.VariableNames)
                T_new.(col) = link.inverse(T_new.(col));
            end
        end
    end

    fprintf('  [link] Applied inverse transform: %s (ystar, global mu)\n', link.name);
end

%% Restore original batch names (#3: per-row restoration)
if ~isempty(reRefBatch)
    batch_col = mnhs_trained.batch{1};
    if ismember('originalBatch', T_new.Properties.VariableNames)
        % Restore original batch names from saved column
        T_new.(batch_col) = T_new.originalBatch;
        % Remove helper column
        T_new.originalBatch = [];
    end
end

%% Summary
n_resp = numel(mnhs_trained.resp);
z_col = mnhs_trained.tregs{end}.opt.resp_z{1};
z_vals = real(T_new.(z_col));
z_clean = z_vals(isfinite(z_vals));
fprintf('  Results: %d rows, %d responses\n', size(T_new, 1), n_resp);
fprintf('  z-score (batch-corrected): mean=%.4f, std=%.4f\n', mean(z_clean), std(z_clean));
fprintf('===================\n\n');

end
