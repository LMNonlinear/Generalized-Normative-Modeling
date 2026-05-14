%% demo_batchModel
% Demonstrates how to configure the batch effect layer via mnhs.batchModel.
%
% The batch layer supports two modes for both mean and variance:
%   'functional' — batch effect varies with specified factors (default: age)
%   'constant'   — batch effect is a simple per-batch constant
%
% Constraint: when batch mean is 'constant', batch variance is
% automatically forced to 'constant' as well.
%
% If mnhs.batchModel is not provided, the default behavior is preserved:
%   batch mean = functional(age), batch variance = functional(age).

%% Example 1: Default — functional batch mean & variance (backward compatible)
% No batchModel field needed; equivalent to the original hardcoded behavior.
mnhs1 = struct();
mnhs1.path_csv = 'path/to/data.csv';   % <-- replace with actual path
% mnhs1 = gnm_build_default_pipeline(mnhs1);
% mnhs1 = gnm_run_pipeline(mnhs1);
disp('--- Example 1: Default (functional/functional) ---');
disp('tregs{2}.factor_mu  = {''age''}');
disp('tregs{2}.factor_sigma = {''age''}');

%% Example 2: Constant batch mean & variance
% When batch mean is constant, variance is automatically constant too.
mnhs2 = struct();
mnhs2.path_csv = 'path/to/data.csv';   % <-- replace with actual path
mnhs2.batchModel = struct( ...
    'meanType', 'constant', ...
    'varType',  'constant' ...
);
% mnhs2 = gnm_build_default_pipeline(mnhs2);
% mnhs2 = gnm_run_pipeline(mnhs2);
disp('--- Example 2: Constant batch effect ---');
disp('tregs{2}.factor_mu  = {}  (constant)');
disp('tregs{2}.factor_sigma = {}  (constant)');

%% Example 3: Functional batch mean, constant batch variance
mnhs3 = struct();
mnhs3.path_csv = 'path/to/data.csv';   % <-- replace with actual path
mnhs3.batchModel = struct( ...
    'meanType',    'functional', ...
    'meanFactors', {{'age'}}, ...
    'varType',     'constant' ...
);
% mnhs3 = gnm_build_default_pipeline(mnhs3);
% mnhs3 = gnm_run_pipeline(mnhs3);
disp('--- Example 3: Functional mean / constant variance ---');
disp('tregs{2}.factor_mu  = {''age''}');
disp('tregs{2}.factor_sigma = {}  (constant)');

%% Example 4: Constraint enforcement — constant mean overrides functional variance
mnhs4 = struct();
mnhs4.path_csv = 'path/to/data.csv';   % <-- replace with actual path
mnhs4.batchModel = struct( ...
    'meanType', 'constant', ...
    'varType',  'functional' ...   % will be forced to 'constant' with a warning
);
% mnhs4 = gnm_build_default_pipeline(mnhs4);
disp('--- Example 4: Constraint enforcement ---');
disp('meanType=constant + varType=functional => warning + varType forced to constant');
