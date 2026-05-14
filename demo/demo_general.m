%% demo_general
% Demonstrates generic normative modeling with real-valued responses.
%
% This demo shows the minimal configuration needed to run GNM-ToolBox
% on arbitrary tabular data (e.g. cortical thickness, brain volume).
%
% Required user inputs:
%   mnhs.path_csv  - path to CSV data file
%   mnhs.resp      - cell array of response variable column names
%   mnhs.factor    - cell array of regression factor column names
%
% Optional user inputs:
%   mnhs.batch     - batch variable column names (default: {'study'})
%   mnhs.inform    - metadata column names       (default: {'name'})
%   mnhs.batchModel - batch layer configuration  (default: functional/functional)

%% Example 1: Single-factor real-valued responses (e.g. cortical thickness)
mnhs1 = struct();
mnhs1.path_csv = 'path/to/cortical_thickness.csv';  % <-- replace with actual path

% REQUIRED: specify response and factor columns from your CSV
mnhs1.resp   = {'thickness_region1', 'thickness_region2', 'thickness_region3'};
mnhs1.factor = {'age'};

% OPTIONAL: metadata and batch
mnhs1.batch  = {'site'};
mnhs1.inform = {'name', 'sex'};

% OPTIONAL: bandwidth parameters (1 factor => 1-row matrix)
mnhs1.hRangeMu    = [0.4, 0.4];
mnhs1.num_hMu     = 1;
mnhs1.hRangeSigma = [0.6, 0.6];
mnhs1.num_hSigma  = 1;

% OPTIONAL: constant batch effect (no factor dependence)
mnhs1.batchModel = struct('meanType', 'constant', 'varType', 'constant');

% Uncomment to run:
% mnhs1 = gnm_build_default_pipeline(mnhs1);
% mnhs1 = gnm_run_pipeline(mnhs1);

disp('--- Example 1: Single-factor, real-valued, constant batch ---');
disp('  mnhs.resp   = {''thickness_region1'', ...}');
disp('  mnhs.factor = {''age''}');
disp('  mnhs.batchModel.meanType = ''constant''');

%% Example 2: Multi-factor responses (e.g. brain volume ~ age + ICV)
mnhs2 = struct();
mnhs2.path_csv = 'path/to/brain_volume.csv';  % <-- replace with actual path

mnhs2.resp   = {'vol_hippocampus', 'vol_amygdala', 'vol_thalamus'};
mnhs2.factor = {'age', 'icv'};   % 2 factors => 2-row bandwidth matrices

mnhs2.batch  = {'scanner'};
mnhs2.inform = {'name', 'sex', 'diagnosis'};

% 2 factors => 2-row bandwidth matrices
mnhs2.hRangeMu    = [0.4, 0.4; 0.4, 0.4];
mnhs2.num_hMu     = [1, 1];
mnhs2.hRangeSigma = [0.6, 0.6; 0.6, 0.6];
mnhs2.num_hSigma  = [1, 1];

% functional batch effect (default behavior)
mnhs2.batchModel = struct( ...
    'meanType',    'functional', ...
    'meanFactors', {{'age'}}, ...       % batch mean varies with age only
    'varType',     'constant' ...       % batch variance is constant
);

% Uncomment to run:
% mnhs2 = gnm_build_default_pipeline(mnhs2);
% mnhs2 = gnm_run_pipeline(mnhs2);

disp('--- Example 2: Multi-factor, functional batch mean ---');
disp('  mnhs.resp   = {''vol_hippocampus'', ...}');
disp('  mnhs.factor = {''age'', ''icv''}');
disp('  mnhs.batchModel.meanType = ''functional'' (age)');

%% Backward compatibility: EEG users migrating from old interface
% Old EEG scripts only need to add these two lines:
%   mnhs.resp   = get_spec_hearder([], 'riemlogm');  % call externally
%   mnhs.factor = {'freq', 'age'};
% All other behavior is preserved.

%% ========================================================================
%  Formula interface (alternative syntax)
%  ========================================================================
% The gnm() function provides a compact formula-based interface.
% It builds the same mnhs struct as the manual approach above.

%% Example 3: Formula — single factor, constant batch
% Equivalent to Example 1 above
% mnhs3 = gnm('path/to/cortical_thickness.csv', ...
%     'thickness_region1+thickness_region2+thickness_region3 ~ s(age) | site[constant]', ...
%     'inform', {'name', 'sex'}, ...
%     'hRangeMu', [0.4, 0.4], 'num_hMu', 1, ...
%     'hRangeSigma', [0.6, 0.6], 'num_hSigma', 1);
% mnhs3 = gnm_fit(mnhs3);

disp('--- Example 3: Formula interface, constant batch ---');
disp('  gnm(csv, ''y1+y2 ~ s(age) | site[constant]'')');

%% Example 4: Formula — multi-factor, functional batch
% Equivalent to Example 2 above
% mnhs4 = gnm('path/to/brain_volume.csv', ...
%     'vol_hippocampus+vol_amygdala+vol_thalamus ~ s(age, icv) | scanner[functional(age)]', ...
%     'inform', {'name', 'sex', 'diagnosis'}, ...
%     'hRangeMu', [0.4, 0.4; 0.4, 0.4], 'num_hMu', [1, 1], ...
%     'hRangeSigma', [0.6, 0.6; 0.6, 0.6], 'num_hSigma', [1, 1]);
% mnhs4 = gnm_fit(mnhs4);

disp('--- Example 4: Formula interface, functional batch ---');
disp('  gnm(csv, ''y ~ s(age,icv) | scanner[functional(age)]'')');

%% Example 5: Formula — no batch variable
% When no batch is specified, gnm_initialize_model uses its default batch = {'study'}
% mnhs5 = gnm('path/to/data.csv', 'y1 ~ s(age)');
% mnhs5 = gnm_fit(mnhs5);

disp('--- Example 5: Formula interface, no batch ---');
disp('  gnm(csv, ''y1 ~ s(age)'')  % uses default batch');
