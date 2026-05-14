function mnhs = gnm_initialize_model(mnhs)
% gnm_initialize_model  Initialize a GNM model struct from configuration.
%
%   mnhs = gnm_initialize_model(mnhs)
%
%   Reads CSV data, validates configuration fields (resp, factor, batch,
%   inform), detects complex responses, sets default bandwidth parameters,
%   and prepares the model metadata struct for pipeline execution.
%
%   This function is called internally by gnm_build_default_pipeline.
%   Users typically do not call it directly.
%
%   See also: gnm_build_default_pipeline, gnm_fit

rng(1)
dbstop if error

%% initialize project information
if ~isfield(mnhs, 'mnhfun')
    mnhs.mnhfun = @gnm_compute_zscores;
end

st = dbstack;
if ~isfield(mnhs, 'tag')
    mnhs.tag = st(2).name;
end
if ~isfield(mnhs, 'batch') || isempty(mnhs.batch)
    error('GNM:config', ...
        'mnhs.batch is required. E.g. mnhs.batch = {''site''}');
end

%% ========================================================================
%  User-configurable fields: resp, factor, inform
%  ========================================================================

% resp (REQUIRED): response variable column names
if ~isfield(mnhs, 'resp') || isempty(mnhs.resp)
    error('GNM:config', ...
        'mnhs.resp is required. E.g. mnhs.resp = {''y1'',''y2''}');
end
if ischar(mnhs.resp)
    mnhs.resp = {mnhs.resp};
end
resp = mnhs.resp(:);

% factor (REQUIRED): regression factor column names
if ~isfield(mnhs, 'factor') || isempty(mnhs.factor)
    error('GNM:config', ...
        'mnhs.factor is required. E.g. mnhs.factor = {''age''}');
end
if ischar(mnhs.factor)
    mnhs.factor = {mnhs.factor};
end
factor = mnhs.factor(:);

% inform (optional): metadata columns to import
if ~isfield(mnhs, 'inform') || isempty(mnhs.inform)
    mnhs.inform = {'name'};
end
if ischar(mnhs.inform)
    mnhs.inform = {mnhs.inform};
end
inform = mnhs.inform(:);

if ~isfield(mnhs, 'clean_individual_outlier')
    mnhs.clean_individual_outlier = false;
end

%% define how many levels of model
mnhs.level(find_char(mnhs.level, 'batch')) = mnhs.batch;
level = mnhs.level;
num_level = length(level);

%% read data
if find_char(level, 'global')
    idglobal = find_char(level, 'global');
    levelraw = mnhs.level;
    level(idglobal) = [];
end
header2import = cat(1, inform(:), level(:), factor(:), resp(:), mnhs.batch(:));

% deduplicate header list (in case inform/factor/batch overlap)
header2import = unique(header2import, 'stable');

opts = detectImportOptions(mnhs.path_csv);
% remove any columns not present in CSV
available_vars = opts.VariableNames;
missing = setdiff(header2import, available_vars);
if ~isempty(missing)
    warning('GNM:config', 'Columns not found in CSV (skipped): %s', ...
        strjoin(missing, ', '));
    header2import = intersect(header2import, available_vars, 'stable');
end
opts.SelectedVariableNames = header2import;
mnhs.T = readtable(mnhs.path_csv, opts);

%% deal with complex cross part
mnhs.parts = {@real, @imag};
mnhs.num_resp = size(resp, 1);

%% Determine mode: training (called from gnm_build_default_pipeline) or predict
% In training mode, set up all options and detect complex responses.
% In predict mode (called from legacy predict path), handle reRefBatch only.
caller_name = st(2).name;
is_training = ~contains(caller_name, 'predict');

if is_training
    % complex responses: auto-detect from actual data
    % sep_complex: abs==0; complex==1; real==2; imagenary=3;
    mnhs.opt.sep_complex = 2 * ones(1, mnhs.num_resp); % default: all real
    for ir = 1:mnhs.num_resp
        if ~isreal(mnhs.T.(resp{ir}))
            mnhs.opt.sep_complex(ir) = 1; % complex
        end
    end
    mnhs.opt.idreal    = find(mnhs.opt.sep_complex == 2);
    mnhs.opt.idcomplex = find(mnhs.opt.sep_complex == 1);
    mnhs.opt.complex_method = 'seperate';

    %% smooth param
    mnhs.opt.standardize = false;
    mnhs.opt.unstandardize = false;
    mnhs.opt.maxibf = 1;
    mnhs.opt.y_keep_all = true; % will be updated after bandwidth config
    mnhs.opt.global_h = false;
    mnhs.opt.ACC = 2;
    mnhs.opt.isgpu = false;
    mnhs.opt.nufft_loop = false;
    mnhs.opt.y_gcv_dstd_mu = 1;     % 1-SE rule for bandwidth selection
    mnhs.opt.y_gcv_dstd_sigma = 1;  % 1-SE rule for variance bandwidth
    mnhs.opt.compact.y = false;
    mnhs.opt.compact.mu = false;
    mnhs.opt.compact.eps = true;
    mnhs.opt.compact.sigma = true;
    mnhs.opt.compact.musigma = false;
    mnhs.opt.compact.z = false;
    mnhs.opt.compact.ystar = false;

    %% bandwidth parameters: dynamic defaults based on factor count
    nf = length(mnhs.factor);
    if ~isfield(mnhs, 'hScale')
        mnhs.hScale = @logspace;
    end
    if ~isfield(mnhs, 'hRangeMu')
        mnhs.hRangeMu = repmat([0.2, 0.8], nf, 1);
    end
    if ~isfield(mnhs, 'num_hMu')
        mnhs.num_hMu = 5 * ones(1, nf);
    end
    if ~isfield(mnhs, 'hRangeSigma')
        mnhs.hRangeSigma = 1.5 * mnhs.hRangeMu;
    end
    if ~isfield(mnhs, 'num_hSigma')
        mnhs.num_hSigma = mnhs.num_hMu;
    end

    % When multiple bandwidths configured, let GCV select the best one
    if max(mnhs.num_hMu) > 1 || max(mnhs.num_hSigma) > 1
        mnhs.opt.y_keep_all = false;
    end
else
    if ~isempty(mnhs.reRefBatch)
        disp(['use ', mnhs.reRefBatch{:,1}, ' as the reference of ', mnhs.reRefBatch{:,2}, ', replace information in T for calculation'])
        for i = 1:size(mnhs.reRefBatch, 1)
            mnhs.T = asnarray2table(mnhs.T, {mnhs.batch{1}, mnhs.reRefBatch{i,1}}, 'orignalBatch', mnhs.reRefBatch(i,1));
            mnhs.T = asnarray2table(mnhs.T, {mnhs.batch{1}, mnhs.reRefBatch{i,1}}, mnhs.batch, mnhs.reRefBatch(i,2));
        end
    end
end

if idglobal
    level = levelraw;
    mnhs.T = initializeTableVar(mnhs.T, 'global', [], {'used'});
end

%% save information
mnhs.inform = inform;
mnhs.factor = factor;
mnhs.resp = resp;
mnhs.num_level = num_level;

% infer num_obs_per_subject from data
subject_col = mnhs.inform{1};
con_name = unique(mnhs.T.(subject_col));
mnhs.num_obs_per_subject = size(get_subtable(mnhs.T, {subject_col, con_name(1)}, resp(:)), 1);

mnhs.resp_ystar = append('ystar', resp);

mnhs.batchCat = cell(mnhs.num_level, 1);
mnhs.num_batch = zeros(mnhs.num_level, 1);
for i_batch = 1:mnhs.num_level
    mnhs.batchCat{i_batch} = unique(mnhs.T.(level{i_batch}));
    mnhs.num_batch(i_batch) = numel(mnhs.batchCat{i_batch});
end

mnhs.cmsite = cell(num_level, 1);
for i = 1:num_level
    if mnhs.num_batch(i) > 1
        mnhs.cmsite{i} = spring(mnhs.num_batch(i));
    else
        mnhs.cmsite{i} = copper(1);
    end
end

%% outlier dimension
mnhs.cluster_dim = 2;

end
