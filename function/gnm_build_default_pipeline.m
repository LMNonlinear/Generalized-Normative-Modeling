function mnhs = gnm_build_default_pipeline(mnhs)
% gnm_build_default_pipeline  Build the default two-level GNM pipeline.
%
%   mnhs = gnm_build_default_pipeline(mnhs)
%
%   Configures a global + batch hierarchical normative model:
%     Level 1: global normative trajectory (all subjects pooled)
%     Level 2: batch-specific effects (per-site correction)
%
%   This function:
%     1. Sets the pipeline function to @gnm_compute_zscores
%     2. Calls gnm_initialize_model() to load data and set defaults
%     3. Builds the tregs (regression specification) for both levels
%
%   The batch layer supports two modes (configured via mnhs.batchModel):
%     'functional' — batch effect varies with factors (default)
%     'constant'   — batch effect is a simple per-batch constant
%
%   See also: gnm_fit, gnm_initialize_model, gnm_compute_zscores

mnhs.mnhfun = @gnm_compute_zscores;
mnhs.level = {'global', 'batch'};
mnhs = gnm_initialize_model(mnhs);

%% define model
tregs = cell(mnhs.num_level, 1);

%% Level 1: global normative model (mu-sigma-z)
tregs{1,1}.condition        = {mnhs.level{1}; unique(mnhs.T.(mnhs.level{1}))};
tregs{1,1}.resp             = append(mnhs.level{1}, mnhs.resp);
tregs{1,1}.factor_mu        = mnhs.factor;
tregs{1,1}.factor_sigma     = mnhs.factor;
ifactor                     = find_char(mnhs.factor, unique([tregs{1,1}.factor_mu, tregs{1,1}.factor_sigma], 'stable'));
tregs{1,1}.hList            = {getHlist(mnhs.num_hMu(ifactor), mnhs.hRangeMu(ifactor,:), mnhs.hScale); ...
    getHlist(mnhs.num_hSigma(ifactor), mnhs.hRangeSigma(ifactor,:), mnhs.hScale)};
tregs{1,1}.opt              = mnhs.opt;
tregs{1,1}.regs             = [];
tregs{1,1}.opt.calc_eps     = true;
tregs{1,1}.opt.calc_sigma   = true;
tregs{1,1}.opt.calc_musigma = true;
tregs{1,1}.opt.calc_z       = true;

%% Level 2: batch layer — read batchModel config or use defaults
%  mnhs.batchModel.meanType    : 'constant' | 'functional' (default: 'functional')
%  mnhs.batchModel.varType     : 'constant' | 'functional' (default: 'functional')
%  mnhs.batchModel.meanFactors : cell of factor names      (default: mnhs.factor)
%  mnhs.batchModel.varFactors  : cell of factor names      (default: mnhs.factor)
%  Constraint: constant mean => constant variance (enforced here).
if isfield(mnhs, 'batchModel')
    bm = mnhs.batchModel;
    if ~isfield(bm, 'meanType'),    bm.meanType    = 'functional'; end
    if ~isfield(bm, 'varType'),     bm.varType     = 'functional'; end
    if ~isfield(bm, 'meanFactors'), bm.meanFactors = mnhs.factor;   end
    if ~isfield(bm, 'varFactors'),  bm.varFactors  = mnhs.factor;  end
    % Constraint: constant mean => constant variance
    if strcmp(bm.meanType, 'constant')
        if strcmp(bm.varType, 'functional')
            warning('GNM:batchModel', ...
                'Batch variance forced to constant because batch mean is constant.');
        end
        bm.varType = 'constant';
    end
    % Translate to factor lists
    if strcmp(bm.meanType, 'constant')
        batch_factor_mu = {};
    else
        batch_factor_mu = bm.meanFactors;
    end
    if strcmp(bm.varType, 'constant')
        batch_factor_sigma = {};
    else
        batch_factor_sigma = bm.varFactors;
    end
    % Store resolved config back
    mnhs.batchModel = bm;
else
    % Default: functional mean and variance by all factors (backward compatible)
    batch_factor_mu    = mnhs.factor;
    batch_factor_sigma = mnhs.factor;
end

tregs{2,1}.condition        = {mnhs.level{2}; unique(mnhs.T.(mnhs.level{2}))};
tregs{2,1}.resp             = append(mnhs.level{2}, mnhs.resp);
tregs{2,1}.factor_mu        = batch_factor_mu;
tregs{2,1}.factor_sigma     = batch_factor_sigma;
% Compute bandwidth grid; empty factors => constant model, no bandwidth needed
batch_factors_all = unique([batch_factor_mu, batch_factor_sigma], 'stable');
if ~isempty(batch_factors_all)
    ifactor = find_char(mnhs.factor, batch_factors_all);
    tregs{2,1}.hList = {getHlist(mnhs.num_hMu(ifactor), mnhs.hRangeMu(ifactor,:), mnhs.hScale); ...
        getHlist(mnhs.num_hSigma(ifactor), mnhs.hRangeSigma(ifactor,:), mnhs.hScale)};
else
    tregs{2,1}.hList = {[]; []};
end
tregs{2,1}.opt              = mnhs.opt;
tregs{2,1}.regs             = [];
tregs{2,1}.opt.calc_eps     = true;
tregs{2,1}.opt.calc_sigma   = true;
tregs{2,1}.opt.calc_musigma = true;
tregs{2,1}.opt.calc_z       = true;

%%
mnhs.tregs = tregs;
mnhs.opt.isremove = [false, true];

end
