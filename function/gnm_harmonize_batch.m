function mnhs = gnm_harmonize_batch(mnhs)
% gnm_harmonize_batch  Harmonize batch effects and compute z-scores.
%
%   mnhs = gnm_harmonize_batch(mnhs)
%
%   Performs two steps:
%     1. Hierarchical kernel regression to estimate normative trajectories
%        and batch effects (via cv_tnureg_zmap)
%     2. Model selection criteria computation (via tnureg_criteria)
%
%   See also: gnm_compute_zscores, cv_tnureg_zmap, tnureg_criteria

%% 1. harmonization
opt = mnhs.opt;
opt.compact.eps = false;
opt.compact.musigma = false;
opt.compact.sigma = false;
[mnhs.T, mnhs.tregs] = cv_tnureg_zmap(mnhs.T, mnhs.resp, mnhs.tregs, opt);

%% 2. get criteria for model selection
opt = mnhs.opt;
[mnhs.T, mnhs.crit] = tnureg_criteria(mnhs.T, mnhs.resp, mnhs.tregs, opt);

end
