function mnhs = gnm_compute_zscores(mnhs)
% gnm_compute_zscores  Main GNM pipeline: outlier detection + harmonization.
%
%   mnhs = gnm_compute_zscores(mnhs)
%
%   Orchestrates the core GNM computation:
%     1. Robust outlier detection (gnm_detect_outliers)
%     2. Batch harmonization and z-score computation (gnm_harmonize_batch)
%
%   This function is set as mnhs.mnhfun by gnm_build_default_pipeline
%   and dispatched by gnm_run_pipeline.
%
%   See also: gnm_detect_outliers, gnm_harmonize_batch, gnm_run_pipeline

%% 1. robust outlier detection
disp('[gnm_detect_outliers]')
mnhs = gnm_detect_outliers(mnhs);

%% 2. harmonization + z-score computation
disp('[gnm_harmonize_batch]')
mnhs = gnm_harmonize_batch(mnhs);

end
