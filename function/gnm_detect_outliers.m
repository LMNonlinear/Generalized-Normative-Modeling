function mnhs = gnm_detect_outliers(mnhs)
% gnm_detect_outliers  Detect and flag individual outliers using robust PCA.
%
%   mnhs = gnm_detect_outliers(mnhs)
%
%   When mnhs.clean_individual_outlier is true, performs:
%     1. Quick global smooth to get residuals
%     2. Robust PCA on z-scores to detect multivariate outliers
%     3. Saves cleaned data and exits
%
%   When false (default), this is a no-op pass-through.
%
%   See also: gnm_compute_zscores, pcaRobustCov

if mnhs.clean_individual_outlier
    %% 0. initialize temp smooth
    opt = mnhs.opt;
    opt.y_keep_all   = false;
    opt.compact.y    = false;
    opt.compact.eps  = false;
    opt.compact.z    = false;
    opt.calc_musigma = true;
    opt.calc_z       = true;

    tregs = mnhs.tregs{1};
    tregs.opt.y_keep_all  = false;
    tregs.opt.compact.y   = false;
    tregs.opt.compact.eps = false;
    tregs.opt.compact.z   = false;
    tregs.opt.calc_musigma = true;
    tregs.opt.calc_z       = true;
    tregs.factor           = mnhs.factor;

    ifactor     = 1:length(mnhs.factor);
    hScale      = mnhs.hScale;
    hRangeMu    = mnhs.hRangeMu;
    num_hMu     = mnhs.num_hMu;
    hRangeSigma = mnhs.hRangeSigma;
    num_hSigma  = mnhs.num_hSigma;

    tregs.hList = {getHlist(num_hMu(ifactor), hRangeMu(ifactor,:), hScale); ...
        getHlist(num_hSigma(ifactor), hRangeSigma(ifactor,:), hScale)};

    %% 1. global smooth
    [mnhs.T, tregs] = cv_tnureg_zmap(mnhs.T, mnhs.resp, tregs, opt);
    tregs = tregs{1};

    %% 2. PCA Robust Cov
    if tregs.opt.calc_z
        resp = tregs.opt.resp_z;
    else
        resp = tregs.opt.resp_eps;
    end
    mnhs.T = pcaRobustCov(mnhs.T, resp, 3, mnhs.batch(1), mnhs.tag, mnhs.clean_individual_outlier, [], mnhs.inform{1});

    %% compact to save memory
    mnhs.T = asnarray2table(mnhs.T, [], tregs.resp, []);
    mnhs.T = asnarray2table(mnhs.T, [], tregs.opt.resp_mu, []);
    mnhs.T = asnarray2table(mnhs.T, [], tregs.opt.resp_eps, []);
    mnhs.T = asnarray2table(mnhs.T, [], tregs.opt.resp_sigma, []);
    mnhs.T = asnarray2table(mnhs.T, [], tregs.opt.resp_musigma, []);
    mnhs.T = asnarray2table(mnhs.T, [], tregs.opt.resp_z, []);
    mnhs.T = asnarray2table(mnhs.T, [], tregs.opt.ystar, []);
    filename = strrep(mnhs.path_csv, '.csv', '_cleanOutlier.csv');
    writetable(mnhs.T, filename)
    diary off
    error('exit after clean Outlier and save')
end

end
