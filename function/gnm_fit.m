function mnhs = gnm_fit(mnhs)
% gnm_fit  Fit a GNM model from a configured mnhs struct.
%
%   mnhs = gnm_fit(mnhs)
%
%   Pipeline:
%     1. gnm_build_default_pipeline(mnhs) — build hierarchical model
%     2. Forward link transform            — g(y) on response columns
%     3. gnm_run_pipeline(mnhs)           — execute pipeline on g(y)
%     4. Back-transform                    — g^{-1} on ystar and global mu
%
%   The link function (if set) transforms the response to approximate
%   Normal scale before running the kernel regression pipeline. Z-scores
%   produced are quantile residuals of the implied distribution.
%
%   Use gnm() to construct mnhs from a formula string, then pass to
%   gnm_fit(). Alternatively, build mnhs manually (struct interface).
%
%   Example:
%     mnhs = gnm('data.csv', 'vol ~ s(age) | site', 'link', 'log');
%     mnhs = gnm_fit(mnhs);
%
%   See also: gnm, gnm_link, gnm_build_default_pipeline, gnm_run_pipeline

arguments
    mnhs (1,1) struct
end

% Validate required fields
assert(isfield(mnhs, 'path_csv') && ~isempty(mnhs.path_csv), ...
    'GNM:gnm_fit', 'mnhs.path_csv is required.');
assert(isfield(mnhs, 'resp') && ~isempty(mnhs.resp), ...
    'GNM:gnm_fit', 'mnhs.resp is required.');
assert(isfield(mnhs, 'factor') && ~isempty(mnhs.factor), ...
    'GNM:gnm_fit', 'mnhs.factor is required.');

%% Step 1: Build hierarchical model (reads CSV, creates tregs)
mnhs = gnm_build_default_pipeline(mnhs);

%% Step 2: Forward link transform — g(y) on response columns
has_link = isfield(mnhs, 'link') && ~strcmp(mnhs.link.name, 'identity');
if has_link
    link = mnhs.link;

    % Link functions are undefined for complex responses
    if any(mnhs.opt.sep_complex == 1)
        error('GNM:gnm_fit', ...
            'Link function "%s" is not supported for complex responses.', link.name);
    end

    for i = 1:numel(mnhs.resp)
        y_col = mnhs.T.(mnhs.resp{i});
        % Domain check (NaN values are transparent)
        valid = link.domain_check(y_col);
        valid(isnan(y_col)) = true;
        assert(all(valid, 'all'), ...
            'GNM:gnm_fit:domainViolation', ...
            'Link "%s" domain violation in column "%s".', link.name, mnhs.resp{i});
        % Apply forward transform
        mnhs.T.(mnhs.resp{i}) = link.forward(y_col);
    end
    fprintf('  [link] Applied forward transform: %s\n', link.name);
end

%% Step 3: Run pipeline on (transformed) data
mnhs = gnm_run_pipeline(mnhs);

%% Step 4: Back-transform — g^{-1} on ystar and global mu
if has_link
    link = mnhs.link;

    % ystar columns → original scale
    for i = 1:numel(mnhs.resp_ystar)
        col = mnhs.resp_ystar{i};
        if ismember(col, mnhs.T.Properties.VariableNames)
            mnhs.T.(col) = link.inverse(mnhs.T.(col));
        end
    end

    % global mu columns → original scale (for normative centile plots)
    if ~mnhs.opt.compact.mu
        mu_cols = mnhs.tregs{1}.opt.resp_mu;
        for i = 1:numel(mu_cols)
            col = mu_cols{i};
            if ismember(col, mnhs.T.Properties.VariableNames)
                mnhs.T.(col) = link.inverse(mnhs.T.(col));
            end
        end
    end

    % NOT back-transformed (by design):
    %   z-scores   — these ARE the quantile residuals, should be ~N(0,1)
    %   batch mu   — additive shift in transform scale, meaningless alone
    %   sigma, musigma — variance params in transform scale
    %   eps        — residual in transform scale
    fprintf('  [link] Applied inverse transform: %s (ystar, global mu)\n', link.name);
end

end
