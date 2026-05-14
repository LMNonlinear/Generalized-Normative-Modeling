function T = tnureg_zmap_predict(T, y, tregs, opt)
% tnureg_zmap_predict  Predict z-scores and ystar for new data using trained models.
%
%   T = tnureg_zmap_predict(T, y, tregs, opt)
%
%   Multi-level prediction engine: iterates through trained regression
%   levels (global, batch), predicts mu/sigma at each, computes z-scores,
%   and reconstructs ystar by reversing the batch effect.
%   Supports complex responses and constant batch models.
%
%   See also: gnm_predict, nureg_table_predict

if isstruct(tregs)
    tregs = {tregs};
end
opt = set_defaults(opt, 'standardize', false);
opt = set_defaults(opt, 'unstandardize', false);
opt = set_defaults(opt, 'maxibf', 1);

dy = numel(y);
n = size(T, 1);

%% Standardize (if training used it)
if opt.standardize
    func = @(x) [mean(x, 1); std(real(x), 1) + 1i * std(imag(x))];
    mu = func(reshape(table2array(get_subtable(T, [], y)), n, [], dy));
    stdy = mu(2, :, :);
    mu = mu(1, :, :);
    for i = 1:dy
        T.(y{i}) = T.(y{i}) - mu(:, :, i);
        T.(y{i}) = real(T.(y{i})) ./ real(stdy(:, :, i)) + ...
            1i * imag(T.(y{i})) ./ imag(stdy(:, :, i));
    end
end

%% Predict through levels
num_level = size(tregs, 1);

for ibf = 1:opt.maxibf
    for i = 1:num_level
        fprintf('  [predict] level %d/%d\n', i, num_level);

        % Clone response columns for this level's input
        if i == 1 && ibf == 1
            T = clone_table_var(T, tregs{i, 1}.resp, y);
        elseif i == 1 && ibf > 1
            if tregs{num_level, 1}.opt.calc_z
                T = clone_table_var(T, tregs{i, 1}.resp, tregs{num_level, 1}.opt.resp_z);
            else
                T = clone_table_var(T, tregs{i, 1}.resp, tregs{num_level, 1}.opt.resp_eps);
            end
        else
            if tregs{i - 1, 1}.opt.calc_z
                T = clone_table_var(T, tregs{i, 1}.resp, tregs{i - 1, 1}.opt.resp_z);
            else
                T = clone_table_var(T, tregs{i, 1}.resp, tregs{i - 1, 1}.opt.resp_eps);
            end
        end

        % Predict for each batch condition using trained regs
        T = tnureg_predict_level(T, tregs{i, 1});
    end
end

%% Compact intermediate columns before ystar reconstruction (frees memory)
if opt.compact.y
    T = asnarray2table(T, [], y, []);
end
for i = 1:num_level
    if opt.compact.musigma
        ms_cols = tregs{i, 1}.opt.resp_musigma;
        if ~isempty(ms_cols) && ismember(ms_cols{1}, T.Properties.VariableNames)
            T = asnarray2table(T, [], ms_cols, []);
        end
    end
    if opt.compact.mu
        mu_cols = tregs{i, 1}.opt.resp_mu;
        if ~isempty(mu_cols) && ismember(mu_cols{1}, T.Properties.VariableNames)
            T = asnarray2table(T, [], mu_cols, []);
        end
    end
    if opt.compact.eps
        eps_cols = tregs{i, 1}.opt.resp_eps;
        if ~isempty(eps_cols) && ismember(eps_cols{1}, T.Properties.VariableNames)
            T = asnarray2table(T, [], eps_cols, []);
        end
    end
end

%% Reconstruct ystar per-response (memory-efficient: recompute mu/musigma from fpp)
ystar = append('ystar', y);

if tregs{num_level}.opt.calc_musigma
    T = clone_table_var(T, ystar, tregs{num_level, 1}.opt.resp_z);
else
    T = clone_table_var(T, ystar, tregs{num_level, 1}.opt.resp_eps);
end

% Build factor vectors once for each level (needed for fpp interpolation)
level_x = cell(num_level, 2);  % {x_mu, x_sigma} per level
for i = 1:num_level
    treg = tregs{i, 1};
    if ~isempty(treg.factor_mu)
        x_mu_i = zeros(n, numel(treg.factor_mu));
        for ix = 1:numel(treg.factor_mu)
            x_mu_i(:, ix) = T.(treg.factor_mu{ix});
        end
        level_x{i, 1} = x_mu_i;
    end
    if ~isempty(treg.factor_sigma)
        x_sig_i = zeros(n, numel(treg.factor_sigma));
        for ix = 1:numel(treg.factor_sigma)
            x_sig_i(:, ix) = T.(treg.factor_sigma{ix});
        end
        level_x{i, 2} = x_sig_i;
    end
end

for i = num_level:-1:1
    if ~opt.isremove(i)
        treg = tregs{i, 1};
        cond_col = treg.condition{1};
        cond_vals = treg.condition{2};
        row_cond = T.(cond_col);
        has_mu_factors = ~isempty(treg.factor_mu);
        has_sigma_factors = ~isempty(treg.factor_sigma);

        for ir = 1:dy
            ys_val = T.(ystar{ir});
            mu_accum = zeros(n, 1);
            ms_accum = ones(n, 1);  % default to 1 (identity for multiplication)
            is_complex = (opt.sep_complex(ir) == 1);

            for ic = 1:numel(cond_vals)
                mask = strcmp(row_cond, cond_vals{ic});
                if ~any(mask), continue; end
                regs_ic = treg.regs{ic};
                n_mask = sum(mask);

                % Predict mu for single response ir
                if has_mu_factors
                    x_sub = level_x{i, 1}(mask, :);
                    xc = mat2cell(x_sub, n_mask, ones(1, size(x_sub, 2)));
                    mu_single = interp_single_resp(regs_ic(1).fpp_yhat, xc, ir, false, []);
                    mu_accum(mask) = mu_single;
                else
                    % Constant mean: not available from fpp, use 0 (no shift)
                    % For constant models, the mean was the grand mean of the condition.
                    % It's stored in regs(1).yhat during training but compacted away.
                    % For ystar reconstruction, this level's contribution is already
                    % embedded in the z-score computation.
                    mu_accum(mask) = 0;
                end

                % Predict musigma for single response ir
                if treg.opt.calc_musigma
                    if has_sigma_factors
                        x_sub_s = level_x{i, 2}(mask, :);
                        xc_s = mat2cell(x_sub_s, n_mask, ones(1, size(x_sub_s, 2)));
                        ms_single = interp_single_resp(regs_ic(2).fpp_yhat, xc_s, ir, true, regs_ic(2).d);
                        ms_single = max(real(ms_single), eps);
                        ms_accum(mask) = ms_single;
                    else
                        % Constant variance: use 1 (identity scaling)
                        ms_accum(mask) = 1;
                    end
                end
            end

            % ystar = mu + musigma * ystar_inner (real + complex parts)
            if treg.opt.calc_musigma
                T.(ystar{ir}) = (real(mu_accum) + real(ms_accum) .* real(ys_val));
                if is_complex
                    T.(ystar{ir}) = T.(ystar{ir}) + ...
                        1i * (imag(mu_accum) + imag(ms_accum) .* imag(ys_val));
                end
            else
                T.(ystar{ir}) = real(mu_accum) + real(ys_val);
                if is_complex
                    T.(ystar{ir}) = T.(ystar{ir}) + ...
                        1i * (imag(mu_accum) + imag(ys_val));
                end
            end
        end
    end
end

%% Unstandardize ystar
if opt.unstandardize
    for i = 1:dy
        T.(ystar{i}) = real(T.(ystar{i})) .* real(stdy(:, :, i)) + ...
            1i * imag(T.(ystar{i})) .* imag(stdy(:, :, i));
        T.(ystar{i}) = T.(ystar{i}) + mu(:, :, i);
    end
end

%% Compact remaining
if opt.compact.sigma
    for i = 1:num_level
        sig_cols = tregs{i, 1}.opt.resp_sigma;
        if ~isempty(sig_cols) && ismember(sig_cols{1}, T.Properties.VariableNames)
            T = asnarray2table(T, [], sig_cols, []);
        end
    end
end
if opt.compact.z
    for i = 1:num_level
        T = asnarray2table(T, [], tregs{i, 1}.opt.resp_z, []);
    end
end
if opt.compact.ystar
    T = asnarray2table(T, [], ystar, []);
end

end

%% ========================================================================
%  Local: interpolate single response from griddedInterpolant
%  ========================================================================
function val = interp_single_resp(fpp_full, xc, ir, is_variance, d)
dx = numel(xc);
ndcolon = repmat({':'}, 1, dx);
vals_slice = fpp_full.Values(ndcolon{:}, ir);

if is_variance
    d_ir = d(ir);
    vals_slice = log(vals_slice * d_ir);
end

fpp_single = griddedInterpolant(fpp_full.GridVectors, vals_slice, ...
    fpp_full.Method, fpp_full.ExtrapolationMethod);
val = fpp_single(xc{:});

if is_variance
    val = exp(val) / d_ir;
end

if ~iscolumn(val), val = val(:); end
end

%% ========================================================================
%  Local: predict for one level across all conditions
%  ========================================================================
function T = tnureg_predict_level(T, tregs)
trained_batches = unique(tregs.condition{2});
new_batches = unique(T.(tregs.condition{1}));
idx = find_char(trained_batches, new_batches);

% Robust check: find_char returns 0 or empty for unmatched batches (#4)
if isnumeric(idx)
    bad = (idx == 0);
else
    bad = cellfun(@isempty, idx);
end
if any(bad)
    if iscell(new_batches)
        missing_str = strjoin(new_batches(bad), ', ');
    else
        missing_str = char(new_batches(bad));
    end
    error('GNM:predict:unknownBatch', ...
        'Batch(es) [%s] not found in trained model.\nAvailable: {%s}\nUse reRefBatch to map new batches to existing ones.', ...
        missing_str, strjoin(trained_batches, ', '));
end

num_cond = length(idx);
for icond = 1:num_cond
    cond = {tregs.condition{1}, tregs.condition{2}{idx(icond)}};
    fprintf('    predict [%s] in [%s]\n', cond{2}, cond{1});
    T = nureg_table_predict(T, ...
        cond, ...
        tregs.factor_mu, ...
        tregs.factor_sigma, ...
        tregs.resp, ...
        tregs.opt, ...
        tregs.regs{idx(icond)});
end

end
