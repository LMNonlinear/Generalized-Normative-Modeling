function T = nureg_table_predict(T, condition, factor_mu, factor_sigma, resp, opt, regs)
% nureg_table_predict  Predict mu/sigma/z for new data using trained griddedInterpolant.
%
%   T = nureg_table_predict(T, condition, factor_mu, factor_sigma, resp, opt, regs)
%
%   Uses the trained griddedInterpolant objects in regs(1).fpp_yhat (mu)
%   and regs(2).fpp_yhat (musigma) to predict on new factor coordinates.
%   Supports both real and complex responses (sep_complex).
%   Supports constant batch models (empty factor_mu/factor_sigma).
%   Processes one response at a time for memory safety.
%
%   See also: tnureg_zmap_predict, gnm_predict

%% Defaults
num_resp = length(resp);
if ~isfield(opt, 'resp_mu'),      opt.resp_mu = append('mu', resp);           end
if ~isfield(opt, 'resp_eps'),     opt.resp_eps = append('eps', resp);         end
if ~isfield(opt, 'resp_sigma'),   opt.resp_sigma = append('sigma', resp);     end
if ~isfield(opt, 'resp_musigma'), opt.resp_musigma = append('musigma', resp); end
if ~isfield(opt, 'resp_z'),       opt.resp_z = append('z', resp);             end

if ~isfield(opt, 'calc_eps'),      opt.calc_eps = false;      end
if ~isfield(opt, 'calc_sigma'),    opt.calc_sigma = false;    end
if ~isfield(opt, 'calc_musigma'),  opt.calc_musigma = false;  end
if ~isfield(opt, 'calc_z'),        opt.calc_z = false;        end

%% Extract subtable for this condition
cols_needed = resp(:);
if ~isempty(factor_mu),    cols_needed = [cols_needed; factor_mu(:)];    end
if ~isempty(factor_sigma), cols_needed = [cols_needed; factor_sigma(:)]; end
[subtable] = get_subtable(T, condition, cols_needed);

%% Prepare factor matrices
if ischar(factor_mu), factor_mu = {factor_mu}; end
if ~isempty(factor_sigma) && ischar(factor_sigma), factor_sigma = {factor_sigma}; end

n_sub = size(subtable, 1);
has_mu_factors = ~isempty(factor_mu);
has_sigma_factors = ~isempty(factor_sigma);

if has_mu_factors
    dx = length(factor_mu);
    x_mu = zeros(n_sub, dx);
    for ix = 1:dx
        x_mu(:, ix) = table2array(get_subtable(subtable, [], factor_mu(ix)));
    end
else
    x_mu = [];
end

if has_sigma_factors
    dx_s = length(factor_sigma);
    x_sigma = zeros(n_sub, dx_s);
    for ix = 1:dx_s
        x_sigma(:, ix) = table2array(get_subtable(subtable, [], factor_sigma(ix)));
    end
else
    x_sigma = [];
end

y = table2array(get_subtable(subtable, [], resp));
subtable = [];  % free memory

%% Precompute interpolation cell arrays (reused across responses)
if has_mu_factors
    xc_mu = mat2cell(x_mu, n_sub, ones(1, size(x_mu, 2)));
end
if has_sigma_factors
    xc_sig = mat2cell(x_sigma, n_sub, ones(1, size(x_sigma, 2)));
elseif has_mu_factors
    xc_sig = xc_mu;  % fallback: same factors for sigma
end

%% Process each response independently
for ir = 1:num_resp
    is_complex = (opt.sep_complex(ir) == 1);
    y_r = y(:, ir);

    %% Predict mu
    if has_mu_factors
        mu_r = interp_single_resp(regs(1).fpp_yhat, xc_mu, ir, false, []);
    else
        % Constant mean: grand mean per condition (matches cv_nureg_table.m:106)
        mu_r = mean(y_r, 1) .* ones(n_sub, 1);
    end

    %% Residual
    eps_r = y_r - mu_r;

    %% Variance of residual (for sigma output)
    sigma_r = [];
    if opt.calc_sigma
        sigma_r = real(eps_r).^2;
        if is_complex
            sigma_r = sigma_r + 1j * imag(eps_r).^2;
        end
    end

    %% Predict musigma and compute z
    z_r = [];
    musigma_r = [];
    if opt.calc_musigma || opt.calc_z
        if has_sigma_factors
            musigma_r = interp_single_resp(regs(2).fpp_yhat, xc_sig, ir, true, regs(2).d);
        else
            % Constant variance: grand variance per condition (matches cv_nureg_table.m:144)
            ysigma_tmp = real(eps_r).^2;
            if is_complex
                ysigma_tmp = ysigma_tmp + 1j * imag(eps_r).^2;
            end
            musigma_r = mean(ysigma_tmp, 1) .* ones(n_sub, 1);
        end
        musigma_r_real = max(real(musigma_r), eps);

        if opt.calc_z
            z_r = real(eps_r) ./ sqrt(musigma_r_real);
            % Complex part (matches cv_nureg_table.m:155-156, HarMNqEEG reference)
            if is_complex
                musigma_r_imag = max(imag(musigma_r), eps);
                z_r = z_r + 1j * imag(eps_r) ./ sqrt(musigma_r_imag);
            end
        end
    end

    %% Save to table
    T = asnarray2table(T, condition, opt.resp_mu(ir), mu_r);

    if opt.calc_eps
        T = asnarray2table(T, condition, opt.resp_eps(ir), eps_r);
    end

    if opt.calc_sigma
        T = asnarray2table(T, condition, opt.resp_sigma(ir), sigma_r);
    end

    if opt.calc_musigma && (~isfield(opt, 'save_musigma') || opt.save_musigma)
        T = asnarray2table(T, condition, opt.resp_musigma(ir), musigma_r);
    end

    if opt.calc_z
        T = asnarray2table(T, condition, opt.resp_z(ir), z_r);
    end
end

end

%% ========================================================================
%  Local: interpolate a single response from a griddedInterpolant
%  ========================================================================
function val = interp_single_resp(fpp_full, xc, ir, is_variance, d)
% Slice fpp to a single response index and interpolate.
%   fpp_full   - original griddedInterpolant with all responses
%   xc         - cell array of factor vectors {x1, x2, ...}
%   ir         - response index to extract
%   is_variance - if true, apply log/exp variance transform
%   d          - variance scaling vector (only used when is_variance=true)

dx = numel(xc);
ndcolon = repmat({':'}, 1, dx);

% Slice Values to single response (avoid full copy — only copy the slice)
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
