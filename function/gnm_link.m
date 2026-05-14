function link = gnm_link(name)
% gnm_link  Create a link function struct for response transformation.
%
%   link = gnm_link('log')
%   link = gnm_link('identity')
%
%   Transforms non-Gaussian responses to approximately Normal scale.
%   The pipeline runs on g(y), producing z-scores that are quantile
%   residuals of the implied distribution.
%
%   Supported links:
%     'identity' - g(y)=y,         default (Normal)
%     'log'      - g(y)=log(y),    y>0  (LogNormal: brain volume, power spectra)
%     'logit'    - g(y)=logit(y),  0<y<1 (Logistic-Normal: FA, proportions)
%     'sqrt'     - g(y)=sqrt(y),   y>=0  (approx Poisson: counts)
%     'probit'   - g(y)=probit(y), 0<y<1 (Probit-Normal: alternative to logit)
%
%   Output struct fields:
%     link.name         - string name
%     link.forward      - @(y) forward transform g(y)
%     link.inverse      - @(eta) inverse transform g^{-1}(eta)
%     link.domain_check - @(y) logical, true if y is in valid domain
%     link.derivative   - @(y) g'(y) (for delta-method SE back-transform)
%
%   Note: probit uses erfinv/erf (no Statistics Toolbox dependency).
%
%   See also: gnm, gnm_fit

arguments
    name {mustBeTextScalar} = 'identity'
end
name = lower(char(name));

switch name
    case 'identity'
        link.name         = 'identity';
        link.forward      = @(y) y;
        link.inverse      = @(eta) eta;
        link.domain_check = @(y) true(size(y));
        link.derivative   = @(y) ones(size(y));

    case 'log'
        link.name         = 'log';
        link.forward      = @(y) log(y);
        link.inverse      = @(eta) exp(eta);
        link.domain_check = @(y) y > 0;
        link.derivative   = @(y) 1 ./ y;

    case 'logit'
        link.name         = 'logit';
        link.forward      = @(y) log(y ./ (1 - y));
        link.inverse      = @(eta) 1 ./ (1 + exp(-eta));
        link.domain_check = @(y) y > 0 & y < 1;
        link.derivative   = @(y) 1 ./ (y .* (1 - y));

    case 'sqrt'
        link.name         = 'sqrt';
        link.forward      = @(y) sqrt(y);
        link.inverse      = @(eta) eta .^ 2;
        link.domain_check = @(y) y >= 0;
        link.derivative   = @(y) 0.5 ./ sqrt(y);

    case 'probit'
        % probit = Phi^{-1}(y), implemented via erfinv to avoid
        % Statistics Toolbox dependency (norminv/normcdf).
        link.name         = 'probit';
        link.forward      = @(y) sqrt(2) * erfinv(2*y - 1);
        link.inverse      = @(eta) 0.5 * (1 + erf(eta / sqrt(2)));
        link.domain_check = @(y) y > 0 & y < 1;
        link.derivative   = @(y) sqrt(2*pi) * exp(erfinv(2*y - 1).^2);

    otherwise
        error('GNM:gnm_link', ...
            'Unknown link "%s". Supported: identity, log, logit, sqrt, probit.', name);
end

end
