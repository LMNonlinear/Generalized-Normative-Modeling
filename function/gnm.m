function mnhs = gnm(path_csv, formulaStr, varargin)
% gnm  Create GNM configuration from formula syntax.
%
%   mnhs = gnm(path_csv, formula)
%   mnhs = gnm(path_csv, formula, 'Name', Value, ...)
%
%   This is a convenience front-end for the struct-based interface.
%   It parses a formula string and builds the mnhs configuration struct
%   that can be passed to gnm_fit().
%
%   Formula syntax:
%     'y1+y2 ~ s(age) | site'               — basic usage
%     'y1 ~ s(age,icv) | site[constant]'     — constant batch model
%     'y1 ~ s(time) | scanner[functional(time)]' — functional with explicit factors
%     'y1 ~ s(age)'                          — no batch (uses gnm_initialize_model default)
%
%   Name-Value parameters:
%     'inform'      - cell array of metadata column names (default: {'name'})
%     'hRangeMu'    - bandwidth range for mean  (nfactors x 2 matrix)
%     'num_hMu'     - number of bandwidths for mean (1 x nfactors vector)
%     'hRangeSigma' - bandwidth range for variance (nfactors x 2 matrix)
%     'num_hSigma'  - number of bandwidths for variance (1 x nfactors vector)
%     'hScale'      - bandwidth scaling function (default: @logspace)
%     'tag'         - project identifier tag
%     'link'        - response transform: 'identity'|'log'|'logit'|'sqrt'|'probit'
%                     or a struct from gnm_link() (default: 'identity')
%
%   Examples:
%     % Minimal
%     mnhs = gnm('data.csv', 'y1+y2 ~ s(age) | site');
%     mnhs = gnm_fit(mnhs);
%
%     % With options
%     mnhs = gnm('data.csv', 'y1 ~ s(age,icv) | site[constant]', ...
%                'inform', {'name','sex'}, 'hRangeMu', [0.3,0.5; 0.3,0.5]);
%     mnhs = gnm_fit(mnhs);
%
%   See also: gnm_fit, gnm_parse_formula, gnm_build_default_pipeline

%% Validate inputs
if ~ischar(path_csv) && ~isstring(path_csv)
    error('GNM:gnm', 'path_csv must be a string.');
end
path_csv = char(path_csv);

if ~ischar(formulaStr) && ~isstring(formulaStr)
    error('GNM:gnm', 'formulaStr must be a string.');
end
formulaStr = char(formulaStr);

if ~isfile(path_csv)
    error('GNM:gnm', 'CSV file not found: %s', path_csv);
end

%% Parse name-value arguments
p = inputParser;
p.addParameter('inform',      {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('hRangeMu',    [], @isnumeric);
p.addParameter('num_hMu',     [], @isnumeric);
p.addParameter('hRangeSigma', [], @isnumeric);
p.addParameter('num_hSigma',  [], @isnumeric);
p.addParameter('hScale',      [], @(x) isa(x, 'function_handle'));
p.addParameter('tag',         '', @(x) ischar(x) || isstring(x));
p.addParameter('link',  'identity', @(x) ischar(x) || isstring(x) || isstruct(x));
p.parse(varargin{:});
opts = p.Results;

%% Parse formula
parsed = gnm_parse_formula(formulaStr);

%% Build mnhs struct
mnhs = struct();
mnhs.path_csv = path_csv;
mnhs.formula  = formulaStr;  % store for debugging / display

% Core fields from parser
mnhs.resp   = parsed.resp;
mnhs.factor = parsed.factor;

% Batch: only set if specified in formula (let gnm_initialize_model use default otherwise)
if ~isempty(parsed.batch)
    mnhs.batch = parsed.batch;
end

% Batch model: build from parser output
if ~isempty(fieldnames(parsed.batchModel))
    bm = parsed.batchModel;
    % For 'functional' without explicit factors, fill in from mnhs.factor
    if strcmp(bm.meanType, 'functional') && ~isfield(bm, 'meanFactors')
        bm.meanFactors = parsed.factor;
    end
    if strcmp(bm.varType, 'functional') && ~isfield(bm, 'varFactors')
        bm.varFactors = parsed.factor;
    end
    mnhs.batchModel = bm;
end

%% Apply name-value overrides
if ~isempty(opts.inform)
    if ischar(opts.inform), opts.inform = {opts.inform}; end
    if isstring(opts.inform), opts.inform = cellstr(opts.inform); end
    mnhs.inform = opts.inform;
end

if ~isempty(opts.hRangeMu),    mnhs.hRangeMu    = opts.hRangeMu;    end
if ~isempty(opts.num_hMu),     mnhs.num_hMu     = opts.num_hMu;     end
if ~isempty(opts.hRangeSigma), mnhs.hRangeSigma = opts.hRangeSigma; end
if ~isempty(opts.num_hSigma),  mnhs.num_hSigma  = opts.num_hSigma;  end
if ~isempty(opts.hScale),      mnhs.hScale      = opts.hScale;      end

if ~isempty(opts.tag)
    mnhs.tag = char(opts.tag);
end

% Link function
if isstruct(opts.link)
    mnhs.link = opts.link;
else
    mnhs.link = gnm_link(char(opts.link));
end

%% Print configuration summary
fprintf('\n=== GNM Configuration ===\n');
fprintf('  Formula : %s\n', formulaStr);
fprintf('  CSV     : %s\n', path_csv);
fprintf('  resp    : {%s}\n', strjoin(mnhs.resp, ', '));
fprintf('  factor  : {%s}\n', strjoin(mnhs.factor, ', '));
if isfield(mnhs, 'batch')
    fprintf('  batch   : {%s}\n', strjoin(mnhs.batch, ', '));
else
    fprintf('  batch   : (default)\n');
end
if isfield(mnhs, 'batchModel')
    bm = mnhs.batchModel;
    fprintf('  batchModel.meanType = %s\n', bm.meanType);
    fprintf('  batchModel.varType  = %s\n', bm.varType);
    if isfield(bm, 'meanFactors')
        fprintf('  batchModel.meanFactors = {%s}\n', strjoin(bm.meanFactors, ', '));
    end
    if isfield(bm, 'varFactors')
        fprintf('  batchModel.varFactors  = {%s}\n', strjoin(bm.varFactors, ', '));
    end
end
if isfield(mnhs, 'inform')
    fprintf('  inform  : {%s}\n', strjoin(mnhs.inform, ', '));
end
if isfield(mnhs, 'link')
    fprintf('  link    : %s\n', mnhs.link.name);
end
fprintf('=========================\n\n');

end
