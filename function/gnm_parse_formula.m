function parsed = gnm_parse_formula(formulaStr)
% gnm_parse_formula  Parse GNM formula string into struct fields.
%
%   parsed = gnm_parse_formula('y1+y2 ~ s(age,icv) | site[constant]')
%
%   Formula syntax:
%     LHS ~ RHS
%     LHS = resp1 + resp2 + ...
%     RHS = smooth_terms [ | batch_spec ]
%     smooth_terms = s(var1,var2,...) [+ s(var3,...)]
%     batch_spec   = batchvar [ [model] ]
%     model        = constant | functional | functional(var1,var2)
%
%   Output struct:
%     parsed.resp       - cell array of response variable names
%     parsed.factor     - cell array of factor variable names (unique, stable)
%     parsed.batch      - cell array with batch variable name (empty if omitted)
%     parsed.batchModel - struct with meanType, varType, meanFactors, varFactors
%                         (empty struct if no batch or no model spec)

arguments
    formulaStr (1,1) string
end

formulaStr = char(strtrim(formulaStr));
if isempty(formulaStr)
    error('GNM:formula', 'Formula string cannot be empty.');
end

%% Split on '~'
parts_tilde = strsplit(formulaStr, '~');
if numel(parts_tilde) ~= 2
    error('GNM:formula', ...
        'Formula must contain exactly one "~". Got: %s', formulaStr);
end
lhs = strtrim(parts_tilde{1});
rhs = strtrim(parts_tilde{2});

%% Parse LHS: resp1 + resp2 + ...
resp_parts = strsplit(lhs, '+');
resp = cellfun(@strtrim, resp_parts, 'UniformOutput', false);
resp = resp(~cellfun(@isempty, resp));
if isempty(resp)
    error('GNM:formula', 'No response variables found in LHS: %s', lhs);
end
for i = 1:numel(resp)
    if ~isvarname(resp{i})
        error('GNM:formula', ...
            'Invalid response variable name: "%s". Must be a valid MATLAB identifier.', resp{i});
    end
end

%% Split RHS on '|' (batch separator)
idx_pipe = strfind(rhs, '|');
if isempty(idx_pipe)
    smooth_part = rhs;
    batch_part  = '';
elseif numel(idx_pipe) == 1
    smooth_part = strtrim(rhs(1:idx_pipe-1));
    batch_part  = strtrim(rhs(idx_pipe+1:end));
else
    error('GNM:formula', ...
        'Formula RHS must contain at most one "|". Got: %s', rhs);
end

%% Parse smooth terms: s(var1,var2,...) + s(var3,...)
tokens = regexp(smooth_part, 's\(\s*([^)]+)\s*\)', 'tokens');
if isempty(tokens)
    error('GNM:formula', ...
        'No s() terms found in RHS: "%s". Use e.g. s(age) or s(age,icv).', smooth_part);
end
factor = {};
for i = 1:numel(tokens)
    vars = strsplit(tokens{i}{1}, ',');
    vars = cellfun(@strtrim, vars, 'UniformOutput', false);
    vars = vars(~cellfun(@isempty, vars));
    for j = 1:numel(vars)
        if ~isvarname(vars{j})
            error('GNM:formula', ...
                'Invalid factor variable name: "%s". Must be a valid MATLAB identifier.', vars{j});
        end
    end
    factor = [factor, vars]; %#ok<AGROW>
end
factor = unique(factor, 'stable');
if isempty(factor)
    error('GNM:formula', 'No factor variables extracted from s() terms.');
end

%% Parse batch part: batchvar[model]
batch = {};
batchModel = struct();

if ~isempty(batch_part)
    batch_spec = strtrim(batch_part);
    open_idx = strfind(batch_spec, '[');
    close_idx = strfind(batch_spec, ']');

    if isempty(open_idx) && isempty(close_idx)
        batch_name = strtrim(batch_spec);
        model_spec = '';
    else
        if numel(open_idx) ~= 1 || numel(close_idx) ~= 1 || open_idx > close_idx
            error('GNM:formula', ...
                'Invalid batch specification: "%s". Expected e.g. site or site[constant].', batch_part);
        end

        batch_name = strtrim(batch_spec(1:open_idx-1));
        model_spec = strtrim(batch_spec(open_idx+1:close_idx-1));
        trailing = strtrim(batch_spec(close_idx+1:end));

        if isempty(batch_name) || ~isempty(trailing) || contains(model_spec, '[') || contains(model_spec, ']')
            error('GNM:formula', ...
                'Invalid batch specification: "%s". Expected e.g. site or site[constant].', batch_part);
        end
    end

    if ~isvarname(batch_name)
        error('GNM:formula', ...
            'Invalid batch variable name: "%s". Must be a valid MATLAB identifier.', batch_name);
    end
    batch = {batch_name};

    if ~isempty(model_spec)
        batchModel = parse_model_spec(model_spec);
    end
end

%% Build output
parsed.resp       = resp;
parsed.factor     = factor;
parsed.batch      = batch;
parsed.batchModel = batchModel;

end


function bm = parse_model_spec(spec)
% parse_model_spec  Parse batch model specification string.
%
%   'constant'            → meanType=constant, varType=constant
%   'functional'          → meanType=functional, varType=functional
%   'functional(age)'     → meanType=functional, varType=functional,
%                           meanFactors={'age'}, varFactors={'age'}
%   'functional(age,icv)' → same with multiple factors

    % try: functional(var1,var2,...)
    tok = regexp(spec, '^functional\(\s*([^)]+)\s*\)$', 'tokens');
    if ~isempty(tok)
        vars = strsplit(tok{1}{1}, ',');
        vars = cellfun(@strtrim, vars, 'UniformOutput', false);
        vars = vars(~cellfun(@isempty, vars));
        for i = 1:numel(vars)
            if ~isvarname(vars{i})
                error('GNM:formula', ...
                    'Invalid variable in batch model spec: "%s".', vars{i});
            end
        end
        bm.meanType    = 'functional';
        bm.varType     = 'functional';
        bm.meanFactors = vars;
        bm.varFactors  = vars;
        return;
    end

    % try: constant or functional (plain keywords)
    spec_lower = lower(strtrim(spec));
    switch spec_lower
        case 'constant'
            bm.meanType = 'constant';
            bm.varType  = 'constant';
        case 'functional'
            bm.meanType = 'functional';
            bm.varType  = 'functional';
        otherwise
            error('GNM:formula', ...
                'Unknown batch model spec: "%s". Use constant, functional, or functional(var1,...).', spec);
    end
end
