function tregs = parseTregsFormula(formulaStr, config)
% PARSE FORMULA WITH MODEL CONFIGURATION
% Inputs:
%   formulaStr: Model formula (e.g., "y ~ f(freq,age) + f(sex,age)")
%   config: Configuration structure from gnm_initialize_model
% Output:
%   tregs: Structured array of term relationships

% Extract formula components
[responseVar, terms] = splitFormulaComponents(formulaStr);

% Initialize tregs structure with smart defaults
tregs = initializeTregsTemplate(config);

% Process each term in the formula
for t = 1:numel(terms)
    termStruct = processTerm(terms{t}, responseVar, config,t);
    tregs = appendTermToTregs(tregs, termStruct, config);
end
end

%% Core Parsing Functions
function [response, terms] = splitFormulaComponents(formulaStr)
% Split formula using pattern-based parsing
formulaParts = split(formulaStr, '~');
response = strtrim(formulaParts{1});
termList = split(formulaParts{2}, '+');
terms = cellfun(@strtrim, termList, 'UniformOutput', false);
end

function termStruct = processTerm(termStr, responseVar, config,ithterm)% Parse individual term using decomposition pattern
termStruct = struct(...
    'variables', {{}}, ...  % Initialize as cell array in scalar struct
    'batch_vars', {{}}, ...
    'continuous_vars', {{}}, ...
    'level', 'global', ...
    'hMu', {{}}, ...
    'hSigma', {{}}, ...
    'complex_handling', config.opt.complex_method ...
);


% Extract variables within f() using direct character matching
startParen = strfind(termStr, '(');
endParen = strfind(termStr, ')');
varStr = termStr(startParen+1:endParen-1);
variables = split(varStr, ',');


% Classify variables
termStruct.variables = variables;
termStruct.discret_vars = variables(config.hRangeMu{ithterm}==0);
termStruct.continuous_vars = setdiff(variables, termStruct.discret_vars);
termStruct.isbatch=intersect(variables, [config.batch]);

% Determine hierarchical level using precedence matching
if ~isempty(termStruct.discret_vars)
    termStruct.level = firstmatching(termStruct.discret_vars, config.level);
else
    termStruct.level = 'global';
end

% Configure bandwidth parameters from model settings
% [termStruct.hMu, termStruct.hSigma] = setBandwidths(...
%     numel(variables), ...
%     termStruct.level, ...
%     config.hRangeMu, ...
%     config.hRangeSigma);

termStruct.hRangeMu=config.hRangeMu{ithterm};
termStruct.hRangeSigma=config.hRangeSigma{ithterm};


end

%% Configuration Integration Helpers
function tregs = initializeTregsTemplate(config)
% Create tregs template from model configuration
tregs = cell(0,7);
tregsHeaders = {
    'level', 'condition', 'factor_mu', 'factor_sigma', ...
    'hList', 'options', 'calc_z'
    };
                                        
% Inherit complex handling options
defaultOptions = struct(...
    'is_plateau', true, ...
    'complex_method', config.opt.complex_method ...
    );%    'sep_complex', config.opt.sep_complex ...


tregs = cell2table(cell(0,7), 'VariableNames', tregsHeaders);
% tregs.options = {defaultOptions};
end

% function [hMu, hSigma] = setBandwidths(nVars, level, hRangeMu, hRangeSigma)
% % Set bandwidths based on hierarchical level and variable count
% switch level
%     case 'global'
%         hBaseMu = mean(hRangeMu(1:2));
%         hBaseSigma = mean(hRangeSigma(1:2));
%     otherwise
%         hBaseMu = mean(hRangeMu(3:4)) * 0.75;
%         hBaseSigma = mean(hRangeSigma(3:4)) * 0.6;
% end
% 
% % Create hList with proper dimensionality
% hMu = repmat({hBaseMu}, 1, nVars);
% hSigma = repmat({hBaseSigma}, 1, nVars);
% end

%% Intelligent Matching Utilities
function match = firstmatching(items, pool)
% Find first item in both arrays
for i = 1:numel(items)
    if any(strcmp(items{i}, pool))
        match = items{i};
        return;
    end
end
match = 'global';
end

function tregs = appendTermToTregs(tregs, termStruct, config)
% Append parsed term to tregs structure with validation
newEntry = cell(1,7);

% Translate term structure to tregs format
newEntry{1} = termStruct.level;
newEntry{2} = 'default';
newEntry{3} = termStruct.variables;
newEntry{4} = termStruct.variables;
newEntry{5} = termStruct.hMu;
newEntry{6} = config.opt;
newEntry{7} = true;

% Handle complex response types
if ~strcmp(config.opt.complex_method, 'seperate')
    newEntry{3} = appendComplexParts(newEntry{3}, config.resp_ystar);
end

tregs = [tregs; newEntry];
end
