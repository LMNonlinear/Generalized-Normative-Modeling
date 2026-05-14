% setup.m - Initialize GNM-ToolBox MATLAB paths
%
% Adds core functions and NUFFT regression engine to the MATLAB search path.
% Also verifies that required MATLAB toolboxes are available.
%
% Usage:
%   >> run('setup.m')
%
% After setup, you can use:
%   >> mnhs = gnm('data.csv', 'y ~ s(age) | site[constant]');
%   >> mnhs = gnm_fit(mnhs);
%   >> T_new = gnm_predict(mnhs, 'new_data.csv');

rootDir = fileparts(mfilename('fullpath'));

% ── Add core paths ──
addpath(fullfile(rootDir, 'function'));
addpath(fullfile(rootDir, 'external', 'fast_nufft_reg'));

% ── Optional: ComBat baseline (if submodule is initialized) ──
combatDir = fullfile(rootDir, 'external', 'ComBatHarmonization', 'Matlab', 'scripts');
if exist(combatDir, 'dir')
    addpath(combatDir);
end

% ── Verify MATLAB version ──
v = ver('MATLAB');
if ~isempty(v)
    matlab_version = v(1).Version;
    if str2double(matlab_version(1:strfind(matlab_version,'.')-1)) < 9
        warning('GNM:setup:oldMatlab', ...
            'MATLAB R2022a or later recommended (detected: %s).', v(1).Release);
    end
end

% ── Verify required toolboxes ──
required_toolboxes = {'Statistics and Machine Learning Toolbox'};
installed = ver;
installed_names = {installed.Name};

missing = {};
for i = 1:numel(required_toolboxes)
    if ~any(strcmp(installed_names, required_toolboxes{i}))
        missing{end+1} = required_toolboxes{i}; %#ok<AGROW>
    end
end

if ~isempty(missing)
    warning('GNM:setup:missingToolbox', ...
        'Missing required toolbox(es): %s', strjoin(missing, ', '));
end

% ── Summary ──
fprintf('GNM-ToolBox initialized.\n');
fprintf('  Root      : %s\n', rootDir);
fprintf('  MATLAB    : %s\n', v(1).Release);
fprintf('  Functions : %d .m files\n', numel(dir(fullfile(rootDir, 'function', '*.m'))));
fprintf('\nQuick start:\n');
fprintf('  >> mnhs = gnm(''demo/test_synthetic_data.csv'', ''y ~ s(age) | site[constant]'');\n');
fprintf('  >> mnhs = gnm_fit(mnhs);\n');
fprintf('\nSee README.md for documentation.\n');

clear rootDir combatDir v matlab_version required_toolboxes installed installed_names missing i

