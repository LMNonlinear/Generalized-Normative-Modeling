function mnhs = gnm_run_pipeline(mnhs)
% gnm_run_pipeline  Execute the configured GNM pipeline.
%
%   mnhs = gnm_run_pipeline(mnhs)
%
%   Dispatches to the pipeline function stored in mnhs.mnhfun (typically
%   @gnm_compute_zscores), manages diary logging, and saves results.
%
%   See also: gnm_fit, gnm_build_default_pipeline, gnm_compute_zscores

close all;
dbstop if error

%% record
diary off
diary(get_save_path(mnhs.path_csv, '_zmap', {mnhs.tag, 'diary'}, '.txt'))
mnhs.time_start = datetime;
mnhfun = mnhs.mnhfun;
mnhs = mnhfun(mnhs);

%% only for model compare
if isfield(mnhs, 'only_criteria') && mnhs.only_criteria
    warning('not save T and tregs, to save time')
    mnhs.T = [];
    mnhs.tregs = [];
end

%% change the batch back to the original from the reference batch
if isfield(mnhs, 'reRefBatch') && ~isempty(mnhs.reRefBatch)
    idx = ~isemptycell(mnhs.Ttest.orignalBatch);
    orignalBatch = setdiff(unique(mnhs.Ttest.orignalBatch(idx)), {''});
    for i = 1:size(mnhs.reRefBatch, 1)
        mnhs.Ttest = asnarray2table(mnhs.Ttest, {'orignalBatch', orignalBatch}, mnhs.batch, mnhs.reRefBatch(i,1));
    end
end

%% save result
path_struct = save_struct(mnhs, mnhs.path_csv, '_zmap', mnhs.tag);

disp(['finish ', mfilename, '-', func2str(mnhfun), '-', mnhs.tag, ' for ', path_struct]);
mnhs.time_finish = datetime;
diary off
