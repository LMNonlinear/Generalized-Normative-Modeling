function tregs=tnureg_create(tregs)

% if ~isfield(tregs,'T')
%     tregs.T=[];
% end

if ~isfield(tregs,'factor_mu')
    tregs.factor_mu=[];
end
if ~isfield(tregs,'factor_sigma')
    tregs.factor_sigma=[];
end
if ~isfield(tregs,'condition') || isempty(tregs.condition) || isempty(tregs.condition{1})
    tregs.condition=[];
    tregs.num_cond=1;
else
    if ischar(tregs.condition)
        tregs.condition={tregs.condition};
    end
    tregs.num_cond=length(tregs.condition{2});
end
if ~isfield(tregs,'resp')
    tregs.resp=[];
end
if ~isfield(tregs,'regs')
    tregs.regs=[];
end
if ~isfield(tregs,'opt')
    tregs.opt=[];
end
% num_resp=length(tregs.resp);
if ~isfield(tregs.opt,'resp_mu') || isempty(tregs.opt.resp_mu)
    tregs.opt.resp_mu=append('mu',tregs.resp);
end
if ~isfield(tregs.opt,'resp_eps') || isempty(tregs.opt.resp_eps)
    tregs.opt.resp_eps=append('eps',tregs.resp);
end
if ~isfield(tregs.opt,'resp_sigma') || isempty(tregs.opt.resp_sigma)
    tregs.opt.resp_sigma=append('sigma',tregs.resp);
end
if ~isfield(tregs.opt,'resp_musigma') || isempty(tregs.opt.resp_musigma)
    tregs.opt.resp_musigma=append('musigma',tregs.resp);
end
if ~isfield(tregs.opt,'resp_z') || isempty(tregs.opt.resp_z)
    tregs.opt.resp_z=append('z',tregs.resp);
end

if ~isfield(tregs.opt,'calc_eps') || isempty(tregs.opt.calc_eps)
    tregs.opt.calc_eps=false;
end
if ~isfield(tregs.opt,'calc_sigma') || isempty(tregs.opt.calc_sigma)
    tregs.opt.calc_sigma=false;
end
if ~isfield(tregs.opt,'calc_musigma') || isempty(tregs.opt.calc_musigma)
    tregs.opt.calc_musigma=false;
end
if ~isfield(tregs.opt,'calc_z') || isempty(tregs.opt.calc_z)
    tregs.opt.calc_z=false;
end

tregs.regs=cell(tregs.num_cond,1);





end

