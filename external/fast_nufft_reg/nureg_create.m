function [regs]=nureg_create(x,y,h,opt)




%% options
if nargin<4||isempty(opt)
    opt=[];
end
if ~isfield(opt,'y_gcv_metric')
    opt.y_gcv_metric=[];
end
if ~isfield(opt,'y_gcv_dstd')
    opt.y_gcv_dstd=0;
end
if ~isfield(opt,'y_grid_method')
    opt.y_grid_method='griddedInterpolant';
end
if ~isfield(opt,'y_grid_opt')
    opt.y_grid_opt.Method='spline';
    opt.y_grid_opt.ExtrapolationMethod='linear';
end
if ~isfield(opt,'num_dof_sample')
    opt.num_dof_sample=10;
end
if ~isfield(opt,'select_h')
    opt.select_h=[];
end
if ~isfield(opt,'y_keep_all')
    opt.y_keep_all=false;
end
if ~isfield(opt,'isgpu')
    opt.isgpu=false;
end
if ~isfield(opt,'nufft_loop')
    opt.nufft_loop=false;
end

% 'mean','residual' %'variance'
if ~isfield(opt,'y_type_out')
    opt.y_type_out='mean';%%'residual' %'variance'
end
if ~isfield(opt,'flag_power2')
    opt.flag_power2='true';%%'residual' %'variance'
end
if ~isfield(opt,'ACC')
    opt.ACC=2;%%'residual' %'variance'
end
if ~isfield(opt,'Nratio')
    opt.Nratio=2;%%'residual' %'variance'
end
if ~isfield(opt,'calc_dof')
    opt.calc_dof=true;
end

if ~isfield(opt,'compact')
    opt.compact=true;
end



%%
[Tx,dx]=size(x);
if ismatrix(y)
    [Ty,dy]=size(y);
    opt.y_corr_bandwith=false;
elseif ndims(y)==3
    [Ty,dhy,dy]=size(y);
    opt.y_corr_bandwith=true;
end
assert(Tx==Ty,'Length of x and y has to be the same');

%% initialize
dtype=class(y);
[x,y]=bit_convert(opt.ACC,x,y);
if strcmp(dtype,'gpuArray') %#ok<STISA>
    x=gpuArray(x);
    y=gpuArray(y);
end
dtype=class(y);
N=ceil(opt.Nratio*(Tx^(1/dx)));
ndcolon(1:dx) = {':'};
% length of samples after padding zeros
if opt.flag_power2
    % zero-padding improves accuracy and speed
    L=2^nextpow2(2*N-1);
%     L=2^nextpow2(1.2*N-1);
else
    L=2*N-1;
end



%% work on scaled space
xraw=x;
[x,x_mean,x_std]=zscore(x);
x_min=min(x);
x_max=max(x);
x_scale=x_max-x_min;
x_delta=x_scale./N;

xraw_min=x_min.*x_std+x_mean;
xraw_max=x_max.*x_std+x_mean;
% xraw_scale=x_scale.*x_std;
xraw_delta=x_delta.*x_std;


%position padding points
qin=zeros(2,dx);
q=(L/2-double((abs(x_min)./x_scale)*N));
qin(1,:)=round(q);
qin(2,:)=(qin(1,:)+N-1);


%% get frequency responce of kernel
% todo: could be improved with symbolic calculation
[kdf,ihbad,dh,dh_input]=nureg_kdf(x,h,opt.flag_power2,opt.Nratio);  %% todo: it is better put outside. for local linear
if opt.y_corr_bandwith && dh~=dhy && dhy~=1
    error('Corresponding reggression dimenssion not match')
end

%% get the knots(sample point position)
% get the index where raw data is, the rest are padded to zeros
idx_knot=[qin(1):qin(2)].';
idx_range=linspace(-1/2,1/2-1/L, L)';
knot_scale=idx_range(idx_knot(end))-idx_range(idx_knot(1));
knot=(x-x_min)./x_scale.*knot_scale+idx_range(idx_knot(1))*ones(1,dx);
Msp = compute_grid_params(Ty, opt.ACC);

%% index for output
qout=(L-N*ones(dx,1))/2;
qout=[ceil(qout)+1,ceil(qout)+N];

idx_mq=get_patch_index(qout,[L*ones(1,dx),dh]);
idx_list(1:dx)={idx_range};
%% using raw space for fitting grid parametric methods
xlist=multispace(min(xraw),max(xraw),N).';
xlist=mat2cell(xlist,N,ones(1,dx));


%% save into regs
regs.xraw=xraw;
regs.x=x;
regs.y=y;
regs.Tx=Tx;
regs.Ty=Ty;
regs.dx=dx;
regs.dy=dy;
regs.dtype=dtype;
regs.Msp=Msp;

% regs.ACC=opt.ACC;
% regs.flag_power2=opt.flag_power2;
regs.L=L;
regs.N=N;
regs.ndcolon=ndcolon;

regs.x_min=x_min;
regs.x_max=x_max;
regs.x_scale=x_scale;
regs.x_mean=x_mean;
regs.x_std=x_std;
regs.x_delta=x_delta;

regs.h=h;
regs.kdf=kdf;
regs.ihbad=ihbad;
regs.dh=dh;
regs.dh_input=dh_input;

regs.qin=qin;
regs.qout=qout;
regs.idx_mq=idx_mq;


regs.knot_scale=knot_scale;
regs.knot=knot;
regs.idx_list=idx_list;

regs.xlist=xlist;
regs.opt=opt;



end