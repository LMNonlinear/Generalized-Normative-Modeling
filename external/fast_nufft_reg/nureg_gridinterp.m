function [fpp,yhat]=nureg_gridinterp(regs,mq,xlist,method,opt)
% mq is value coresponding to the mesh grid
% mq=f(get_nd_graid_scatter(regs.xraw));
if nargin<3 || isempty(xlist)
    xlist=regs.xlist;
end
if nargin<4 || isempty(method)
    method ='griddedInterpolant';
end
if nargin<5 || isempty(opt)
    opt.Method='spline';
    opt.ExtrapolationMethod ='linear';
end
if isa(mq,'gpuArray')
    mq=gather(mq);
end

% fit grid value
sz=size(mq);
% figure(102)
% plot(xlist{1},reshape(mq,[regs.N*ones(1,regs.dx),prod(sz(regs.dx+1:end))]));

switch method
    case{'griddedInterpolant'}
        fpp = griddedInterpolant(xlist,reshape(mq,[regs.N*ones(1,regs.dx),prod(sz(regs.dx+1:end))]),opt.Method,opt.ExtrapolationMethod);
    otherwise
        error('unrecognized method')
end
% predict scatter value
if nargout>1
    xrawc=mat2cell(regs.xraw,regs.Tx,ones(1,regs.dx));
    yhat=squeeze(fpp(xrawc{:}));
    yhat=reshape(yhat,[regs.Tx,sz(regs.dx+1:end),1]);
end

end