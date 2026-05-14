function [regs]=update_nureg(regs,y,opt)

if ~isfield(opt,'y_corr_bandwith')
    opt.y_corr_bandwith=false;
end


%%
regs.opt=copy_field(opt,regs.opt);
regs.y=y;
dimy=ndims(y);
if dimy==2
    regs.dy=size(y,2);
else
    regs.dy=size(y,3);
    opt.y_corr_bandwith=true;
end
%% update h
[regs]=nureg_select_h(regs);
%% get smooth y
[regs]=nureg_y(regs);
%% select variables for the best bandwidth
[regs]=nureg_select_h(regs);
end