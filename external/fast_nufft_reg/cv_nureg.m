function [regs]=cv_nureg(x,y,h,opt)
if nargin<4||isempty(opt)
    opt=[];
end

%% initialization
[regs]=nureg_create(x,y,h,opt);

%% get density
[regs.s]=nureg_conv(regs,ones(regs.Tx,1,regs.dtype))+eps;

%% get degree of freedom(dof) with randmomized method
[regs]=nureg_dof(regs);

%% get smooth y
[regs]=nureg_y(regs);

%% select variables for the best bandwidth
[regs]=nureg_select_h(regs);

%% compact to save memory
[regs]=nureg_compact(regs);
% warning('could uncomment this line to save memory')
end





