function [mq]=nureg_reg(regs,y)
assert(isfield(regs,'s'),'need precompute density');


[m_ft]=nureg_conv(regs,y);
mq = m_ft...
    ./(regs.s+eps);
% % lambda=0;
% lambda=0;
% % lambda=1e-3;
% % lambda=1e-4;
% mq = m_ft...
%     ./(regs.s+lambda*mean(regs.s,[1:regs.dx]));
% for k=1:size(y,2)
%     mq(:,:,:,k)= m_ft(:,:,:,k)...
%     ./(regs.s+lambda*mean(regs.s,[1:regs.dx]));
% end



end