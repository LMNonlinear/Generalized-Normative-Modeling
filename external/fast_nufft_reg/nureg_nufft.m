function y_ft=nureg_nufft(regs,y)

% 
%% FGG_1d_type1
% y_ft=zeros([regs.L*ones(1,regs.dx),size(y,2)]);
% for iy=1:size(y,2)
%     if regs.dx==1
%         [yy]=FGG_1d_type1(y(:,iy),regs.knot,regs.L,regs.opt.ACC,regs.idx_list{1});
%         [y_ft(:,iy)] =fftshift((yy));%ifftshift fftshift
% %         [y_ft(:,iy),A] =fftshift((FGG_1d_type1(y(:,iy),regs.knot,regs.L,regs.ACC,regs.idx_list{1})));%ifftshift fftshift
%     elseif regs.dx==2
%         y_ft(:,:,iy) = fftshift((FGG_2d_type1(y(:,iy),regs.knot,[regs.L,regs.L,],regs.opt.ACC,regs.idx_list{1},regs.idx_list{2})));%,tau*ones(dx,1) x1grid(:),x2grid(:)
%     elseif regs.dx==3
%         y_ft(:,:,:,iy) = fftshift((FGG_3d_type1(y(:,iy),regs.knot,[regs.L,regs.L,regs.L],regs.opt.ACC,regs.idx_list{1},regs.idx_list{2},regs.idx_list{3})));
%     end
% end
%%
y_ft=nufftn_type1(regs.knot,y,regs.idx_list,-1,regs.opt.ACC,false,regs.opt.nufft_loop);
%%
for i=1:regs.dx
    y_ft=ifftshift(y_ft,i);
end
% y_ft=ifftshift(nufftn_type1(regs.knot,y,regs.idx_list,-1,regs.ACC,false));

% plot([abs(y_ft1(:,1)),abs(y_ft(:,1))])
end
