function [m]=nureg_conv(regs,y)
sz=size(y);
dy=sz(end);
Ny=sz(1);
if Ny~=regs.L%since L must large than N in the setting in nureg_create
    flag_transformed=false;
else
    flag_transformed=true;
    dy=size(y,regs.dx+1);
    %     kdf=gpuArray(single(regs.kdf));
end
% msz=bytesize(numel(regs.kdf),'b',regs.kdf(1));
% isusegpu=test_gpu([],2*msz.b);
% if isusegpu
%     y=gpuArray(single(y));
% else
kdf=regs.kdf;
% end
%% nureg_conv
% nufft version for gridless samples convolution
% if ~regs.opt.y_corr_bandwith || ~flag_transformed
%     if flag_transformed
%         m_ft = kdf.*permute(y,[1:regs.dx,regs.dx+2,regs.dx+1]);
%     else
%         %     if ~regs.opt.y_corr_bandwith
%         m_ft = kdf.*permute(nureg_nufft(regs,y),[1:regs.dx,regs.dx+2,regs.dx+1]);
%         %     else
%         %         m_ft = kdf.*reshape(nureg_nufft(regs,y),[regs.L*ones(1,regs.dx),regs.dh,regs.dy]);
%         %     end
%     end
%     sz=size(m_ft);
%     m=reshape(m_ft,[regs.L*ones(1,regs.dx),regs.dh*dy]);
% else
%     if flag_transformed
%         m_ft = kdf.*reshape(y,[regs.L*ones(1,regs.dx),regs.dh,regs.dy]);
%     else
%         error('not supoort, only for y')
%     end
%     m=reshape(m_ft,[regs.L*ones(1,regs.dx),regs.dh*dy],[regs.L*ones(1,regs.dx),regs.dh*dy]);
% end
if  ~flag_transformed
    m_ft = kdf.*permute(nureg_nufft(regs,y),[1:regs.dx,regs.dx+2,regs.dx+1]);
    sz=size(m_ft);
    m=reshape(m_ft,[regs.L*ones(1,regs.dx),regs.dh*dy]);    
else
    if ~regs.opt.y_corr_bandwith
        m_ft = kdf.*permute(y,[1:regs.dx,regs.dx+2,regs.dx+1]);
        sz=size(m_ft);
        m=reshape(m_ft,[regs.L*ones(1,regs.dx),regs.dh*dy]);
    else
        m_ft = kdf.*reshape(y,[regs.L*ones(1,regs.dx),regs.dh,dy/regs.dh]);
        m=reshape(m_ft,[regs.L*ones(1,regs.dx),dy]);
    end
end

%%
clear kdf y m_ft
for ix=regs.dx:-1:1
    m = ifft(m, [], ix);
end
if ~flag_transformed || ~regs.opt.y_corr_bandwith
    m=reshape(real(m),[regs.L*ones(1,regs.dx),regs.dh,dy]);%real(m)
else
    m=reshape(real(m),[regs.L*ones(1,regs.dx),regs.dh,dy/regs.dh]);%real(m)
end
%%
idx_mq=get_patch_index(regs.qout,sz);
m=m(idx_mq);

if ~flag_transformed || ~regs.opt.y_corr_bandwith
    m=reshape(m,[regs.N*ones(1,regs.dx),regs.dh,dy]);    
else
    m=reshape(m,[regs.N*ones(1,regs.dx),regs.dh,dy/regs.dh]);
end
% if isusegpu
%     m=bit_convert(regs.ACC,gather(m));
% end
end



