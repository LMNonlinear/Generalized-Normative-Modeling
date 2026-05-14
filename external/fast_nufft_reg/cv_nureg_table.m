function [T,regs]=cv_nureg_table(T,condition,factor_mu,factor_sigma,resp,hList,opt,regs)%,mnhs

%%
num_resp=length(resp);
if ~isfield(opt,'resp_mu')
    opt.resp_mu=append('mu',resp);
end
if ~isfield(opt,'resp_eps')
    opt.resp_eps=append('eps',resp);
end
if ~isfield(opt,'resp_sigma')
    opt.resp_sigma=append('sigma',resp);
end
if ~isfield(opt,'resp_musigma')
    opt.resp_musigma=append('musigma',resp);
end
if ~isfield(opt,'resp_z')
    opt.resp_z=append('z',resp);
end

if ~isfield(opt,'calc_eps')
    opt.calc_eps=false;
end
if ~isfield(opt,'calc_sigma')
    opt.calc_sigma=false;
end
if ~isfield(opt,'calc_musigma')
    opt.calc_musigma=false;
end
if ~isfield(opt,'calc_z')
    opt.calc_z=false;
end
if ~isfield(opt,'global_h')
    opt.global_h=false;
end
%%

if (opt.calc_musigma||opt.calc_z)
    if isnumeric(hList)
        hList={hList};
        hList{2}=hList{1};
    elseif iscell(hList) && numel(hList)==1
        hList{2}=hList{1};
    end
    dh(2)=size(hList{2},1);
end
% dh(1)=size(hList{1},1);



[subtable]=get_subtable(T,condition,[resp(:);factor_mu(:);factor_sigma(:)]);

if ischar(factor_mu)
    factor_mu={factor_mu};
end
dx=length(factor_mu);

%% get data
if ~isempty(factor_mu)
    x_mu=zeros(size(subtable,1),dx);
    for ix=1:dx
        x_mu(:,ix) = table2array(get_subtable(subtable,[],factor_mu(ix)));
    end
else
    x_mu=[];
end
if ~isempty(factor_sigma)
    x_sigma=zeros(size(subtable,1),dx);
    for ix=1:dx
        x_sigma(:,ix) = table2array(get_subtable(subtable,[],factor_sigma(ix)));
    end
else
    x_sigma=[];
end
y=table2array(get_subtable(subtable,[],resp));% y=gpuArray(table2array(get_subtable(subtable,[],resp)));



%% calculate
% mu
sz=size(y);
if opt.y_keep_all && size(y,2)~=num_resp
    y=reshape(y,sz(1),[],num_resp);
elseif opt.y_keep_all && size(y,2)==num_resp
    y=reshape(y,sz(1),1,num_resp);
end


opt.y_type_out='mean';
opt.y_gcv_dstd=opt.y_gcv_dstd_mu;
if ~isempty(x_mu)
    if nargin<6 || isempty(regs)
        clear regs
        [regs(1)]=cv_nureg_complex(x_mu,y,hList{1},opt);
    else
        if prod(ismember(hList{1},regs(1).h))
            opt.select_h=hList{1};
            [regs(1)]=update_nureg_complex(regs(1),y,opt);
        else
            [regs(1)]=cv_nureg_complex(x_mu,y,hList{1},opt);
        end
    end
    yhat=regs(1).yhat;
    regs(1).yhat=[];
else % constant mean: empty factor_mu => use grand mean per condition
    yhat=mean(y,1).*ones(size(subtable,1),1);
    regs(1).df_m=1;
end



% residual
yeps=[];
if opt.calc_eps || opt.calc_sigma|| opt.calc_musigma || opt.calc_z
    if opt.y_keep_all && ismatrix(y) && num_resp>1
        yeps=permute(y,[1,3,2])-yhat;
    else
        yeps=y-yhat;
    end
end

% variance
if opt.calc_sigma|| opt.calc_musigma || opt.calc_z
    ysigma=real(yeps).^2+1j*(imag(yeps)).^2;
end
if opt.global_h && ~opt.y_keep_all
    hList{2}=regs(1).gcv_yhat.h1se';
    opt.y_corr_bandwith=true;
end
% mean of variance
if opt.calc_musigma || opt.calc_z
    opt.y_type_out='variance';
    opt.y_gcv_dstd=opt.y_gcv_dstd_sigma;
    if ~isempty(x_sigma)
        if prod(ismember(hList{2},regs(1).h))
            opt.select_h=hList{2};
            [regs(2)]=update_nureg_complex(regs(1),ysigma,opt);
        else
            [regs(2)]=cv_nureg_complex(x_sigma,ysigma,hList{2},opt);
        end
        ymusigma=regs(2).yhat;
        regs(2).yhat=[];
    else % constant variance: empty factor_sigma => use grand variance per condition
        ymusigma=mean(ysigma,1).*ones(size(subtable,1),1);
        regs(2).df_m=1;
    end
    if any(real(ymusigma(:))<0)||any(imag(ymusigma(:))<0)
        error('ymusigma may not smooth well')
    end
end

% fisher's z score
if opt.calc_z
    if ~opt.y_keep_all
        yz=real(yeps)./sqrt(real(ymusigma))+...
            1j*(opt.sep_complex~=2).*imag(yeps)./(sqrt(imag(ymusigma))+eps);
    else
        yz=real(yeps)./sqrt(real(ymusigma))+...
            1j*permute((opt.sep_complex~=2),[1,3,2]).*imag(yeps)./(sqrt(imag(ymusigma))+eps);
    end
end
%% save to the table

T=asnarray2table(T,condition,opt.resp_mu,yhat);
if opt.calc_eps
    T=asnarray2table(T,condition,opt.resp_eps,yeps);
end
if opt.calc_sigma
    T=asnarray2table(T,condition,opt.resp_sigma,ysigma);
end
if opt.calc_musigma
    T=asnarray2table(T,condition,opt.resp_musigma,ymusigma);
end
if opt.calc_z
    T=asnarray2table(T,condition,opt.resp_z,yz);
end


end



%{

%% save into table
% mu
if ~opt.y_keep_all
    yhat=regs.yhat(:,1:num_resp);
    yhat(:,opt.idcomplex)=yhat(:,opt.idcomplex)+1j*regs.yhat(:,num_resp+[1:length(opt.idcomplex)]);
else
    yhat=regs.yhatall(:,:,1:num_resp);
    regs.yhatall(:,:,1:num_resp)=[];
    yhat(:,:,opt.idcomplex)=yhat(:,:,opt.idcomplex)+1j*regs.yhatall(:,1:length(opt.idcomplex));
    regs.yhatall=[];
end
T=asnarray2table(T,condition,opt.resp_mu,yhat);

% residual
yeps=[];
if opt.calc_eps
    yeps=y-yhat;
    T=asnarray2table(T,condition,opt.resp_eps,yeps);
end

% variance
ysigma=[];
if opt.calc_sigma
    if isempty(yeps)
        yeps=y-yhat;
    end
    ysigma=real(yeps).^2+1j*(imag(yeps)).^2;
    T=asnarray2table(T,condition,opt.resp_sigma,ysigma);
end

% mean of variance
if opt.calc_musigma
    if isempty(ysigma)
        if isempty(yeps)
            if ~opt.y_keep_all
                yeps=y-yhat;
            else
                yeps=permute(y,[1,3,2])-yhat;
            end
        end
        ysigma=real(yeps).^2+1j*(imag(yeps)).^2;
    end
    % smooth musigma
    opt.y_type_out='variance';
    opt.y_gcv_dstd=opt.y_gcv_dstd_sigma;
    if prod(ismember(hList{2},regs.h))
        opt.select_h=hList{2};
        [regs(2)]=update_nureg_complex(regs,ysigma,opt);
    else
        [regs(2)]=cv_nureg_complex(x,ysigma,hList{2},opt);
    end
    
    % mu
    if ~opt.y_keep_all
        yhat=regs.yhat(:,1:num_resp);
        yhat(:,opt.idcomplex)=yhat(:,opt.idcomplex)+1j*regs.yhat(:,num_resp+[1:length(opt.idcomplex)]);
    else
        yhat=regs.yhatall(:,:,1:num_resp);
        regs.yhatall(:,:,1:num_resp)=[];
        yhat(:,:,opt.idcomplex)=yhat(:,:,opt.idcomplex)+1j*regs.yhatall(:,1:length(opt.idcomplex));
        regs.yhatall=[];
    end
    T=asnarray2table(T,condition,opt.resp_musigma,ymusigma);
end
%}
% else
%     yhat=regs.yhat(:,1:num_resp);
%     yhat(:,opt.idcomplex)=yhat(:,opt.idcomplex)+1j*regs.yhat(:,num_resp+[1:length(opt.idcomplex)]);
%     % mu
%     [icol,T]=get_varidx(T,opt.resp_mu,regs.dtype);
%     T(irow,icol)= array2table(yhat);
%
%     % residual
%     yeps=[];
%     if opt.calc_eps
%         yeps=y-yhat;
%         icol=get_varidx(T,opt.resp_eps,regs.dtype);
%         T(irow,icol)= array2table(yeps);
%         regs.mse_all=norm(yeps,'fro');
%     end
%
%     % variance
%     ysigma=[];
%     if opt.calc_sigma
%         if isempty(yeps)
%             yeps=y-yhat;
%             regs.mse_all=norm(yeps,'fro');
%         end
%         ysigma=real(yeps).^2+1j*(imag(yeps)).^2;
%         [icol,T]=get_varidx(T,opt.resp_sigma,regs.dtype);
%         T(irow,icol)= array2table(ysigma);
%     end
%     % mean of variance
%     if opt.calc_musigma
%         if isempty(ysigma)
%             if isempty(yeps)
%                 yeps=y-yhat;
%                 regs.mse_all=norm(yeps,'fro');
%             end
%             ysigma=real(yeps).^2+1j*(imag(yeps)).^2;
%         end
%         [icol,T]=get_varidx(T,opt.resp_musigma,regs.dtype);
%         T(irow,icol)= array2table(ysigma);
%     end
% end

%{
    %% if save all y
    %mu
    yhatall=regs.yhatall(:,:,1:num_resp);
    regs.yhatall(:,:,1:num_resp)=[];
    yhatall(:,:,opt.idcomplex)=yhatall(:,:,opt.idcomplex)+1j*regs.yhatall(:,1:length(opt.idcomplex));
    regs.yhatall=yhatall;
    clear yhatall
    
    % eps
    if opt.calc_eps
        regs.yepsall=permute(y,[1,3,2])-regs.yhatall;
    end
    
    % sigma
    if opt.calc_sigma
        if ~isfield(regs,'yepsall') || isempty(regs.yepsall)
            regs.yepsall=permute(y,[1,3,2])-regs.yhatall;
        end
        regs.yepsall=real(regs.yepsall).^2+1j*(imag(regs.yepsall)).^2;
    end
%}








% if save all y
% if opt.y_keep_all
% %     opt.resp_mu_all=append('all',opt.resp_mu);
% %     opt.resp_res_all=append('all',opt.resp_res);
% %     opt.resp_sigma_all=append('all',opt.resp_sigma);
% %
%     yhatall=regs.yhatall(:,:,1:num_resp);
%     regs.yhatall(:,:,1:num_resp)=[];
%     yhatall(:,:,opt.idcomplex)=yhatall(:,:,opt.idcomplex)+1j*regs.yhatall(:,1:length(opt.idcomplex));
%     for i=1:num_resp
%         T.(opt.resp_mu_all{i})=yhatall(:,:,1);
%     end
% end











