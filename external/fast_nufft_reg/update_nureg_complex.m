function [regs]=update_nureg_complex(regs,y,opt)

if ~isfield(opt,'y_corr_bandwith')
    opt.y_corr_bandwith=false;
end


dimy=ndims(y);
if dimy==2
    [~,dy]=size(y);
    y=[real(y),imag(y)];
    y(:,dy+find(opt.sep_complex==2))=[];
%     if opt.y_corr_bandwith
%         y=permute(y,[1,3,2]);
%     end    
else
    [~,~,dy]=size(y);
    y=cat(3,real(y),imag(y));
    y(:,:,dy+find(opt.sep_complex==2))=[];
    opt.y_corr_bandwith=true;
end
sz=size(y);
dh=regs.dh;
if isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && size(opt.inherit_bandwidth,1)==dh
    y=reshape(y(:,opt.inherit_bandwidth(:),:),sz(1),dh,size(opt.inherit_bandwidth,2)*sz(3));
elseif  isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && size(opt.inherit_bandwidth,1)~=dh
    y=reshape(y,sz(1),[],sz(3));    
end
%%
[regs]=update_nureg(regs,y,opt);
% if opt.y_keep_all
%     regs.yhat=reshape(regs.yhat,sz(1),dh,[]);
% end
if isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && size(opt.inherit_bandwidth,1)==dh
    regs.yhat=reshape(regs.yhat,sz(1),numel(opt.inherit_bandwidth),sz(3));
    regs.yhat=regs.yhat(:,opt.inherit_bandwidth(:),:);
elseif  isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && size(opt.inherit_bandwidth,1)~=dh
    regs.yhat=reshape(regs.yhat,sz(1));
end

%%
dimy_output=ndims(regs.yhat);
if dimy_output==2
    regs.yhat(:,opt.idcomplex)=regs.yhat(:,opt.idcomplex)+1j*regs.yhat(:,dy+[1:length(opt.idcomplex)]);
    regs.yhat(:,dy+[1:length(opt.idcomplex)])=[];
else
    regs.yhat(:,:,opt.idcomplex)=regs.yhat(:,:,opt.idcomplex)+1j*regs.yhat(:,:,dy+(1:length(opt.idcomplex)));
    regs.yhat(:,:,dy+(1:length(opt.idcomplex)))=[];
end

end