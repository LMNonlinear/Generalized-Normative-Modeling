function [regs]=cv_nureg_complex(x,y,h,opt)
dh=size(h,1);

%%
dimy_input=ndims(y);
if dimy_input==2
    [~,dy]=size(y);
    y=[real(y),imag(y)];
    y(:,dy+find(opt.sep_complex==2))=[];
else
    [~,dhy,dy]=size(y);
    y=cat(3,real(y),imag(y));
    y(:,:,dy+find(opt.sep_complex==2))=[];
end
sz=size(y);
if isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && dimy_input~=2 && size(opt.inherit_bandwidth,1)==dh
    y=reshape(y(:,opt.inherit_bandwidth(:),:),sz(1),dh,size(opt.inherit_bandwidth,2)*sz(3));
elseif  isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && dimy_input~=2 && size(opt.inherit_bandwidth,1)~=dh
    y=reshape(y,sz(1),[],sz(3));
end

%%
[regs]=cv_nureg(x,y,h,opt);
% if opt.y_keep_all
%     regs.yhat=reshape(regs.yhat,sz(1),dh,sz(3));
% end



if isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && dimy_input~=2 && size(opt.inherit_bandwidth,1)==dh
    regs.yhat=reshape(regs.yhat,sz(1),numel(opt.inherit_bandwidth),sz(3));
    regs.yhat=regs.yhat(:,opt.inherit_bandwidth(:),:);
elseif  isfield(opt,'inherit_bandwidth')&& ~isempty(opt.inherit_bandwidth) && dimy_input~=2 && size(opt.inherit_bandwidth,1)~=dh
    regs.yhat=reshape(regs.yhat,sz(1));
end




%%
dimy_output=ndims(regs.yhat);
if dimy_output==2
%     yhat=regs(1).yhat(:,1:num_resp);
    regs.yhat(:,opt.idcomplex)=regs.yhat(:,opt.idcomplex)+1j*regs.yhat(:,dy+[1:length(opt.idcomplex)]);
    regs.yhat(:,dy+[1:length(opt.idcomplex)])=[];
else
%     yhat=regs(1).yhat(:,:,1:num_resp);
%     regs(1).yhat(:,:,1:num_resp)=[];
%     yhat(:,:,opt.idcomplex)=yhat(:,:,opt.idcomplex)+1j*regs(1).yhat(:,:,1:length(opt.idcomplex));
%     regs(1).yhat=[];
    regs.yhat(:,:,opt.idcomplex)=regs.yhat(:,:,opt.idcomplex)+1j*regs.yhat(:,:,dy+(1:length(opt.idcomplex)));
    regs.yhat(:,:,dy+(1:length(opt.idcomplex)))=[];
end




% % deal with complex mu
% if ismatrix(regs(1).yhat) && num_resp~=1
%     yhat=regs(1).yhat(:,1:num_resp);
%     yhat(:,opt.idcomplex)=yhat(:,opt.idcomplex)+1j*regs(1).yhat(:,num_resp+[1:length(opt.idcomplex)]);
% else
%     yhat=regs(1).yhat(:,:,1:num_resp);
%     regs(1).yhat(:,:,1:num_resp)=[];
%     yhat(:,:,opt.idcomplex)=yhat(:,:,opt.idcomplex)+1j*regs(1).yhat(:,:,1:length(opt.idcomplex));
%     regs(1).yhat=[];
% end

% deal with complex mu
% plot(x(:,1),y(:,1,1),'.');hold on;plot(x(:,1),regs.yhat(:,1,1))



end




