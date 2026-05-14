function [kd,ihbad,dh,dh_input]=nureg_kdf(x,h,flag_power2,Nratio)

[Tx,dx]=size(x);
x_min=min(x);
x_max=max(x);
x_scale=x_max-x_min;
N=ceil(Nratio*(Tx^(1/dx)));
% length of samples after padding zeros
if flag_power2
    % zero-padding improves accuracy and speed
    L=2^nextpow2(2*N-1);
%     L=2^nextpow2(1.2*N-1);
else
    L=2*N-1;
end
qin=zeros(2,dx);
qin(1,:)=round(L/2-double((abs(x_min)./x_scale)*N));
qin(2,:)=(qin(1,:)+N-1);

% bandwidth
if nargin<3 || isempty(h)
    h=(4/(dx+2)/Tx)^(1/(dx+4))*ones(1,dx,class(y));
elseif (isvector(h)|| isscalar(h) )&& length(h)~=dx    
    h=repmat(h(:),[1,dx]);
elseif size(h,2)~=dx && size(h,1)==dx
    h=h.';
elseif (isscalar(h) && dx==1) || (ismatrix(h) && dx==size(h,2))    
else
    error('error h dimenssion')
end

[dh]=size(h,1);
dh_input=dh;

%% get frequency response of kernel
% generate ndgrid with the zscored data range
xgrid=get_ndgrid_scatter(x,'array',N);
% get the time domain kernel, apply kernel function to the zscored ndgrid
kd=epan_kernel(xgrid,h);

ihbad=squeeze(~sum(kd,(1:dx)));
ihbad=ihbad(:);
ndcolon(1:dx) = {':'};
if sum(ihbad) && sum(~ihbad)>0
    kd(ndcolon{:},ihbad)=[];
    dh=sum(~ihbad);
    %     h(ihbad,:)=[];
elseif prod(ihbad)
    % assert(sum(kd,[1:dx]),['bandwidth h = ', num2str(h) ,'  is too small'])
    warning(['all bandwidth h  is too small. set as empirical h = ',num2str((4/(dx+2)/Tx)^(1/(dx+4)))])%= ', num2str(h(ihbad,:)) ,'
    h=(4/(dx+2)/Tx)^(1/(dx+4));
    dh=1;
    kd=epan_kernel(xgrid,h);
    %     kd_cv=epan_kernel(x,h);
end

% padding zeros
kd=padarray(kd,qin(1,:).*ones(1,dx),0,'pre');
kd=padarray(kd,(L-qin(2,:)-1).*ones(1,dx),0,'post');
% kdf=zeros(size(kd),class(x));
% n dimenssion fft
% for i=1:dh
% %     kdf(ndcolon{:},i)=fftshift(fftn(ifftshift(kd(ndcolon{:},i))));
%     kdf(ndcolon{:},i)=(fftn((kd(ndcolon{:},i))));
% end
for i=1:dx
    kd=fft(kd,[],i);
end



end