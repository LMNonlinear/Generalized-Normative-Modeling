function Yq = nufftn_type1(x, y, xlist,iflag,acc,isdeconv,isloop)
%nufftn_type1
% x: the N*d dimenssion non-equispaced locations (knots) where y lies on
% y: a N*dy dimenssion non-equispaced function values
% M: a dx1 vector denoting number of FFT grid in each dimenssion
% df: frequency spacing
% eps: tolerance for NUFFT
% iflag: sign for the exponential (0 means :math:`-i`, greater than 0 means :math:`+i`)

[M,dx]=size(x);
[~,dy]=size(y);
% ndcolon(1:dx) = {':'};

if nargin<3 || isempty(xlist)
    N=ceil(M^(1/dx))*ones(dx,1);
    xlist=[];
else
    if dx>1 || iscell(xlist)
        N=cellfun(@(x) length(x),xlist);
    else
        N=length(xlist);
    end
end
% N=permute(N,[3,1,2]);

% [fgrid] = nufftfreqs(N);
if nargin<4 || isempty(iflag)
    iflag=1;
end
if nargin<5 || isempty(acc)
    acc=6;
end
if nargin<6 || isempty(isdeconv)
    isdeconv=true;
end
if nargin<7 || isempty(isloop)
    isloop=false;
end
% Fast Non-Uniform Fourier Transform with MATLAB
% Msp denote the number of grid points to which the "spreading" will be
%   accounted for in each direction, depends on the descire accuracy
% Mr is the length of the FFTs(oversampled) Mr=ratio*M
% tau is spacing of grid(ftau)?


[Msp, Mr, tau,R] = compute_grid_params(N, acc);

% x=permute(x,[1,3,2]);
% Msp=permute(Msp,[3,2,1]);
% Mr=permute(Mr,[3,2,1]);


% Construct the convolved grid
% ftau is the grid data after convolved
% hx is spacing of grid(ftau)
% xmod is the the source x location in 0~2pi
% mm is spreading grid(2Msp *1)
% m is the grid index where source xmod should be
% mpmm is m+mm with broadcasting N*Msp matrix denoting the spredding neigbourhood for each source
% indx is make sure all mpmm don't exceed Mr

hx   = 2*pi./Mr;
xmod = scale_knots(x,N,xlist);
m    = round(xmod ./ hx);
Msp2=Msp*2+1;
%for broadcasting
m=permute(m,[1,3,2]);
Mr=permute(Mr(:),[3,2,1]);
hx=permute(hx(:),[3,2,1]);
xmod=permute(xmod,[1,3,2]);
mpmm=m;
for i=1:dx
    m=mpmm;
    mm=zeros(1,Msp2(i),dx);
    mm(1,:,i) = round(linspace(-Msp(i),Msp(i),Msp2(i)).');
    mpmm = reshape(m + mm,[],1,dx);
end
mpmm=reshape(mpmm,M,prod(Msp2),dx);




spread = heat_kernel(xmod-hx.* mpmm,tau);
y=permute(y,[1,3,2]);
ysp = y .* spread;
clear y
ysp=reshape(ysp,[],dy);
%% xindx
if isloop
    dtype=min_numerical_type(0,max(Mr));
    %     dtype='int64';
    fun_dtype=str2func(dtype);
    xindx = fun_dtype(mod(mpmm, Mr)+1);
    clear mpmm
    xindx = reshape(xindx,[],dx);
else
    dtype=min_numerical_type(0,max(Mr));
    %     dtype='int64';
    fun_dtype=str2func(dtype);
    xindx = fun_dtype(mod(mpmm, Mr)+1);
    clear mpmm
    xindx = reshape(xindx,[],dx);
    if dy>1
        xindx(:,dx+1)=1;
    end
    xindx = repmat(xindx,[dy,1]);
    if dy>1
        xindx(:,dx+1)=reshape(fun_dtype(1:dy).*ones(M*prod(Msp2),1,dtype),[],1);%kron((1:dy).',ones(M*prod(Msp2),1));
    end
end
%% accumarray
mu=bit_convert(underlyingType(ysp),0);
if isloop
    
    %%
    sz=Mr(:)';
    Ftau=zeros([sz,dy]);
    if dx==1
        sz=[Mr(:)',1];
    end
    for iy=1:dy
        Ftau(:,:,iy) = accumarray(xindx, ysp(:,iy), sz,[],mu);
    end

    %%
    %     if dx==1
    %         sz=[Mr(:)',1];
    %     else
    %         sz=Mr(:)';
    %     end
    %     Ftau = arrayfun(@(iy) accumarray(xindx, ysp(:,iy), sz,[],mu),1:dy,'UniformOutput',false);
    %     % Ftau = arrayfun(@(iy) accumarray_source(xindx, ysp(:,iy), sz,[],mu),1:dy,'UniformOutput',false);
    %     Ftau = cell2mat(permute(Ftau,[1,3,2]));
    %
    %%
    %     sz=fun_dtype([Mr(:)',dy]);
    %     Ftau = addat_mex(xindx, ysp, sz);
    
    clear ysp xindx
    
    
else
    if dx==1 || dy>1
        sz=[Mr(:);dy]';
    else
        sz=[Mr(:)]';
    end
    Ftau = (accumarray(xindx, ysp(:), sz,[],mu));
    clear ysp xindx
end
%% fft
if iflag < 0
    for ix=dx:-1:1
        Ftau = fftshift(fft(Ftau, [], ix),ix);
    end
    Ftau=Ftau/(M*(R^dx));
    
else
    for ix=dx:-1:1
        Ftau = ifft(ifftshift(Ftau,ix), [], ix);
    end
end

% plot([ifftshift(A.f_tau(:,1))-(ftau(:,1))])
% plot([A.f_tau,ifftshift(ftau)])
% figure(101);clf;plot(real(Ftau)*(M*(R^dx)));
% figure(101);hold on;plot(real(A.F_tau))
% figure(101);clf;plot(real(Ftau));
% figure(101);hold on;plot(real(A.F_tau))

%% index for output
q=(Mr(:)-N(:))/2;
q=[ceil(q)+1,ceil(q)+N(:)];

sz=[Mr(:);dy]';
idx=get_patch_index(q,sz);
Ftau=Ftau(idx);
Ftau=reshape(Ftau,[N(:);dy].');
% plot([A.F_tau,ifftshift(Ftau)])
% plot(abs([(A.F_tau(:,1))/(M*(R^dx)),(Ftau(:,1))]))
%% Deconvolve the grid using convolution theorem, 2004 formula 11
if isdeconv
    %         Kninv=sqrt((pi^dx)/prod(tau(:)))* exp(sum(permute(tau(:),[2:dx+1,1]).* fgrid.^2,dx+1));
    %         Kninv=Kninv./(eps*mean(Kninv,'all'));
    %         Yq =  Kninv .* Ftau;
    Kn=sqrt(prod(tau(:))/(pi^dx))* exp(-sum(permute(tau(:),[2:dx+1,1]) .* nufftfreqs(N).^2,dx+1));
    %     Kninv=1./(Kn+1e8*max(Kn,[],1:dx));
    Kninv=1./Kn;
    %     Kninv=1./(Kn+mean(Kn,'all'));
    Yq =  Kninv .* Ftau;
    
else
    Yq =  Ftau;
end

end

function x=heat_kernel(x,h)

h=permute(h(:),[3:-1:1]);
x=exp( -sum((x.^2)./(4*h),3));
end