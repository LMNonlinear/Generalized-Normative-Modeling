function [T,crit]=tnureg_criteria(T,y,tregs,opt)
%%
% if no musigma or zmap, require: resp eps
% if no musigma or zmap, require: resp eps musigma

num_level=size(tregs,1);
dy=numel(y);
dyreal=dy;
dycomplex=length(tregs{1}.opt.idcomplex);
n=size(T,1);
dh=size(T.(tregs{num_level}.opt.resp_eps{1})(1,:),2);% dh=tregs{num_level}.regs{1}.dh;

%% sigma accumulate
ysigmaaccum=cell(num_level,1);
for i=1:num_level
    ysigmaaccum{i}=append('sigmaaccum',tregs{i}.resp);
end
for i=1:num_level
    T=asnarray2table(T,[],ysigmaaccum{i},repmat(1+1i,[1,dh,dy]));
end

for i=1:num_level
    for j=1:i % For zmap based.
        if tregs{j}.opt.calc_musigma
            func= @(sigma,sigmaccum) sigma.*sigmaccum;
            T = complexfun(@tablefun,opt.complex_method,func,T,ysigmaaccum{i},[],...
                tregs{j}.opt.resp_musigma,...
                ysigmaaccum{i});
        end
    end
end

for i=1:num_level
    if ~tregs{i}.opt.calc_musigma % For residual based.
        T=clone_table_var(T,tregs{i}.opt.resp_musigma,tregs{i}.opt.resp_sigma);
        func=@(sigma) mean(sigma,1);
        T=tablefun(func,T,tregs{i}.opt.resp_musigma,[],tregs{i}.opt.resp_musigma);
    end
end


%% compact to save memory
for i=1:num_level
    T=asnarray2table(T,[],tregs{i,1}.resp,[]);
end
if opt.compact.musigma
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_sigma,[]);
    end
end
%% sse (standardrized square error)(sigmaaccumpre*sigma/sigmaaccum=sigma/sigma or simga/musigma->n)
sse=zeros(num_level,dh);
ysse=cell(num_level,1);
for i=1:num_level
    ysse{i}=append('sse',tregs{i}.resp);
end
for i=1:num_level
    func=@(yeps,sigma) (real(yeps).^2)./(real(sigma))+...
        1i*permute((tregs{i,1}.opt.sep_complex~=2),[1,3,2]).*(imag(yeps).^2)./(imag(sigma)+eps);
    sse_temp=func(reshape(table2array(get_subtable(T,[],tregs{i}.opt.resp_eps)),[],dh,dy),...
        reshape(table2array(get_subtable(T,[],tregs{i}.opt.resp_musigma)),[],dh,dy));
    sse(i,:)=sum(sse_temp,[1 3]);
end
%% compact to save memory
if opt.compact.eps
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_eps,[]);
    end
end

%% df
df=zeros(num_level,dh);
for i=1:num_level
    if i>1
        df(i,:)=df(i-1,:)+df(i,:);
    end
    for iregs=1:length(tregs{i}.regs)
        df_m_raw = tregs{i}.regs{iregs}(1).df_m(:)';
        if ~isfield(tregs{i}.opt,'inherit_bandwidth') || isempty(tregs{i}.opt.inherit_bandwidth)
            % When df_m has more entries than dh (GCV selected subset),
            % use only the GCV-selected value (index from gcv_yhat if available)
            if numel(df_m_raw) > dh && isfield(tregs{i}.regs{iregs}(1),'gcv_yhat')
                sel = tregs{i}.regs{iregs}(1).gcv_yhat.id1se;
                df_m_raw = df_m_raw(sel);
            elseif numel(df_m_raw) > dh
                df_m_raw = df_m_raw(1:dh);
            end
            df(i,:)=df(i,:)+df_m_raw;
        else
            df_m_temp=zeros(1,dh);
            for ih=1:size(tregs{i}.opt.inherit_bandwidth,2)
                if numel(df_m_raw) >= max(tregs{i}.opt.inherit_bandwidth(:,ih))
                    df_m_temp(tregs{i}.opt.inherit_bandwidth(:,ih))=df_m_raw(tregs{i}.opt.inherit_bandwidth(:,ih));
                else
                    df_m_temp = df_m_temp + mean(df_m_raw);
                end
            end
            df(i,:)=df(i,:)+df_m_temp;
        end
    end
    if tregs{i}.opt.calc_musigma
        for iregs=1:length(tregs{i}.regs)
            df_m_raw = tregs{i}.regs{iregs}(2).df_m(:)';
            if ~isfield(tregs{i}.opt,'inherit_bandwidth') || isempty(tregs{i}.opt.inherit_bandwidth)
                if numel(df_m_raw) > dh && isfield(tregs{i}.regs{iregs}(2),'gcv_yhat')
                    sel = tregs{i}.regs{iregs}(2).gcv_yhat.id1se;
                    df_m_raw = df_m_raw(sel);
                elseif numel(df_m_raw) > dh
                    df_m_raw = df_m_raw(1:dh);
                end
                df(i,:)=df(i,:)+df_m_raw;
            else
                df_m_temp=zeros(1,dh);
                for ih=1:size(tregs{i}.opt.inherit_bandwidth,2)
                    if numel(df_m_raw) >= max(tregs{i}.opt.inherit_bandwidth(:,ih))
                        df_m_temp(tregs{i}.opt.inherit_bandwidth(:,ih))=df_m_raw(tregs{i}.opt.inherit_bandwidth(:,ih));
                    else
                        df_m_temp = df_m_temp + mean(df_m_raw);
                    end
                end
                df(i,:)=df(i,:)+df_m_temp;
            end
        end
    end
end
%% log sigma
sigma=zeros(num_level,dh);
logsigma=zeros(num_level,dh);
for i=1:num_level
    if i==1
        sigmaaccum=reshape(table2array(get_subtable(T,[],tregs{i}.opt.resp_musigma)),[],dh,dy);
    else
        sigmaaccum=real(reshape(table2array(get_subtable(T,[],ysigmaaccum{i-1})),[],dh,dy)).*real(reshape(table2array(get_subtable(T,[],tregs{i}.opt.resp_musigma)),[],dh,dy))+...
            1i*imag(reshape(table2array(get_subtable(T,[],ysigmaaccum{i-1})),[],dh,dy)).*imag(reshape(table2array(get_subtable(T,[],tregs{i}.opt.resp_musigma)),[],dh,dy));
    end
    sigma(i,:)=sum((real(sigmaaccum))+1i*permute((tregs{i,1}.opt.sep_complex~=2),[1,3,2]).*(imag(sigmaaccum)),[1,3]);
    logsigma(i,:)=sum(log(real(sigmaaccum))+1i*permute((tregs{i,1}.opt.sep_complex~=2),[1,3,2]).*log(imag(sigmaaccum)+eps),[1,3]);
%     logsigma(i,:)=log(sum(real(sigmaaccum),[1,3]))+1i*log(sum(permute((tregs{i,1}.opt.sep_complex~=2),[1,3,2]).*(imag(sigmaaccum)+eps),[1,3]));
end
%% compact to save memory
for i=1:num_level
    T=asnarray2table(T,[],ysigmaaccum{i},[]);
end

if opt.compact.musigma
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_musigma,[]);
    end
end

%% aic
dfuni=df;
df=dyreal*df+1i*dycomplex*df;
% aic=2*dy*df+logsigma+sse+dy*n*log(2*pi);   %% multi variate
aic=2*df+logsigma+sse+(dyreal+1i*dycomplex)*n*log(2*pi);   %% multi variate *10.^2
% aic=real(aic)+imag(aic);   %% multi variate
%% aicc
penalty = n*(1+2*(real(df)+1)./(n-real(df)-2))+1i*n*(1+2*(imag(df)+1)./(n-imag(df)-2));
aicc = penalty + logsigma + sse +(dyreal+1i*dycomplex)*n*log(2*pi);

%% bic
% df=10.^2*df;
bic=log(n).*df+logsigma+sse+(dyreal+1i*dycomplex)*n*log(2*pi);   %% multi variate
% bic=log(n).*df+sigma+sse+(dyreal+1i*dycomplex)*n*log(2*pi); 
% bic=log(n*dy*df)+logsigma+sse;   %% multi variate
% bic=real(bic)+imag(bic);
%% bic
gamma=0.8;kappa=1/(2-2*gamma);
% ebic=log(n).*df*2*gamma*log(n)+logsigma+sse+(dyreal+1i*dycomplex)*n*log(2*pi);
ebic=(2*gamma*log(n.^kappa)+log(n)).*df+logsigma+sse+(dyreal+1i*dycomplex)*n*log(2*pi); 

%% gcv
% gcv=sigma./((1./n.*(n-df*dy)).^2);
% gcv=real(gcv)+imag(gcv);
%%
crit.gamma=gamma;
crit.kappa=kappa;
crit.n=n;
crit.dyreal=dyreal;
crit.dycomplex=dycomplex;

crit.dfuni=dfuni;
crit.df=df;
crit.sigma=sigma;
crit.logsigma=logsigma;
crit.sse=sse;
crit.aic=aic;
crit.aicc=aicc;
crit.bic=bic;
crit.ebic=ebic;


% crit.gcv=gcv;

end
% scatter3(A(:,1),A(:,2),T.(tregs{1}.opt.resp_sigma{1})(:,1))