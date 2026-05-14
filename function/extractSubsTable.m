function [Y,idx_complex,X]=extractSubsTable(T,resp,reture_complex,d,isstd,factor,subject_col)
if iscell(resp)
    num_resp=length(resp);
elseif ischar(resp)
    num_resp=1;
end

if nargin<5 || isempty(isstd)
    isstd=false;
end

if nargin<4 || isempty(d)
    d=2;
end
if nargin<3 || isempty(reture_complex)
    reture_complex=false;
end
if nargin<6
    factor={};
end
if nargin<7 || isempty(subject_col)
    subject_col='name';  % backward-compatible default
end

% infer num_obs_per_subject from data (replaces hardcoded num_freq=47)
con_name=unique(T.(subject_col));
num_name=length(con_name);
num_obs_per_subject = size(get_subtable(T,{subject_col,con_name(1)},resp(:)),1);

idx_complex=find(isComplexTableVar(T,[resp(:)]));
num_complex=length(idx_complex);

if ~isempty(factor)
    X=zeros(num_name,num_obs_per_subject,length(factor));
else
    X=[];
end
Y=zeros(num_name,num_obs_per_subject,num_resp+num_complex);
for i=1:length(con_name)
    y= (table2array(get_subtable(T,{subject_col,con_name(i)},[resp(:)])));
    Y(i,:,1:num_resp)=real(y);
    Y(i,:,num_resp+1:num_resp+num_complex)=imag(y(:,idx_complex));
    if ~isempty(factor)
        X(i,:,:)=table2array(get_subtable(T,{subject_col,con_name(i)},factor));
    end
end

Y=reshape(Y,[],num_resp+num_complex);
if isstd
    sd=std(Y);
    sd(idx_complex)=sqrt(sd(idx_complex).^2+sd(num_resp+1:num_complex+num_resp).^2)/2;
    sd(num_resp+1:num_complex+num_resp)=sd(idx_complex);
    Y=Y./sd;
end
if reture_complex
    Y(:,idx_complex)=Y(:,idx_complex)+1i*Y(:,num_resp+1:num_complex+num_resp);
    Y=Y(:,1:num_resp);
end
if d==2
    Y=reshape(Y,num_name,[]);
else
    Y=reshape(Y,num_name,num_obs_per_subject,[]);
end




end


