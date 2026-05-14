function [gcv]=nureg_gcv(regs,y,yhat,fpp,metric,dstd)
if nargin<6 || isempty(dstd)
    dstd=regs.opt.y_gcv_dstd;
end

if nargin<5 || isempty(metric)
    metric='res';
end
[n,dy]=size(y);
dh=regs.dh;
ihbad=regs.ihbad;


yhat=permute(yhat,[1,3,2]);
switch metric
    case {'res'}
        rss=reshape(sum((y-yhat).^2),dy,regs.dh).';
        mse=reshape((1/regs.Ty)*sum((y-yhat).^2),dy,regs.dh).';
        gcv_m=mse.*regs.pdof_inv_m;
        gcv_sd=mse.*regs.pdof_inv_sd;
        %         aic=2*regs.df_m+regs.Ty*log(mse);
    case {'monotone'}
        ytemp=fpp(regs.xlist);
        ihnomono=false(regs.dh,dy);
        for iy=1:dy
            ihnomono(:,iy)=~check_monotonicity(ytemp(regs.ndcolon{:},(1:regs.dh)+regs.dh*(iy-1)),dy-iy+1,regs.dx+1,false);% only for xi
        end
        
        if sum(sum(ihnomono)==size(ihnomono,1))
            error('no bandwidth make monotone')
        end
        yhat(:,ihnomono')=inf;
        mse=reshape((1/regs.Ty)*sum((y-yhat).^2),dy,regs.dh).';
        %         gcv_m=rss.*regs.pdof_inv_m;
        gcv_m=abs(mse.*regs.pdof_inv_m);
        gcv_sd=mse.*regs.pdof_inv_sd;
    otherwise
        error('only [res],[monotone] support yet')%skewness
end
gcv_up=gcv_m+gcv_sd;
gcv_lo=gcv_m-gcv_sd;
% aic=2*regs.df_m+n*log(mse);   %% 2k+nln(RSS) Here RSS is sigma.^2  (sum(yi-ybar).^2)./n
% bic=log(n)*regs.df_m+n*log(mse);          %% nln(RSS)+k*ln(n)
if nargin<4 || isempty(fpp)
    fpp=[];
end
% ihgood=cell(dy,1);
% for iy=1:size(ihbad,2)
%     ihgood{iy}=find(~ihbad(:,iy));
% end
ihgood=find(~ihbad);

%% minimu error
[~,idmingood]=min(gcv_m,[],1);
%% largest in the std, for smooth
id1segood=zeros(1,dy);
for iy=1:dy
    semin=gcv_m(idmingood(iy),iy)+dstd*gcv_sd(idmingood(iy),iy).';
    [id1sey]=find(gcv_m(:,iy)==max(gcv_m(gcv_m(:,iy)<=semin,iy)));
    id1segood(iy)=id1sey(end);
end
if ~isempty(fpp)
    fpp.Values=fpp.Values(regs.ndcolon{:},id1segood+dh*((1:dy)-1));
end

yhat_min=zeros(regs.Tx,dy);
yhat_1se=zeros(regs.Tx,dy);
for iy=1:dy
    yhat_min(:,iy)=squeeze(yhat(:,iy,idmingood(iy)));
    yhat_1se(:,iy)=squeeze(yhat(:,iy,id1segood(iy)));
end

% for iy=1:size(ihbad,2)
%     idmin(iy)=ihgood{iy}(idmin(iy)).';
%     id1se(iy)=ihgood{iy}(id1se(iy)).';
% end
idmin=idmingood;
id1se=id1segood;
for iy=1:dy
    idmin(iy)=ihgood(idmingood(iy)).';
    id1se(iy)=ihgood(id1segood(iy)).';
end
hmin=regs.h(idmin,:).';
h1se=regs.h(id1se,:).';
% aicmin=aic(idmin);
% aic1se=aic(id1se);
% bicmin=bic(idmin);
% bic1se=bic(id1se);
msemin=mse(idmingood);
mse1se=mse(id1segood);
rssmin=rss(idmingood);
rss1se=rss(id1segood);
pdof_m_min=reshape(regs.pdof_m(idmingood),1,[]);
pdof_m_1se=reshape(regs.pdof_m(id1segood),1,[]);
df_m_min=reshape(regs.df_m(idmingood),1,[]);
df_m_1se=reshape(regs.df_m(id1segood),1,[]);
for iy=1:dy
    disp(['[nw_nufft]: [',num2str(regs.dx),'d] regression bandwidth - h1se [',num2str(h1se(:,iy).'),'] is [',num2str(find(ismember(regs.h, h1se(:,iy).','row'))),'th] of h'])
end


%%
switch metric
    case {'res'}
    case {'monotone'}
        gcv.ihnomono=ihnomono;
    otherwise
        error('only [res],[monotone] support yet')%skewness
end

%%
gcv.metric=metric;
gcv.yhat_min=yhat_min;
gcv.yhat_1se=yhat_1se;
gcv.h1se=h1se;
gcv.hmin=hmin;
gcv.gcv_m=gcv_m;
gcv.gcv_sd=gcv_sd;
gcv.gcv_up=gcv_up;
gcv.gcv_lo=gcv_lo;
gcv.idmin=idmin;
gcv.id1se=id1se;
gcv.fpp=fpp;
% gcv.aicmin=aicmin;
% gcv.aic1se=aic1se;
% gcv.bicmin=bicmin;
% gcv.bic1se=bic1se;
gcv.rss=rss;
gcv.msemin=msemin;
gcv.mse1se=mse1se;
gcv.rssmin=rssmin;
gcv.rss1se=rss1se;
gcv.pdof_m_min=pdof_m_min;
gcv.pdof_m_1se=pdof_m_1se;
gcv.df_m_min=df_m_min;
gcv.df_m_1se=df_m_1se;

end