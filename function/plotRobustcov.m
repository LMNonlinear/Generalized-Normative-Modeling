function plotRobustcov(mappedX,T,mah,outInd,group_base,tag,para,scale)
if nargin<7 || isempty(para)
    para={[];[]};
end
if nargin<8 ||isempty(scale)
    scale='linear';
end

num_group_base=length(group_base);
k=1;
d_classical = pdist2(mappedX, mean(mappedX),'mahal');
if strcmp(scale,'log')
    d_classical=log(d_classical);
    mah=log(mah);
end
p = size(mappedX,2);
chi2quantile = sqrt(chi2inv(0.975,p));

for igroup=1:num_group_base
    
    group=group_base{igroup};
    
    f=figure(600+k);
    clf
    
    % subplot(122)
    if ~exist('group','var') || isempty (group)
        scatter(d_classical,mah,30);
    else
        g=set_gscatter_para(T,group,2);
        gscatter(d_classical,mah,g.group,g.color,g.shape,15);%,
    end
    lgd=legend();lgdname={lgd.String{1:end}};
    if  exist('outInd','var') && ~isempty (outInd)
        hold on
        plot(d_classical(outInd), mah(outInd), 'r+')
        line([chi2quantile, chi2quantile], [0, max(mah)+0.3], 'color', 'r')
        line([min(d_classical)-0.3, max(d_classical)+0.5], [chi2quantile, chi2quantile], 'color', 'r')
        hold off
        xlim([0 chi2quantile+0.3])
        lgdname={lgd.String{1:end-2}};
        lgdname{end}='outliers';
    end
    xlabel('Mahalanobis Distance')
    ylabel('Robust Distance')
    
    
    
    legend(lgdname);
    
    figname=[num2str(igroup),'DD Plot,FMCD method',para{1},num2str(para{2})];
    title(figname);
    saveinfig(f,figname,tag)
    k=k+1;
end
end

% for kk=1:num_label*2
% saveinfig(f{kk},figname{kk},tag)
% end