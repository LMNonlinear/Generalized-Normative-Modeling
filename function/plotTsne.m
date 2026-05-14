function plotTsne(mappedX,T,group_base,tag,outIdx,para)
if nargin<5 || isempty(outIdx)
    outIdx=[];
end
    
if nargin<6 || isempty(para)
    para={[];[]};
end


num_group_base=length(group_base);
k=1;
for igroup=1:num_group_base
    f=figure(500+igroup);
    clf
    group=group_base{igroup};
    switch size(mappedX,2)
        case 2
            if ~exist('group','var') || isempty (group)
                scatter(mappedX(:,1),mappedX(:,2),30);
            else
                g=set_gscatter_para(T,group_base{igroup},size(mappedX,2));
                gscatter(mappedX(:,1),mappedX(:,2),g.group,g.color,g.shape,7);
            end
            if exist('outIdx','var') || isempty (outIdx)
                hold on
                scatter(mappedX(outIdx,1),mappedX(outIdx,2),30,'r+');
            end
        otherwise
            if ~exist('group','var') || isempty (group)
                scatter3(mappedX(:,1),mappedX(:,2),mappedX(:,3),30);
            else
                g=set_gscatter_para(T,group,size(mappedX,2));
                gscatter3(mappedX(:,1),mappedX(:,2),mappedX(:,3),g.group,g.color,g.shape,7);%,'auto'
            end
            if exist('outIdx','var') || isempty (outIdx)
                hold on
                scatter3(mappedX(outIdx,1),mappedX(outIdx,2),mappedX(outIdx,3),30,'r+');
            end
            lgd=legend();
            lgdname=lgd.String;
            lgdname{end}='outliers';
            legend(lgdname);
    end
    figname=[num2str(igroup),'t-sne',para{1},num2str(para{2})];
    title(figname)
    saveinfig(f,figname,tag)
    k=k+1;
end
end
% for kk=1:num_label*2
% saveinfig(f{kk},figname{kk},tag)
% end