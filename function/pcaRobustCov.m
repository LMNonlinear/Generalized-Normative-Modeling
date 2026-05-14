function [T,idx_outlier]=pcaRobustCov(T,resp,cluster_dim,group_base,tag,calc_robustcov,label_outluers,subject_col)
if nargin<6|| isempty(calc_robustcov)
    calc_robustcov=true;
end
if nargin<7|| isempty(label_outluers)
    label_outluers=true;
end
if nargin<8 || isempty(subject_col)
    subject_col='name';  % backward-compatible default
end
rng(1)

if nargin<3 || isempty(cluster_dim)
    cluster_dim=2;
end

Xtsne=extractSubsTable(T,resp,[],[],[],[],subject_col);
[mappedX]= tsne(Xtsne,'NumDimensions',cluster_dim,'Exaggeration',2);%,'Distance','mahalanobis'Distance','mahalanobis''NumPCAComponents',100
if calc_robustcov
    [~, ~, mah, outInd] =robustcov(mappedX);%,'OutlierFraction',0.05,'Method','ogk'
    % [mappedX]=tsne_drtoolbox(Xtsne,[],cluster_dim);
    if ~label_outluers
        outInd=[];
    end
    plotTsne(mappedX,T,group_base,tag,outInd);
    plotRobustcov(mappedX,T,mah,outInd,group_base,tag);
    if label_outluers
        t=writeOulter(T.(subject_col),outInd,'outlier');
        T.outlier=t.outlier;
        idx_outlier=T.outlier;
        T=T(~T.outlier,:);
    end
else
    idx_outlier=[];
    plotTsne(mappedX,T,group_base,tag,[]);
    
end






%% pre select data
% idxinter=true(length(T.name),1);
% idxinter(idouter)=false;
% T.outlier(idxinter)='N';
% T.outlier(~idxinter)='Y';
% T=T(idxinter,:);
% save(['./result/RobustCovTsne',tag])
%% add label to MultiSataInfo.csv
%     tablePath='"H:\PROCESSED_DATA\QMEEG\RiemanMultiData4Gam\MultiDataInfo_uniInfo.csv"';
%     M=readtable(tablePath);
%     subtableM=get_subtable(M,{'dataset','accept'},{'name'});
%     index_inter=find_char(M.name,subtableM.name);
%     M.outlier(index_inter)=repmat({'inter'},length(index_inter),1);
%     index_outer=find_char(M.name,con_name(Outfmcd));
%     M.outlier(index_outer)=repmat({'outer'},length(index_outer),1);
%     writetable(M,tablePath);

end











