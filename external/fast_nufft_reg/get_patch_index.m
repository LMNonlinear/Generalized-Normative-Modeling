function idx_patch=get_patch_index(q,M)
% q is a d*2 matrix represent the subpatch, q(:,1) is start index, q(:,2) is end index
% M is a d*1 vector denote the whole space with each direction of length
dx=size(q,1);
dM=length(M);


N=q(:,2)-q(:,1)+1;
idx_patch=multispace(q(:,1),q(:,2),N,[],'cell',true);
if dM>dx
    for i=1:dM-dx
        idx_patch(dx+i)={1:M(dx+i)};
    end
end
idx_patch=get_ndgrid(idx_patch,'cellList');

% idx_patch=cellfun(@(x) reshape(x,[],1),idx_patch,'UniformOutput',false);
idx_patch=sub2ind([M(:);1].',idx_patch{:});
end