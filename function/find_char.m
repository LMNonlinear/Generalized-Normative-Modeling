function index=find_char(char_cell1,char_cell2,opt,ignorecase)
% Authors:
% - Min Li
% - Ying Wang
% https://github.com/LMNonlinear
% https://github.com/rigelfalcon
% Date: May 22, 2021
if nargin<3
    opt=true;
end
if nargin<4
    ignorecase=false;
end
if ischar(char_cell2)
    char_cell2={char_cell2};
end

if ignorecase
    cmp=@strcmpi;
else
    cmp=@strcmp;
end
if opt
    try
        index=cellfun(@(x) find(cmp(char_cell1(:),x)),char_cell2(:),'UniformOutput',true);
    catch
        index=cellfun(@(x) find(cmp(char_cell1(:),x)),char_cell2(:),'UniformOutput',false);
%         nonindex=cellfun(@(x) isempty(x),index(:),'UniformOutput',true);
%         index(nonindex)={0};
        index=cell2mat(index);
    end
else
    index=cellfun(@(x) find(cmp(char_cell1(:),x)),char_cell2(:),'UniformOutput',opt);
end
end