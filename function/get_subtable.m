function [datatable,row_indx_orig]=get_subtable(datatable,condition,items)
num_condition=numel(condition);
row_indx_orig=1:size(datatable,1);
for i=1:2:num_condition
    item=condition{i};
    property=condition{i+1};
    if ischar(property)
        row_indx=strcmp(datatable.(item),property);
    elseif iscell(property)
        row_indx=find_char(datatable.(item),property,true);
    elseif istable(property)
        row_indx=find_char(datatable.(item),table2cell(property),true);
    elseif isnumeric(property)
%         error('not support numeric property yet')
        row_indx=find(ismember(datatable.(item),property));
    end
    datatable=datatable(row_indx,:);
    row_indx_orig=row_indx_orig(row_indx);
end
if nargin>2 && ~isempty(items)
    % num_items=numel(items);
    col_idx=find_char(datatable.Properties.VariableNames,items);
    datatable=datatable(:,col_idx);
end
    

end




