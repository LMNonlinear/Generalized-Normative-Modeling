function name=add_tag(name,tag,pos)
if nargin<3||isempty(pos)
    pos='after';
end
% regexprep(name, tag, 'Wang','ignorecase')
% TF = contains(str,tag)
if ischar(tag)
    tag={tag};
end
for i=1:length(tag)
    if ~contains(name,tag{i},'IgnoreCase',true)
        if ischar(pos) && strcmpi(pos,'after')
            name=strcat(name,'_',tag{i});
        elseif ischar(pos) && strcmpi(pos,'before')
            name=strcat(tag{i},'_',name);
        else
            error('not suport ')
        end
    end
end
end