function T=asnarray2table(T,cond,varname,arr)
% [icol,T]=get_varidx(T,opt.resp_mu,regs.dtype);
% [icol,T]=get_varidx(T,opt.resp_mu);
[~,irow]=get_subtable(T,cond);
% [T,idx]=get_varidx(T,varname,arr(1));

if ischar(varname)
    varname={varname};
end
num_var=length(varname);
if istable(arr)
    arr=table2array(arr);
end
dim=ndims(arr);
if size(arr,dim)~=num_var && num_var==1
    ndcolon(1:dim) = {':'};
elseif size(arr,dim)==num_var || isempty(arr) || isscalar(arr)
    ndcolon(1:dim-1) = {':'};
else
    error('number of variable is not equals to the data dimenssion')
end

for i=1:num_var
    if isscalar(arr)
        if find_char(T.Properties.VariableNames,varname{i})
            if size(arr,2)>length(T.(varname{i})(1,:))
                T.(varname{i})=[];
            end
%         else
%             T=initializeTableVar(T,varname{i},class(arr));
        end
        T.(varname{i})(irow,:)=arr;
    elseif isempty(arr)
        if find_char(T.Properties.VariableNames,varname{i})
            T.(varname{i})=[];
        end
    elseif size(arr,1)==1
        if find_char(T.Properties.VariableNames,varname{i})
            if size(arr,2)>length(T.(varname{i})(1,:))
                T.(varname{i})=[];
            end
        end
        T.(varname{i})(irow,:)=repmat(arr(ndcolon{:},i),numel(irow),1);
    else
        if find_char(T.Properties.VariableNames,varname{i})
            if size(arr,2)>length(T.(varname{i})(1,:))
                T.(varname{i})=[];
            end
        end
        T.(varname{i})(irow,:)=arr(ndcolon{:},i);
    end
end


end