function ismono=check_monotonicity(X,dim_check,dim_trial,isstrongly)
% dx: number of data dimension, the rest is the number of trial
% dim: dims to check

% example:
% assert(prod(check_monotonicity(repmat([1:100].',[1,10]),1,2)==true(10,1)))
% assert(check_monotonicity(repmat([1:100],[10,1]),2)==true)
% assert(check_monotonicity([ones(10,1),repmat([1:100],[10,1])],2,[],true)==false)
% assert(check_monotonicity([ones(10,1),repmat([1:100],[10,1])],2,[],false)==true)

sz=size(X);
if nargin<2 || isempty(dim_check)
    dim_check=length(sz);
end
if nargin<3 || isempty(dim_trial)
    dim_trial=length(sz)+1;
end 
if nargin<4 || isempty(isstrongly)
    isstrongly=false;
end 
if isstrongly
    comp=@(x,y) x>y;
else
    comp=@(x,y) x>=y;
end
dim_check=dim_check(:)';
dim_trial=dim_trial(:)';


dim_notcheck=setdiff(1:length(sz),[dim_check(:);dim_trial(:)].');
ismono=true([size(X,dim_trial),1]);
num_dim_check=length(dim_check);

for idim=1:num_dim_check
    idim_rest_check=setdiff(1:length(dim_check),idim);
    temp_ismono=sum( comp(diff(X,[],dim_check(idim)),0) ,[dim_check,dim_notcheck])==...
        prod(sz(dim_check(idim_rest_check)))*prod(sz(dim_notcheck))*(sz(dim_check(idim))-1);
    ismono=ismono & temp_ismono(:); 
end


%         sz(end+1)=1;
%     dim=1:length(sz);
% sz=size(X);
% % if nargin<3 || isempty(dx)
% %     dx=length(sz)-1;
% % end
% % if nargin<2 || isempty(dim)
% %     dim=1:length(sz)-1;
% % end 
% if nargin<3 || isempty(dx)
%     dx=length(sz);
%     sz(end+1)=1;
% end
% if nargin<2 || isempty(dim)
%     dim=1:length(sz);
% end 
% 
% % idx=1:length(sz);
% idim_notcheck=true(length(sz),1);
% idim_notcheck(dim)=false;
% ismono=true([sz(idim_notcheck),1]);
% for idim=1:length(dim)
%     idx=true(length(sz),1);
%     idx(dim(idim))=false;
%     idx(dim(idim))=false;
%     ismono=ismono & sum( (diff(X,idim)>0) ,[1:dx])==prod(sz(idx))*(sz(~idx)-1); 
% end
% %             ihbad(ihgood)=ihbad(ihgood)|...
% %                 reshape(~(...
% %                  sum( (diff(ytemp,idim)>0) ,[1:regs.dx])==...
% %                 (regs.N)^(regs.dx-1)*(regs.N-1)...
% %                 ),[],1); 