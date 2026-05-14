function [Hlist]=getHlist(M,range,space)

dh=length(M);
if nargin<2 || isempty(range)
    range=repmat([0.01,1],dh);
end
if nargin<3 || isempty(space)
    space=@linspace;
end

%  if isscalar(H)
Hlist=multispace(range(:,1),range(:,2),M,space)';
Hlist=get_ndgrid(Hlist,'list');
% Hlist=reshape(Hlist,[],dh);
%  else
%         HList=H;
% end
end