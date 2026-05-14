function knots=scale_knots(knots,N,xlist)
% xlist should be uniform grid list in a cell
dx=length(N);
kmin=min(knots);
kmax=max(knots);
% kmin=permute(kmin,[3,1,2]);
% kmax=permute(kmax,[3,1,2]);

if nargin>2 && ~isempty(xlist)
    if dx>1 || iscell(xlist)
        hx=cellfun(@(x) x(2)-x(1),xlist);
    else
        hx=xlist(2)-xlist(1);
    end
%     hx=permute(hx,[3,1,2]);
%     N=permute(N,[3,1,2]);
    kmean=(kmin+kmax)/2;
    kmin=kmean(1)-(N/2).*hx;
    kmax=kmin+(N-1).*hx;
end



bw=kmax-kmin;
scale=(N-1)./bw;
shift=-N/2-kmin.*scale;
knots=scale.*knots + shift;
%Switch knot locations to [0,2*pi] convention (used by Greengard):
knots=mod(2*pi*knots./N,2*pi);%Makes NUFFT implementation
%more straightforward when notation is the same!






