function x=epan_kernel(x,h)
% x: support grid x (Nx1*Nx2*...Nxd*d)
%       or vector x (Tx*d)
% h: a scaller h, or a vector h

[dh,dx]=size(h);
d=ndims(x);
if d>dx %is ndgrid
    h=permute(h,[dx+2:-1:1]);
else
    h=permute(h,[dx+1:-1:1]);
end
% x=x./reshape(h,[ones(1,d-1),dx,dh]);
x=x./h;
x=squeeze(max(0,3/4*(1-sum((x.^2),d))));
end