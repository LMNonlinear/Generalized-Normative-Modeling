function [grid,df]=nufftfreqs(N, df)
dx=length(N);
if nargin<2 || isempty(df)
    df=1*ones(1,dx);
end


% Compute the frequency range used in nufft for M frequency bins
grid=df .* linspacemulti(fix(-N/2),N - fix(N/2)-1,N).';
grid=get_ndgrid(grid);
grid=permute(grid,[dx+1:-1:1]);
grid=cell2mat(grid);
end