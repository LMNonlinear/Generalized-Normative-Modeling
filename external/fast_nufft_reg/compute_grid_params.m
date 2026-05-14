function [Msp, Mr, tau,ratio]=compute_grid_params(N,accuracy)
%   
%     if nargin < 2 || isempty(eps)
%         eps = 1e-15;
%     end
    if nargin < 2 || isempty(accuracy)
        accuracy=8;
    end
%     Choose Msp & tau from eps following Dutt & Rokhlin (1993)
%     if eps <= 1e-33 || eps >= 1e-1
%         error(' must satisfy 1e-33 < eps < 1e-1.')
%     end
%     if eps > 1e-11 
%         ratio = 2;
%     else
%         ratio = 3;
%     end
%     Msp = round(-log(eps) / (pi * (ratio - 1) / (ratio - 0.5)) + 0.5)*ones(length(M),1);
%     Mr = max(ratio * M, 2 * Msp);
%     lambda = Msp / (ratio * (ratio - 0.5));
%     tau = pi * lambda./ M.^ 2;




ratio = 2;
% 
Msp=accuracy*ones(1,length(N));
tau = (pi*Msp./(N.*N*ratio*(ratio-.5)));
Mr = ratio*N;


end