%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [R] = cTSVD(Ca, C, lambda, m)
%
% Circulant Truancated Singular Value Decomposition (cTSVD)
%
% Input:
%   Ca  - Concentration at artery input (AIF) [Tx1]
%   C   - Tissue Time Enhancement Curve (TEC) [TxXxYxZ]
%   lambda  - parameter for truncation (% of the signal) e.g. lambda = 0.06
%   m   -   extend the Ca matrix m times for block circulant
%
% Output:
%   R   - Residue function * CBF estimated by cTSVD algorithm
%   [TxXxYxZ]
% 
% Example:
% t = (0 : 60)';
% Ca = aif(t,0,1,3,1.5);
% R = irf(t, 10, 1.2);
% BF = 20;
% C = tec(Ca, R, BF);
% X = 10; Y = 10; Z = 10;
% Cv = repmat(C,[1 X Y Z]);
% lambda = 0.06;
% [R_cTSVD] = cTSVD(Ca, Cv, lambda);
%
% Ruogu Fang
% 9/29/2011
%
% Revised 6/3/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R] = cTSVD(Ca, C, lambda, m)

if nargin < 4 % no block circulant
    m = 1;
end

[T X Y Z] = size(C);
L = m * T;

cAIF = pct_bcaif(Ca,m); % Block-circulant AIF Matrix

[U,S,V] = svd(cAIF);
invS = inv(S);
idx = diag(S) < max(diag(S)) * lambda;
d = diag(invS);
d(idx) = 0;
tinvS = diag(d);
invAIF = V * tinvS * U';

cC = cat(1, C, zeros(L-T,X,Y,Z));
R = tprod(invAIF,[1 -1],cC,[-1 2 3 4]);
R = R(1:T,:,:,:);

