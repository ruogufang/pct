%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PCT_TEC generates a map of tissue Time Enhancement Curves (TEC)
%
%   Ruogu Fang 09/29/2011
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  TEC = PCT_TEC(AIF, R, CBF)
%
%   PRE:
%       AIF     - Arterial input function [T x 1]
%       R       - A map of impulse residue functions [T x X x Y x Z]
%       CBF     - A map of cerebral blood flow [X x Y x Z] (optional)
%
%   POST:
%       C       - A map of time enhancement curves [T x X x Y x Z]
%
%   This function generates the impulse residue function as
%   C = CBF * conv(AIF,R)
%
% Example:
% t = (0 : 60)';
% Ca = aif(t,0,1,3,1.5);
% R = IRF(t, 10, 1.2);
% BF = 20;
% C = TEC(AIF, R, BF);
%

function TEC = pct_tec(AIF, RIF, CBF)

[T,X,Y,Z] = size(RIF);

if nargin < 3
    CBF = ones(X,Y,Z);
end

%Create A, the convolution matrix of the AIF
A = tril(toeplitz(AIF));

TEC = A*reshape(RIF,T,[]).*repmat(reshape(CBF,1,[]),[T 1]);
TEC = reshape(TEC,[T X Y Z]);

% TEC = zeros(T,X,Y,Z);
% for i = 1 : X
%     for j = 1 : Y
%         for  k = 1 : Z
%             c = CBF(i,j,k) * conv(AIF, RIF(:,i,j,k));
%             TEC(:,i,j,k) = c(1:T);
%         end
%     end
% end

end