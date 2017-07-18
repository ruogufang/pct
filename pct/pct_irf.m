function IRF = pct_irf(t, alpha, beta)
%IRF generates a map of Impulse Residue Functions (IRF)
%
%   Ruogu Fang 09/29/2011
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  IRF = PCT_IRF(T, ALPHA, BETA)
%
%   PRE:
%       T       - A time series stems [T x 1]
%       ALPHA   - A map of parameters [Y x X]
%       BETA    - A map of sparameter [Y x X]
%
%   POST:
%       IRF     - A map of impulse residue functions [T x Y x X]
%
%   This function generates the impulse residue function as
%
%   R(t) = 1 - integal_0^t h(tau)d_tau = 1 - gamcdf(t,alpha,beta)
%   where h(t;alpha,beta) =
%           1/(beta^alpha*Gamma(alpha))*t^(alpha-1)*e(-t/beta)
%
%   Example:
%   MTT = 4; % sec
%   alpha = ones(40) * 10;
%   beta = MTT / alpha;
%   t = 0 : 60;
%   R = IRF(t, alpha, beta)
%

[X,Y] = size(alpha);
T = length(t);
IRF = 1 - gamcdf(repmat(t,[1 X Y]), permute(repmat(alpha,[1 1 T]),[3 1 2]), permute(repmat(beta,[1 1 T]),[3 1 2]));

end
