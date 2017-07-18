%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PCT_AIF simulates the Aeterial Input Function
%
%   Ruogu Fang Revised 09/29/2011
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  AIF = PCT_AIF(T,T0,A,B,C);
%
%   PRE:
%       t       - time stems of the aif signal [T x 1]
%       t0      - delay of aif in seconds [Scalar]
%       a,b,c   - AIF simulation parameters (default a=1,b=3,c=1.5 s for
%                 adults, a=1,b=3.26,c=1.02 s for children)
%
%   POST:
%       AIF     - Arterial input function [T x 1]
%
% PCT_AIF simulates the arterial input function using 
%   AIF(t)  =   0                           t<=t0
%           =   a(t-t0)^b*e^(-(t-t0)/c)     t>t0
%
%   Example:
% t = 0 : 60;
% t0 = 0;
% a = 1;
% b = 3;
% c = 1.5;
% AIF = pct_aif(t,t0,a,b,c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AIF = pct_aif(t, t0, a, b, c)

% analyze input parameters
if nargin < 2
    t0 = 0;
end
if nargin < 3
    a = 1;
end
if nargin < 4
    b = 3;
end
if nargin < 5
    c = 1.5;
end

AIF = a * (t - t0).^b .* exp(-(t - t0) ./ c);
AIF(1 : t0 + 1) = 0;

end