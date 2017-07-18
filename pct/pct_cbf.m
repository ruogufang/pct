function [ cbfmap ] = pct_cbf(rmap, rho, mask)
%PCT_CBF Calculates a Cerebral Blood Flow map from a map of residue
%functions
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  CBFMAP = PCT_CBF(RMAP, RHO)
%
%   PRE:
%       rmap    - A map of residue functions [T x Y x X]
%       rho     - Brain tissue density in g/mL [Scalar]
%
%   POST:
%       cbfmap  - A map of regional cerebral blood flow in 
%       mL/100g tissue/min [Y x X]
%
%   This function calculates the cerebral blood flow from a map of residue
%   functions using the Indicator-Dilutor model. The CBF is calculated as
%   the maximum of the residue function corrected for average brain tissue
%   density.

if nargin<2
    rho = 1;
end

if nargin<3
    dim=size(rmap);
    mask = ones(dim(2:end));
end

%Convert to mL/100g/min
scaling_factor = 60*100/rho;

%Calculate the CBF map
cbfmap = scaling_factor * squeeze(max(rmap));

cbfmap(cbfmap<0)=0;
cbfmap(~mask) = 0;

end

