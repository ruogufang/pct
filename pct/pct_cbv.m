function [ cbvmap ] = pct_cbv(rmap, rho, mask)
%PCT_CBV Calculates a CBV map from the residue functions.
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE: CBVMAP = PCT_CBV(CTMAP, RHO)
%
%   PRE:
%       RHOMAP  - A map of residue functions [T x Y x X]
%       RHO     - Brain tissue density in mL/g. [Scalar]
%   
%   POST:
%       CBVMAP  - A map of regional cerberal blood volume in 
%                 mL/100g [Y x X]
%
%   This function takes in a map of residue functions for each pixel
%   and outputs a CBV map. The residue function is the r(t) function
%   in the Dilutor-Indicator model. The CBV is calculated as the area under
%   the curve corrected for average brain tissue density.
%

if nargin<2
    rho = 1;
end
if nargin < 3
    dim = size(rmap);
    mask = ones(dim(2:end));
end

scaling_factor = 100/(rho);
cbvmap = scaling_factor * squeeze(sum(rmap,1));
cbvmap(cbvmap<0)=0;
cbvmap(~mask) = 0;

end

