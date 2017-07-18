function [ mttmap ] = pct_mtt(rmap)
%PCT_MTT Calculates a Mean Transit Time map from a map of residue
%functions
%
%   Ruogu Fang 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  MTTMAP = PCT_MTT(RMAP, RHO)
%
%   PRE:
%       rmap    - A map of residue functions [T x Y x X]
%
%   POST:
%       mttmap  - A map of regional mean transit time in 
%       sec [Y x X]
%
%   This function calculates the mean transit time from a map of residue
%   functions using the Indicator-Dilutor model. The MTT is calculated as
%   the first moment of the residue function corrected for average brain tissue
%   density.

[T Y X] = size(rmap);

%Calculate the MTT map
mttmap = squeeze(sum(repmat((1:T)',[1 Y X]).*rmap) ./ (sum(rmap)+eps));

end

