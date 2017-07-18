function [ mttmap ] = pct_mtt(rmap, mask, first, last)
%PCT_MTT Calculates an MTT map from a map of residue functions
%
%   USAGE: MTTMAP = PCT_MTT(RMAP);
%
%   PRE:
%       RMAP    - A map of residue functions [T x Y x X]
%       FIRST   - The first frame to include in the calculations [Scalar]
%                 (optional)
%       LAST    - The last frame to include in the calculations [Scalar]
%                 (optional)
%
%   POST:
%       MTTMAP  - A map of Mean Transit Time in seconds [Y x X]
%
%   This function calculates the mean transit time (MTT) from a map of
%   residue functions from the Indicator-Dilutor model. The MTT is
%   calculated using the Central Volume Theorem, i.e. as the ratio of the
%   CBV and the CBF. If FIRST and LAST are not supplied, then the MTT is
%   calculated from the whole sequence.
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Truncate the data if necessary
if nargin < 2
    dim = size(rmap);
    mask = ones(dim(2:end));

end
if nargin == 4
    rmap = rmap(first:last,:,:);
end

mttmap = squeeze(sum(rmap,1)) ./ (squeeze(max(rmap))+eps);
mttmap(mttmap<0)=0;
mttmap(~mask) = 0;

end

