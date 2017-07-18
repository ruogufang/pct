function [ mttmap ] = pct_mttmoments(rmap, dt, first, last)
%PCT_MTTMOMENTS Calculates an MTT map from a map of residue functions
%
%   USAGE: MTTMAP = PCT_MTTMOMENTS(RMAP, DT, FIRST, LAST);
%
%   PRE:
%       RMAP    - A map of residue functions [T x Y x X]
%       DT      - Time interval between samples in seconds [Scalar]
%       FIRST   - The first frame to include in the calculations [Scalar]
%                 (optional)
%       LAST    - The last frame to include in the calculations [Scalar]
%                 (optional)
%
%   POST:
%       MTTMAP  - A map of Mean Transit Time in seconds [Y x X]
%
%   This function calculates the mean transit time (MTT) from a map of
%   residue functions from the Indicator-Dilutor model. The MTT is computed as
%   the first moment of the residue function. If FIRST and LAST are not
%   supplied, then the MTT is calculated from the whole sequence.
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Truncate the data if necessary
if nargin >= 4
    rmap = rmap(first:last,:,:);
end

[T H W] = size(rmap);

mttmap = zeros(H,W);

t = dt*((1:T)-1);

for i = 1:H
    for j = 1:W
        mttmap(i,j) = sum(rmap(:,i,j) .* t(:))/sum(dt*rmap(:,i,j));
    end
end


end

