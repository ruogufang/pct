function [ ttpmap ] = pct_ttp(inmap,dt,mask)
%PCT_TTP Calculates a Time-To-Peak map from a CT sequence
%
%   USAGE:  TTPMAP = PCT_TTP(INMAP);
%
%   PRE:
%       INMAP   - A CT / Tracer concentration map sequence [T x Y x X]
%       DT      - Time interval between samples in seconds [Scalar]
%       MASK    - Mask of valid pixels [Y x X]
%
%   POST:
%       TTPMAP  - A Time-To-Peak map [Y x X]
%
%   This function calculates a TTP map from a CT / tracer concentration 
%   map. It basically shows you the time index of
%   the peak of each pixel.
%
%   Kolbeinn Karlsson 06/12/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Get the inmap dimensions
[time,height,width] = size(inmap);

%Pre-allocate the output variable
ttpmap = zeros(height, width);

if nargin < 3
    mask = ones(height,width);
end

%Do the calc's
for i = 1:height
    for j = 1:width
        [V,IND] = max(inmap(:,i,j));
        ttpmap(i,j) = IND;
    end
end

ttpmap = ttpmap-1;

%Convert to seconds
ttpmap = dt * ttpmap;

ttpmap(~mask) = 0;

end

