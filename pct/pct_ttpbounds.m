function [ ttpmap ] = pct_ttpbounds(inmap,dt,minttp,maxttp)
%PCT_TTP Calculates a Time-To-Peak map from a CT sequence within given bounds
%
%   USAGE:  TTPMAP = PCT_TTPBOUNDS(INMAP, DT, MINTTP, MAXTTP);
%
%   PRE:
%       INMAP   - A CT / Tracer concentration map sequence [T x Y x X]
%       DT      - Time interval between samples in seconds [Scalar]
%       MINTTP  - The minimum value of the TTP [Scalar]
%       MAXTTP  - The maximum value of the TTP [Scalar]
%
%   POST:
%       TTPMAP  - A Time-To-Peak map [Y x X]
%
%   This function calculates a TTP map from a CT / tracer concentration 
%   map. It estimates the TTP of each pixel with the maximum within the bounds
%   of MINTTP and MAXTTP.
%
%   Kolbeinn Karlsson 06/12/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University


%Get the inmap dimensions
[time height width] = size(inmap);

%Pre-allocate the output variable
ttpmap = zeros(height, width);

%Convert the bounds from time units to index units
lo = minttp/dt;
hi = maxttp/dt;

%Do the calc's
for i = 1:height
    for j = 1:width
        [V IND] = max(inmap(lo:hi,i,j));
        ttpmap(i,j) = IND+lo;
    end
end

%Convert to seconds
ttpmap = dt * ttpmap;

end

