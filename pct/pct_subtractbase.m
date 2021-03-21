function [ out ] = pct_subtractbase(ctmap, frames, k)
%PCT_SUBTRACTBASE Converts CT sequence to a contrast concentration sequence
%   
%   USAGE:  OUT = PCT_SUBTRACTBASE(IN, FRAMES, K);
%
%   PRE:    
%       CTMAP   - A CT sequence [T x Y x X]
%       FRAMES  - No. of frames to average [Scalar]
%       K       - The contrast conversion factor in g/mL/HU [Scalar]
%
%   POST:
%       OUT     - A contrast concentration map sequence [T x Y x X]
%
%   This function calcualtes the average of each pixel in the first
%   FRAMES frames and subtracts it from every frame.
%
%   Kolbeinn Karlsson 06/04/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%For backwards compatibility
if nargin < 3
    k = 1;
end

%Calculate the average - the base image
avg = mean(ctmap(1:frames,:,:),1);

%Get input dimensions
[t, Y, X] = size(ctmap);

%Pre-allocate output variable
out = zeros([t,Y,X]);

%Convert to contrast
for i = 1:t
    out(i,:,:) = k*(ctmap(i,:,:) - avg);
end

end

