function [ outmap ] = pct_timefilter(inmap)
%PCT_TIMEFILTER Filters each pixel in time in an image sequence
%   
%   USAGE:  OUTMAP = PCT_TIMEFILTER(INMAP);
%
%   PRE:
%       INMAP   - A CT / image sequence [T x Y x X]
%
%   POST:
%       OUTMAP  - A filtered image sequence [T x Y x X]
%
%   Filters each pixel individually by convolving it with a Gaussian
%   FIR filter.
%   NOTE: Not yet completed. Do not use without reading code first.
%
%   Kolbeinn Karlsson 06/12/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

% Get image dimensions
[time height width] = size(inmap);

% Pre-allocate the output variable
outmap = zeros([time height width]);

% Create the Gaussian kernell
h = gaussfir(0.3,2,4);
tail = 8;

for i = 1:height
    for j = 1:width
        R = inmap(:,i,j);
        R = conv(R,h);
        outmap(:,i,j) = R(tail+1:end-tail);
    end
end

end

