function [ outmap ] = pct_gaussfilter(inmap, l, stddev)
%PCT_GAUSSFILTER Filters each pixel in time in an image sequence
%   
%   USAGE:  OUTMAP = PCT_GAUSSFILTER(INMAP);
%
%   PRE:
%       INMAP   - A CT / image sequence [T x Y x X]
%       L       - Length of the filter [Scalar]
%       STDDEV  - The standard deviation of the filter [Scalar]
%
%   POST:
%       OUTMAP  - A filtered image sequence [T x Y x X]
%
%   Filters each pixel individually by convolving it with a windowed Gaussian
%   function.
%
%   Kolbeinn Karlsson 06/12/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

% Get image dimensions
[time height width] = size(inmap);

% Pre-allocate the output variable
outmap = zeros([time height width]);

if nargin >= 3
    % Create the Gaussian kernell
    alpha = 1/stddev; 
    h = gausswin(l,alpha);
else
    h = gausswin(l);
end
h = h/sum(h); %normalize the filter
tail = floor(l/2);

for i = 1:height
    for j = 1:width
        R = inmap(:,i,j);
        R = conv(R,h);
        outmap(:,i,j) = R(1+tail:end-tail);
    end
end

end

