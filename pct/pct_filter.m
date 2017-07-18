function [ outmap ] = pct_filter(ctmap, fsize, stddev)
%PCT_FILTER Spatially filters each frame of a CT series.
%
%   Kolbeinn Karlsson 06/07/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUTMAP = PCT_FILTER(CTMAP, FSIZE);
%
%   PRE:
%       CTMAP   - A CT Map in HU [T x Y x X]
%       FSIZE   - Size of the filter [Scalar]
%       STDDEV  - The standard deviation of the filter [Scalar]
%
%   POST:
%       OUTMAP  - A filtered CT map in HU [T x Y x X]
%
%   This function spatially filters each frame individually using the
%   Gaussian filter kernell of size FSIZE. 
%

if nargin < 3
    stddev = 1;
end

%Get the image dimensions
[time,height,width] = size(ctmap);

%Create the Gaussian filter kernell
h = fspecial('gaussian',fsize,stddev);

%Pre-allocate the ouput variable
outmap = zeros(size(ctmap));

%Filter each frame individually
for i = 1:time
    outmap(i,:,:) = imfilter(squeeze(ctmap(i,:,:)),h);
end

end

