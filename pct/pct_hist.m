function [bins vals] = pct_hist(img, nbins, range)
%PCT_HIST Calculates a histogram of a 2D image
%
%   USAGE:  [BINS VALS] = PCT_HIST(IMG, NBINS, RANGE);
%
%   PRE:
%       IMG         - A 2D matrix [M x N]
%       NBINS       - Number of bins [Scalar] (optional)
%       RANGE       - A two-element vector specifying the range of the
%                     histogram [2 x 1] (optional).
%
%   POST:
%       BINS        - An array of bin center values [1 x NBINS]
%       VALS        - The number of pixels in each bin [1 x NBINS]
%
%   Computes the histogram of a 2D image. If NBINS is not supplied,
%   256 bins will be computed. If RANGE is not specified, then a
%   histogram for the full range of the image is computed.
%
%   Kolbeinn Karlsson 06/18/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%

%Apply standard arguments if less than 3 arguments are passed
if nargin <= 2
    range = [min(min(img)) max(max(img))];
end
if nargin == 1
    nbins = 256;
end

%Get the thresholds
minval = range(1);
maxval = range(2);

%Get image dimensions
[M N] = size(img);

%Flatten image
img = reshape(img,1,M*N);

%Calculate the edges
if minval == maxval
    incr = 1/nbins;
    bins = (minval-0.5):incr:(maxval+0.5);
else
    incr = (maxval-minval)/nbins;
    bins = minval:incr:maxval;
end

%Get the historgram
vals = histc(img, bins);

end

