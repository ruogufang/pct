function [bins vals] = pct_showhist(img, nbins, range)
%PCT_SHOWHIST Computes the histogram of an image and displays it
%
%   USAGE:  [BINS VALS] = PCT_SHOWHIST(IMG, NBINS, RANGE);
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
%   This function uses pct_hist to compute the histogram and immediately
%   displays the results. See pct_hist for more details.
%
%   Kolbeinn Karlsson 06/18/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if nargin > 2
    [bins vals] = pct_hist(img, nbins, range);
else
    [bins vals] = pct_hist(img, nbins);
end

stem(bins,vals,'MarkerSize',1);

end

