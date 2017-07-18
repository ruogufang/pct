function [ outmap ] = pct_mapresize(inmap, M, N)
%PCT_MAPRESIZE Resizes a map/image to [M x N] using bilinear interpolation
%
%   USAGE:  OUTMAP = PCT_MAPRESIZE(INMAP, M, N);
%
%   PRE:
%       INMAP   - A map/image [Y x X]
%       M       - The desired output height (nrows) [Scalar]
%       N       - The desired output width (ncols) [Scalar]
%
%   POST:
%       OUTMAP  - A rescaled version of INMAP [M x N]
%
%   Kolbeinn Karlsson 06/14/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%This is actually quite simple
outmap = imresize(inmap, [M N], 'bicubic');


end

