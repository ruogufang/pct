function [ outmap ] = pct_downsample(inmap, k)
%PCT_DOWNSAMPLE Downsamples an image / map by k
%
%   USAGE:  OUTMAP = PCT_DOWNSAMPLE(INMAP, K);
%
%   PRE:
%       INMAP   - The input image [M x N]
%       K       - The downsample factor [Scalar]
%
%   POST: 
%       OUTMAP  - The downsampled image [M/K x N/K]
%
%   Downsamples an image by k. If M/K or N/K are not integers, they will be
%   rounded up to the nearest integer.
%
%   Kolbeinn Karlsson 06/14/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

inmap = downsample(inmap,k);
inmap = downsample(inmap',k);
outmap = inmap';


end

