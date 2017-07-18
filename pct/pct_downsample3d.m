function [ outmap ] = pct_downsample3d(inmap, k)
%PCT_DOWNSAMPLE3D Spatially downsamples an image sequence by k
%
%   USAGE:  OUTMAP = PCT_DOWNSAMPLE3D(INMAP, K);
%
%   PRE:
%       INMAP   - An image sequence [T x M x N]
%       K       - The downscale factor [Scalar[
%
%   POST:
%       OUTMAP  - The downsampled sequence [T x M/K x N/K]
%
%   If M/K or N/K are not integers, they are rounded up.
%
%   Kolbeinn Karlsson 06/14/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Get the input dimensions
[T M N] = size(inmap);

%Pre-allocate the output variable
outmap = zeros(T, ceil(M/k), ceil(N/k));

%Downsample each frame individually
for i = 1:T
    outmap(i,:,:) = pct_downsample(squeeze(inmap(i,:,:)),k);
end


end

