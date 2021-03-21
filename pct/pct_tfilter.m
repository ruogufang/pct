function [ outmap ] = pct_tfilter(ctmap, tsize)
%PCT_TFILTER Temporally filters each pixel of a CT series.
%
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUTMAP = PCT_TFILTER(CTMAP, FSIZE);
%
%   PRE:
%       CTMAP   - A CT Map in HU [T x Y x X]
%       TSIZE   - Size of the moving average filter [Scalar]
%
%   POST:
%       OUTMAP  - A filtered CT map in HU [T x Y x X]
%
%   This function temporally filters each pixel individually using the
%   moving average filter kernell of size TSIZE. 
%

%Get the image dimensions
[time height width] = size(ctmap);

%Create the Gaussian filter kernell
h = ones(1,tsize)/tsize;

%Pre-allocate the ouput variable
outmap = zeros(size(ctmap));

%Filter each frame individually
for i = 1:height
    for j = 1:width
        outmap(:,i,j) = conv(ctmap(:,i,j),h,'same');
    end
end

end

