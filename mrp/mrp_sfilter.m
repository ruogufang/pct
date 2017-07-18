function [ outmap ] = mrp_sfilter(inmap, fsize)
%MRP_SFILTER Spatially filters each frame of a CT series.
%
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%
%   USAGE:  OUTMAP = MRP_SFILTER(INMAP, FSIZE);
%
%   PRE:
%       INMAP   - A MRP Data [T x X x Y x Z]
%       FSIZE   - Size of the filter [Scalar]
%
%   POST:
%       OUTMAP  - A filtered MRP data [T x X x Y x Z]
%
%   This function spatially filters each frame individually using the
%   Gaussian filter kernell of size FSIZE. 
%

%Get the image dimensions
[T,X,Y,Z] = size(inmap);

%Create the Gaussian filter kernell
h = fspecial('gaussian',fsize);

%Pre-allocate the ouput variable
outmap = zeros(size(inmap));

%Filter each frame individually
for t = 1:T
    outmap(t,:,:,:) = imfilter(squeeze(inmap(t,:,:,:)),h);
end

end

