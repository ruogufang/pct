function [ out ] = mrp_subtract(in, frames)
%MRP_SUBTRACT generates subtraction images by subtracting a fixed precontrast
%images from the contrast concentration images
%   
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%
%   USAGE:  OUT = MRP_SUBTRACT(IN, FRAMES);
%
%   PRE:    
%       IN     - A series of CT maps [T x X x Y x Z]
%       FRAMES - No. of frames to average
%
%   POST:
%       OUT    - A series of contrast concentration maps [T x X x Y x Z]
%
%   This function calcualtes the average of each pixel in the first
%   FRAMES frames and subtracts it from every frame
%

%Calculate the average - the CT values before the contrast is introduced
avg = mean(in(1:frames,:,:,:),1);

%Get input dimensions
[T,X,Y,Z] = size(in);

%Pre-allocate output variable
out = zeros(size(in));

%Convert to contrast
for t = 1:T
    out(t,:,:,:) = in(t,:,:,:) - avg;
end

end

