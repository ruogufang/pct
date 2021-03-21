function [out,v] = pct_velim(CBV,mask,thresh)
%PCT_VELIM detects and remove the vascular pixels from the input image
%
%   USAGE:  K = PCT_VELIM(IN,MASK,THRESH)
%
%   INPUT:
%       CBV     - Input CBV map [Y x X]
%       MASK    - A logical [Y x X] mask. The computation will only be performed
%                 for the voxels that are logical 1.
%       THRESH  - Threshold for the vascular. [Scalar] Default value is 1.5
%                   times the average CBV of the unaffected hemisphere.
%
%   OUTPUT:
%       OUT     - Mask image with vascular pixels removed. [Y x X] (logical)
%       V       - The vascular image. [Y x X] (grayscale)
%   
%
%   Ruogu Fang, 10/07/13
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if nargin < 2
    mask = ones(size(CBV));
end

if nargin < 3
    thresh = mean(CBV(mask))*1.5;
end

v = CBV;
v(v <= thresh) = 0;
out = mask;
out(find(v)) = 0;

end
