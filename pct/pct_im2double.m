function [ out ] = pct_im2double(img, lo, hi)
%PCT_IM2DOUBLE Scales the color values of an image to [0,1]
%
%   This function takes in an image/map of any scale, windows it from LO,HI to
%   [0,1]. The result can easily be saved as a JPEG or any image format of
%   choice.
%
%   USAGE: OUT = PCT_IM2DOUBLE(IMG, LO, HI);
%
%   PRE:
%       IMG     - An image/map [M x N]
%       LO      - The bottom of the contrast window [Scalar]
%       HI      - The top of the contrast window [Scalar]
%
%   POST:
%       OUT     - The windowed output image on a double format [M x N]
%
%   Kolbeinn Karlsson, 06/26/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%

%Convert to double
img = double(img);

%Get image size
[M N] = size(img);
len = M*N;

%Pre-allocate output
out = zeros(M,N);

%Window
for i = 1:len
    if img(i) > hi
        out(i) = 1.0;
    elseif img(i) < lo
        out(i) = 0.0;
    else
        out(i) = (img(i)-lo) * 1/(hi-lo);
    end
end

end

