function [ out ] = bilateral1(img, n, std_d, std_r)
%BILATERAL1 A 1D bilateral filter
%
%   USAGE:  OUT = BILATERAL1(IMG, N, STD_D, STD_R);
%
%   INPUT:
%       IMG     - A 1D array. [T x 1] or [1 x T]
%       N       - Filter size. Must be an odd number [Scalar]
%       STD_D   - The standard deviation of the domain filter [Scalar]
%       STD_R   - The standard deviation of the range filter [Scalar]
%
%   OUTPUT:
%       OUT     - The output array [T x 1] or [1 x T]
%
%   This function applies a 1D bilateral filter on the last dimension of a 4D
%   array. 
%
%   Kolbeinn Karlsson, 07/25/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Input validation
%TODO: Validate input

%Get the image size
T = length(img);

%Precalculate constants with respect to the loop
Sr = 2*std_r^2;

%Create the domain Gaussian weighting function
hs = (n-1)/2;
ind = -hs:hs;
Gd = exp( -(ind.^2)/(2*std_d^2));

%Pre-allocate output image
out = zeros(size(img));

%Apply the bilateral filter
for i = 1:T
    %Calculate bounds for the window
    imin = max(1,i-hs);
    imax = min(T,i+hs);
    Gdt = Gd((imin:imax)-i+hs+1);
    %Get the 1D temporal function
    I = squeeze(img(imin:imax));
    %Compute the temporal weighting function
    Gr = exp(- (I - img(i)).^2 / Sr);
    %Apply the temporal filter
    W = Gr(:) .* Gdt(:);
    out(i) = sum(W(:).*I(:))/sum(W(:));
end

end

