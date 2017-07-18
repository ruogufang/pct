function y = aosiso1(x, d, t)
% AOSISO1   Aditive Operator Splitting Isotropic Interation (1D)
%
%    y = AOSISO1(x, d, t) calculates the new image "y" as the result of an
%    isotropic (scalar) diffusion iteration on image "x" with diffusivity 
%    "d" and steptime "t" using the AOS scheme.
%
%  - If "d" is constant the diffusion will be linear, if "d" is
%    a matrix the same size as "x" the diffusion will nonlinear.
%  - The stepsize "t" can be arbitrarially large, in contrast to the explicit
%    scheme, where t < 0.25. Using larger "t" will only affect the quality
%    of the diffused image. (Good choices are 5 < t < 20)
%  - "x" must be a 2D image.

% initialization
y = zeros(size(x));
p = zeros(size(d));

% Start  ======================================================================
% Operating on Columns
q = ( d(:,1:end-1)+d(:,2:end) );
p(:,1) = q(:,1);
p(:,end) = q(:,end);
p(:,2:end-1) = ( q(:,1:end-1) + q(:,2:end) );

a = 1 + t.*p';
b = -t.*q';

y = thomas(a,b,b,x')';
% End  ========================================================================


