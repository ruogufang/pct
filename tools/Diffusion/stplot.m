function [varargout] = stplot(u, sigma, rho, varargin)
%
% STPLOT    Plots the structure tensor.
%
%    STPLOT( u, sigma, rho ) plots the structure tensor of image u. The structure
%    tensor (S) is given by:
%
%       S(G) =  / R(Gx*Gx)  R(Gx*Gy) \    
%               \ R(Gx*Gy)  R(Gy*Gy) /       
%
%       where:
%             G = [Gx Gy]' is grad(us) 
%             R = 2D rho standard deviation gaussian smoothing
%             grad(us) = the gradient of the image u smoothed by a sigma std.dev.
%               2D gaussian.
%
%    If used with left arguments, STPLOT returns [Vx, Vy] the components of the 
%    eigenvectors of S. (the other eigenvector has components [-Vy, Vx])
%

% Gradient of smoothed u
[gradx, grady] = gsderiv(u, sigma, 1);

% Structure tensor
S11 = gsderiv(gradx.^2,rho,0);
S12 = gsderiv(gradx.*grady,rho,0);
S22 = gsderiv(grady.^2,rho,0);

% Struc. Direction
vx = 2*S12;
vy = S22 - S11 + sqrt( (S11-S22).^2+4*S12.^2 );

ang_st = 180/pi*atan(vx./(vy+eps));
imagesc(ang_st)
colorbar
title('Structure orientation (degrees)')

switch nargout
case 1
   varargout{1} = vx;
case 2
   varargout{1} = vx;
   varargout{2} = vy;
end