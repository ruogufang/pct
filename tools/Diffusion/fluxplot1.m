function [varargout] = fluxplot1(dif, grad, varargin)
%
% FLUXPLOT1    Plots the diffusiity and flux.
%
%    FLUXPLOT1( dif, grad ) plots de diffusivity (dif) and the flux (dif.*grad)
%    for a diffusing image.
%
%    "dif" and "grad" must have the same size and must be 1 dimensional.
%
%    FLUXPLOT1(...,dy ) plots also the image actualization dy (must be calculated outside)
%

axisoff = 0;
nimage = 2;

for i = 1 : length(varargin)
   if strcmp(varargin{i},'axisoff')
      axisoff = 1;
   elseif nimage == 2
      dy = varargin{i};
      nimage = 3;
   else
      error('Too many inputs')
   end
end

subplot(nimage,1 , 1);
plot(dif)
title('Diffusivity')


subplot(nimage, 1, 2)
plot(dif.*grad)
title('Flux')

if nimage == 3
   subplot(3, 1, 3)
   plot(dy)
   title('dy/dt')
end