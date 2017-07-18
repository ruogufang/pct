function [varargout] = fluxplot(dif, grad, varargin)
%
% FLUXPLOT    Plots the diffusiity and flux.
%
%    FLUXPLOT( dif, grad ) plots de diffusivity (dif) and the flux (dif.*grad)
%    for a diffusing image.
%
%    "dif" and "grad" must have the same size and must be 2 dimensional.
%
%    FLUXPLOT(...,dy ) plots also the image actualization dy (must be calculated outside)
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

subplot(1, nimage, 1);
imagesc(dif)
title('Diffusivity')
colorbar
if axisoff
   axis off
end

subplot(1, nimage, 2)
imagesc(dif.*grad)
title('Flux')
colorbar
if axisoff
   axis off
end

if nimage == 3
   subplot(1,3,3)
   imagesc(dy)
   title('dy/dt')
   colorbar
   if axisoff
      axis off
   end
end