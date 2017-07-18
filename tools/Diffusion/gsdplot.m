function [varargout] = gsdplot(u, sigma, order, varargin)
%
% GSDPLOT    Plots the gaussian smoothed derivates.
%
%    GSDPLOT( u, sigma, order ) plots the "order"-th smoothed derivates of
%    the image "u". The smoothing is performed prior to the derivation by a 2D
%    convolution with a "sigma" std.dev. gaussian. "u" must have two dimensions.
%
%    If used with left arguments, GSDPLOT returns the same output as GSDERIV.
%    The number of outputs depends on the "order". For "order"=0 there is only one
%    output corresponding to the smoothed input. For "order"=1 there are two outputs,
%    respectivelly the x and y derivates. For "order"=2 there are three outputs,
%    corresponding respectivelly to the 2º x derivate, the crossed x y derivate and
%    the 2º y derivate.
%
%    GSDPLOT(..,'abs') plots the absolute values.
%    GSDPLOT(x,sigma,1,'mod_angle') plots modulus and angle instead of components.
%         for order = 2 the modulus and angle are taken from d2u/dx2 and d2u/dy2.
%         the parameter 'mod_angle' cannot be used with order = 0.
%
%
%    Examples:
%       GSDPLOT(x,1,1)
%       [gx gy] = GSDPLOT(x,1,1, 'mod_angle')
%       [gx2 gxy gy2] = GSDPLOT(x,1,2,'abs')
%

[absolute, mod_angle, axisoff] = parse_inputs(varargin{:});

switch order
case 0
   if mod_angle
      error('Modulus / Angle images are only availables to order 1 & 2')
   end
   y = gsderiv(u, sigma, 0);
   nimage(y, absolute, 0, axisoff);

case 1
   [y(:,:,1) y(:,:,2)] = gsderiv(u, sigma, 1);
   modulus = nimage(y, absolute, mod_angle, axisoff);

case 2
   [y(:,:,1) y(:,:,2) y(:,:,3)] = gsderiv(u, sigma, 2);
   modulus = nimage(y, absolute, mod_angle, axisoff);
   
end
drawnow;

if nargout
   for i = 1 : size(y,3)
      varargout(i) = {y(:,:,i)};
   end
end

if nargout-size(y,3) == 1 & mod_angle
   varargout(nargout) = {modulus};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modulus = nimage(x, absolute, mod_angle, axisoff)
modulus = 0;

switch size(x,3)
case 1 % order zero
   absimage(x, absolute);
   title('Smoothed')
   if axisoff
      axis off
   end
case 2 % order 1
   if mod_angle
      modulus = sqrt(x(:,:,1).^2+x(:,:,2).^2);
      absimage(modulus, absolute);
      title('Grad. Magnitude')
      colorbar
      if axisoff
         axis off
      end
      if mod_angle == 1
         figure(gcf+1);
         ang = 180/pi*atan(x(:,:,2)./(x(:,:,1)+eps));
         imagesc(ang);
         colormap(colormapc(4));
         title('Grad. Angle')
         colorbar
         if axisoff
            axis off
         end
      end
   else
      subplot(1,2,1);
      absimage(x(:,:,1), absolute);
      title('Grad x')
      colorbar
      if axisoff
         axis off
      end
      subplot(1,2,2);
      absimage(x(:,:,2), absolute);
      title('Grad y')
      colorbar
      if axisoff
         axis off
      end
   end
case 3  % order 2
   if mod_angle
      subplot(1,2,1);
      modulus = sqrt(x(:,:,1).^2+x(:,:,3).^2);
      absimage(modulus, absolute);
      title('Grad. Magnitude')
      colorbar
      if axisoff
         axis off
      end
      if mod_angle == 1
         figure(gcf+1);
         ang = 180/pi*atan(x(:,:,3)./(x(:,:,1)+eps));
         imagesc(ang);
         colormap(colormapc(4));
         title('Grad. Angle')
         colorbar
         if axisoff
            axis off
         end
      end
   else
      subplot(2,2,1);
      absimage(x(:,:,1), absolute);
      title('dx2')
      colorbar
      if axisoff
         axis off
      end
      subplot(2,2,2);
      absimage(x(:,:,2), absolute);
      title('dxy')
      colorbar
      subplot(2,2,3);
      nimage(x(:,:,3), absolute);
      title('dy2')
      colorbar
      if axisoff
         axis off
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function absimage(x, absolute)
if absolute
   imagesc(abs(x));
else
   imagesc(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [absolute, mod_angle, axisoff] = parse_inputs(varargin)

axisoff = 0;
absolute = 0;
mod_angle = 0;

for i = 1 : length(varargin)
   if strcmp(varargin{i}, 'abs')
         absolute = 1;
      elseif strcmp(varargin{i}, 'mod_angle')
         mod_angle = 1;
      elseif strcmp(varargin{i}, 'modulus')
         mod_angle = 2;
      elseif strcmp(varargin{i}, 'axisoff')
         axisoff = 1;
      else
         error('Unknown parameter.')
      end
end