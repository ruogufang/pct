function difplot(x, time, step, tit, fig, varargin)
%  DIFPLOT    Plots the diffused image
%
%     DIFPLOT(x, time, step, tit, fig) plots the image "x" at the
%     subplot(1,2,2) of figure "fig". "tit" specifies the plot title.
%     "time" and "step" are the diffusion time and diffusion step show
%     as the xlabel of the plot.
%     DIFPLOT(..., 'imscale') plots the image using imagesc(x) intead
%     of image(x).

[imscale] = parse_inputs(varargin{:});

figure(fig);
subplot(1,2,2);

if ndims(x) == 3
   image(scale(x,[0 1]))
else
   if imscale
      imagesc(x);
      colorbar
   else
      image(x);
   end
end
xlabel(['Time =  ', num2str(time),',  Step = ', int2str(step)]);
title(tit);
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [imscale] = parse_inputs(varargin)

flag = zeros(1,nargin); % zero = par. nao usado / 1 = param . usado

imscale = 0;
for i = 1 : length(varargin)
   if strcmp(varargin(i), 'imscale')
      imscale = 1;
      flag(1,i) = 1;
   end
   if strcmp(varargin(i), 'null')
      flag(1,i) = -1;
   end
end

if any(flag==0)
   error('Unknown parameter.')
end

