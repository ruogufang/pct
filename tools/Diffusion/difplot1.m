function difplot1(x, time, step, tit, fig, ax, varargin)
%  DIFPLOT1    Plots the diffused image - 1D
%
%     DIFPLOT1(x, time, step, tit, fig) plots the image "x" at the
%     subplot(1,2,2) of figure "fig". "tit" specifies the plot title.
%     "time" and "step" are the diffusion time and diffusion step show
%     as the xlabel of the plot.
%     DIFPLOT1(..., 'imscale') plots the image using imagesc(x) intead
%     of image(x).

[imscale] = parse_inputs(varargin{:});

figure(fig);
subplot(2,1,2);

plot(x)

if ~imscale
    axis(ax)
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

