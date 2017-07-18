function y = ldif( u,stepsize, nosteps, varargin)
%
% LDIF  Linear Diffusion
%
%    y = LDIF( u, stepsize, steps, verbose, drawstep ) returns the image y as the 
%    result of the application of the linear diffusion to the image u.
%
%    dy/dt = div( 1 . grad(y) ),  y(0) = u,  y(t+T) = y(t) + T*dy/dt
%
%    -  The stepsize parameter can be a scalar (for constant stepsize) or a row vector for
%       variable stepsize. If stepsize is a row vector, length(stepsize) = number of steps.
%    -  The steps parameter indicates the number of iterations to be performed.
%    -  The verbose parameter is a positive integer indicating the figure number 
%       where the output will be plotted.
%    -  The drawstep parameter indicates the number of steps between updating the 
%       displayed image.
% 
%    y = LDIF(..., 'imscale') uses imagesc instead of image to plot the diffusion.
%
%    See also: PMDIF, NLDIF, EEDIF, CEDIF, DIFFUSIVITY, GSDPLOT, FLUXPLOT
[verbose, drawstep, imscale] = parse_inputs(varargin{:});

% Variable initialization
dif_time = 0;
if ndims(u) > 2
   error('Input image must be grayscale.')
   return
end
if strcmp(class(u),'double')
   y = u;
else
   y = double(u);
end

% Verifying inputs
[stepsize] = verify_inputs(stepsize, nosteps);

% Initial Drawing
if verbose
   figure(verbose);
   subplot(1,2,1); 
   if strcmp(imscale,'imscale')
      imagesc(y);
   else
      image(y);
   end
   colorbar
   title('Original Image'); drawnow;
   difplot(u, 0, 0, 'Linear Diffusion', verbose, imscale);
end

% Diffusion 
for i=1:nosteps
   
   % Calculate Kernel
   sd = sqrt(2*stepsize(i));
   kern_index = -ceil(3*sd):ceil(3*sd);
   kern = normpdf(kern_index,0,sd);
   
   % Calculate new image
   y = conv2br(kern,kern,y);
   
   % Calculate diffusion time
   dif_time = dif_time + stepsize(i);  
      
   % Plot actualization   
   if verbose & ~mod(i,drawstep)
      difplot(y, dif_time, i, 'Linear Diffusion', verbose, imscale);
   end
end

% Last plot
if verbose
difplot(y, dif_time, i, 'Linear Diffusion', verbose, imscale);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [verbose, drawstep, imscale] = parse_inputs(varargin)
verbose = -1;
drawstep = -1;
imscale = 'null';

for i = 1 : length(varargin)
   flag = 0;
   if strcmp(varargin{i},'imscale')
      imscale = 'imscale';
      flag = 1;
   end
   if flag == 0 & verbose == -1
      verbose = varargin{i};
      flag = 1;
   end
   if flag == 0 & drawstep == -1
      drawstep = varargin{i};
      flag = 1;
   end
   if flag == 0
      error('Too many parameters !')
      return
   end
end

if verbose == -1
   verbose = 0;
end

if drawstep == -1
   if verbose == 0
      drawstep = 0;
   else
      drawstep = 1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nstepsize] = verify_inputs(stepsize, nosteps)

% Verifying stepsize
if sum(size(stepsize)>1) == 0 % constant stepsize
   nstepsize = linspace(stepsize,stepsize,nosteps);
else
   if sum(size(stepsize)>1) > 1
      error('stepsize must be a row vector')
      return
   end
   if length(stepsize)~=nosteps
      error('length(stepsize) must be equal to number of steps')
      return
   end
   nstepsize = stepsize;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = normpdf(x,m,s)
y = exp( -.5*(( (x-m)./s).^2) ) ./ (sqrt(2*pi)*s);