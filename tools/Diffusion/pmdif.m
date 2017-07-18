function y = pmdif( u, lambda, sigma, stepsize, nosteps, varargin )
%
% PMDIF  Perona and Malik Diffusion
%
%    y = PMDIF( x, lambda, sigma, stepsize, steps, verbose, drawstep ) returns the 
%    image y as the result of the application of the Perona and Malick diffusion
%    to the image x. 
%
%    - Lambda is the contrast parameter (if the contrast is inferior to lambda the 
%      flux is increasing with the contrast and if the contrast is larger then lambda
%      the flux decreases as the contrast grows.
%      For time-variable lambda use lambda as a row vector and length(lambda)=number of steps.
%   -  Sigma is the space regularization parameter (the stardard deviation of the 
%      gaussian which will be convolve with the image before the gradient is calculated.
%      For time-variable sigma use sigma as a row vector and length(sigma)=number of steps.
%    -  The stepsize parameter can be a scalar (for constant stepsize) or a row vector for
%       variable stepsize. If stepsize is a row vector, length(stepsize) = number of steps.
%       For stability use stepsize < 0.25. If using 'aos' option stepsize can be any positive number.
%   -  The verbose parameter is a positive integer indicating the figure number 
%      where the output will be plotted.
%   -  The drawstep parameter indicates the number of steps between updating the 
%      displayed image.
%
%    y = PMDIF(..., 'imscale') uses imagesc instead of image to plot the diffusion.
%    y = PMDIF(..., 'grad','flux') plots also the gradient and flux images.
%    y = PMDIF(..., 'dfstep',n) only recalculates the diffusivities after n steps (increase speed)
%    y = PMDIF(..., 'aos') uses the AOS scheme to update the image. 
%         !!!!!   This allows the steptime to be arbitrarially big (t<inf) !!!!!
%    
%    See also: NLDIF, EEDIF, CEDIF, DIFFUSIVITY, GSDPLOT, FLUXPLOT
%

[verbose, drawstep, imscale, plotgrad, plotflux, dfstep, aos] = parse_inputs(varargin{:});

% Variable initialization
if ndims(u) > 2
   error('Input image must be grayscale.')
   return
end
dif_time = 0;
if strcmp(class(u),'double')
   y = u;
else
   y = double(u);
end

% Verifying inputs
[lambda, sigma, stepsize] = verify_inputs(lambda, sigma, stepsize, nosteps);

% Initial drawing
if verbose
   figure(verbose);
   subplot(1,2,1);
   if strcmp(imscale,'imscale')
      imagesc(y);
   else
      image(y);
   end
   colorbar
   difplot(u, 0, 0, 'Perona & Malick Diffusion', verbose, imscale);
end

for i=1:nosteps
   
   if mod(i-1,dfstep) == 0 % diffusivity recalc step
      % Calculate gradient of smoothed image
      if plotgrad
         figure(verbose+1)
         [gradx, grady] = gsdplot(y, sigma(i), 1, 'mod_angle');
      else
         [gradx, grady] = gsderiv(y, sigma(i), 1);
      end
      grad2 = gradx.^2 + grady.^2;
      grad = sqrt(grad2);
      % Calculate difusivity
      g = 1./(1+grad2./lambda(i)^2);
   end
   
   % Calculate dy/dt
   if aos
      if plotflux
         yo = y;
      end
      y = aosiso(y,g,stepsize(i)); % updating
   else
      dy = isodifstep(y, g);
      y = y + stepsize(i) * dy;  % updating
   end
   
   % Calculate diffusion time
   dif_time = dif_time + stepsize(i);
   
   % Plot actualization   
   if verbose & ~mod(i,drawstep)
      difplot(y, dif_time, i, 'Perona & Malick Diffusion', verbose, imscale);
      if plotflux
         figure(verbose+2*plotgrad+1)
         if aos
            fluxplot(g,grad, abs(y-yo), 'axisoff')
         else
            fluxplot(g,grad, abs(stepsize(i).*(dy)), 'axisoff')
         end
      end
   end
end % for i

% Last plot
if verbose
   difplot(y, dif_time, i, 'Perona & Malick Diffusion', verbose, imscale);
   if plotflux
      figure(verbose+2*plotgrad+1)
      fluxplot(g,grad, stepsize(i).*(dy))
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [verbose, drawstep, imscale, plotgrad, plotflux, dfstep, aos] = parse_inputs(varargin)
aos = 0;
dfsteppos = -1;
dfstep = 1;
verbose = -1;
drawstep = -1;
imscale = 'null';
plotgrad = 0;
plotflux = 0;

for i = 1 : length(varargin)
   flag = 0;
   if i == dfsteppos
      flag = 1;
   end
   if strcmp(varargin{i},'imscale')
      imscale = 'imscale';
      flag = 1;
   elseif strcmp(varargin{i},'grad') 
      plotgrad = 1;
      flag = 1;
   elseif strcmp(varargin{i},'flux')
      plotflux = 1;
      flag = 1;
   elseif strcmp(varargin{i},'dfstep')
      dfstep = varargin{i+1};
      flag = 1;
      dfsteppos = i+1;
   elseif strcmp(varargin{i},'aos')
      aos = 1;
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
function [nlambda, nsigma, nstepsize] = verify_inputs(lambda, sigma, stepsize, nosteps)

% Verifying lambda
if sum(size(lambda)>1) == 0 % lambda is constant
   nlambda = linspace(lambda,lambda,nosteps);
else
   if sum(size(lambda)>1) > 1
      error('lambda must be a row vector')
      return
   end
   if length(lambda)~=nosteps
      error('length(lambda) must be equal to number of steps')
      return
   end
   nlambda = lambda;
end

% Verifying simga
if sum(size(sigma)>1) == 0 % sigma is constant
   nsigma = linspace(sigma,sigma,nosteps);
else
   if sum(size(sigma)>1) > 1
      error('sigma must be a row vector')
      return
   end
   if length(sigma)~=nosteps
      error('length(sigma) must be equal to number of steps')
      return
   end
   nsigma = sigma;
end

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

