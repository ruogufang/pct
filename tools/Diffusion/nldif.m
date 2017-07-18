function y = nldif( u, lambda, sigma, m, stepsize, nosteps, varargin)
%
% NLDIF  Nonlinear Diffusion
%
%    y = NLDIF( u, lambda, sigma, m, stepsize, steps, verbose, drawstep ) returns the 
%    image y as the result of the application of the Non Linear diffusion
%    to the image u.
%
%    dy/dt = div( d(grad(y)).grad(y) ),  y(0) = u,  y(t+T) = y(t) + T*dy/dt
%
%    The diffusivity function (d) is given by:
%
%       d(g) = 1 - exp( -Cm / (g/lambda)^m ), g > 0
%              1                            , g <=0
%
%    -  The constant Cm is calculated to make the flux (g*d(g)) ascending for g < lambda 
%       and descending for g >= lambda.
%    -  Lambda is the contrast parameter (if the gradient is inferior to lambda the 
%       flux is increasing with the gradient and if the gradient is larger then lambda
%       the flux decreases as the gradient grows.
%       For time-variable lambda use lambda as a row vector and length(lambda)=number of steps.
%    -  Sigma is the space regularization parameter (the stardard deviation of the 
%       gaussian which will be convolved with the image before the gradient is calculated.
%       For time-variable sigma use sigma as a row vector and length(sigma)=number of steps.
%    -  'm' defines the speed the difusivity (and the flux) changes for a variation in the
%       gradient. Big values of 'm' make the flux the change quickly. 'm' must be bigger than 1.
%       A good choice for m is 8 < m < 16.
%    -  The stepsize parameter can be a scalar (for constant stepsize) or a row vector for
%       variable stepsize. If stepsize is a row vector, length(stepsize) = number of steps.
%       For stability use stepsize < 0.25. If using 'aos' option stepsize can be any positive number.
%    -  The steps parameter indicates the number of iterations to be performed.
%    -  The verbose parameter is a positive integer indicating the figure number 
%       where the output will be plotted.
%    -  The drawstep parameter indicates the number of steps between updating the 
%       displayed image.
% 
%    y = NLDIF(..., 'imscale') uses imagesc instead of image to plot the diffusion.
%    y = NLDIF(..., 'grad','flux') plots also the gradient and flux images.
%    y = NLDIF(..., 'dfstep',n) only recalculates the diffusivities after n steps (increase speed)
%    y = NLDIF(..., 'aos') uses the AOS scheme to update the image. 
%         !!!!!   This allows the steptime to be arbitrarially big (t<inf) !!!!!
%    y = NLDIF(..., 'cedif',{p, lambda, sigma, rho, m, stepsize, steps}) performs an coherence enhancing diffusion 
%         after every p steps using the given parameters. To use in the cedif the same parameters used in the
%         NLDIF type 'sp' in the corresponding place (only works for lambda, sigma and m).
%
%    See also: PMDIF, CEDIF, DIFFUSIVITY, GSDPLOT, FLUXPLOT
[verbose, drawstep, imscale, plotgrad, plotflux, dfstep, aos, ch, eed, eedpar] = parse_inputs(varargin{:});

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
[lambda, sigma, stepsize] = verify_inputs(lambda, sigma, stepsize, nosteps);

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
   if ch ~= 0
      xlabel( ['Channel = ', int2str(ch)] )
   end
   difplot(y, 0, 0, 'Non Linear Diffusion', verbose, imscale);
end

% Calculate Cm constant
Cm = Cmcalc(m);
for i=1:nosteps
   
   if eed ~= 0
      eestep = ( mod(i-1,eedpar{1})==0 & i~=1 ); % Current step will use Coherence enhancing diffusion.
   else
      eestep = 0;
   end
      
   if eestep
      if plotflux
         yo = y;
      end
      [y, eedtime] = eedfun(y,eedpar,lambda(i), sigma(i),m);
   else
      if mod(i-1,dfstep) == 0 % diffusivity recalc step
         % Calculate gradient of smoothed image
         if plotgrad
            figure(verbose+1)
            [gradx, grady, grad] = gsdplot(y, sigma(i), 1, 'modulus','axisoff');
        else
            [gradx, grady] = gsderiv(y, sigma(i), 1);
            grad = sqrt(gradx.^2 + grady.^2);
         end
         
         % Calculate difusivity
         g = 1 - exp(-Cm./( (grad+eps)./lambda(i)).^m);  % grad + eps to avoid division by zero (eps -> 0)
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
      
   end % else eestep
   
   % Calculate diffusion time
   if eestep
      dif_time = dif_time + eedtime;
   else
      dif_time = dif_time + stepsize(i);
   end
         
   % Plot actualization   
   if verbose & ~mod(i,drawstep)
      difplot(y, dif_time, i, 'Non Linear Diffusion', verbose, imscale);
      if plotflux
         figure(verbose+plotgrad+1)
         if aos | eestep
            fluxplot(g,grad, abs(y-yo), 'axisoff')
         else
            fluxplot(g,grad, abs(stepsize(i).*(dy)), 'axisoff')
         end
      end
   end
   
end % for i = 1 : nosteps

% Last plot
if verbose
   difplot(y, dif_time, i, 'Non Linear Diffusion', verbose, imscale);
   if plotflux
      figure(verbose+plotgrad+1)
      if aos | eestep
         fluxplot(g,grad, abs(y-yo), 'axisoff')
      else
         fluxplot(g,grad, abs(stepsize(i).*(dy)), 'axisoff')
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [verbose, drawstep, imscale, plotgrad, plotflux, dfstep, aos, ch, eed, eedpar] = parse_inputs(varargin)
eed = 0;
eedpar = 0;
eedpos = -1;
ch = 0;
chpos = -1;
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
   if i == dfsteppos | i == chpos | i == eedpos
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
   elseif strcmp(varargin{i},'ch')
      ch = varargin{i+1};
      flag = 1;
      chpos = i+1;
   elseif strcmp(varargin{i},'cedif')
      eed = 1;
      eedpar = varargin{i+1};
      if ~length(eedpar)==7 | ~iscell(eedpar)
         error('CEDIF parameters incorrect. CEDIF parameters must be passed as a 1x7 cell array.')
      end
      eedpos = i+1;
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
function Cm = Cmcalc(m)
if m <= 1
   error('Use m > 1')
   return
else
   Cm = fzero(strcat('1-exp(-x)-x*exp(-x)*',num2str(m)),[1e-10 1e100]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, eedtime] = eedfun(x, eedpar, lambda2, sigma2, m2)
p = eedpar{1};
lambda = eedpar{2};
sigma = eedpar{3};
rho = eedpar{4};
m = eedpar{5};
stepsize = eedpar{6};
steps = eedpar{7};

if strcmp(lambda,'sp')
   lambda = lambda2;   
end
if strcmp(sigma,'sp')
   sigma = sigma2;   
end
if strcmp(m,'sp')
   m = m2;   
end

if length(stepsize)==1
   eedtime = stepsize*steps;
else
   eedtime = sum(stepsize);
end

y = cedif(x, lambda, sigma, rho, m, stepsize, steps);