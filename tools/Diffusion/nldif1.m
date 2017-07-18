function y = nldif( u, lambda, sigma, m, stepsize, nosteps, varargin)
%
% NLDIF  1D Nonlinear Diffusion
%
%    y = NLDIF1( u, lambda, sigma, m, stepsize, steps, verbose, drawstep ) returns the 
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
%    y = NLDIF(..., 'scale') re-scales the diffused signal for it to best fit the graph (no change in the real values).
%    y = NLDIF(..., 'grad','flux') plots also the gradient and flux images.
%    y = NLDIF(..., 'dfstep',n) only recalculates the diffusivities after n steps (increase speed)
%    y = NLDIF(..., 'aos') uses the AOS scheme to update the image. 
%         !!!!!   This allows the steptime to be arbitrarially big (t<inf) !!!!!
%
%    See also: PMDIF, CEDIF, DIFFUSIVITY, GSDPLOT, FLUXPLOT
[verbose, drawstep, imscale, plotgrad, plotflux, dfstep, aos] = parse_inputs(varargin{:});

% Variable initialization
dif_time = 0;
if sum(size(u)>1) > 1 
    error('NLDIF1 works only with 1D signals.')
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
    subplot(2,1,1); 
    plot(y)
    ini_axis = axis;
    ini_axis = axis.*[1 1 1.1 1.1];
    axis(ini_axis);
    title('Original Signal'); drawnow;
    difplot1(y, 0, 0, 'Non Linear Diffusion', verbose, ini_axis, imscale);
end

% Calculate Cm constant
Cm = Cmcalc(m);
for i=1:nosteps
    
    if mod(i-1,dfstep) == 0 % diffusivity recalc step
        % Calculate gradient of smoothed image  
        smoothsize = floor(sigma(i));
        if smoothsize > 1; 
            coef = ones(1,smoothsize)/smoothsize;
            ys = conv2(y,coef,'same');
            %ys = conv2br(1,coef,y);
        else
            ys = y;
        end
        grad = gradient(ys);
        
        if plotgrad
            figure(verbose+1)
            clf
            plot(grad);
            Title('Gradient')
        end
        % Calculate difusivity
        g = 1 - exp(-Cm./( (grad+eps)./lambda(i)).^m);  % grad + eps to avoid division by zero (eps -> 0)
    end
    
    % Calculate dy/dt
    if aos
        if plotflux
            yo = y;
        end
        y = aosiso1(y,g,stepsize(i)); % updating
    else
        dy = isodifstep1(y, g);
        y = y + stepsize(i) * dy;  % updating
    end
    
    
    % Calculate diffusion time
    dif_time = dif_time + stepsize(i);

    % Plot actualization   
    if verbose & ~mod(i,drawstep)
        difplot1(y, dif_time, i, 'Non Linear Diffusion', verbose, ini_axis, imscale);
        if plotflux
            figure(verbose+plotgrad+1)
            if aos 
                fluxplot1(g,grad, abs(y-yo), 'axisoff')
            else
                fluxplot1(g,grad, abs(stepsize(i).*(dy)), 'axisoff')
            end
        end
    end
    
end % for i = 1 : nosteps

% Last plot
if verbose
    difplot1(y, dif_time, i, 'Non Linear Diffusion', verbose, ini_axis, imscale);
    if plotflux
        figure(verbose+plotgrad+1)
        if aos
            fluxplot1(g,grad, abs(y-yo), 'axisoff')
        else
            fluxplot1(g,grad, abs(stepsize(i).*(dy)), 'axisoff')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
   if strcmp(varargin{i},'scale')
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
      error('Too many \ Unknown parameters !')
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
