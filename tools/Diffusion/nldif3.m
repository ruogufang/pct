function y = nldif( u, lambda, sigma, m, stepsize, nosteps, varargin)
%
% NLDIF3  3D Nonlinear Diffusion
%
%    y = NLDIF3( u, lambda, sigma, m, stepsize, steps, verbose, drawstep ) returns the 
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
%       For stability use stepsize < 0.16667. 
%    -  The steps parameter indicates the number of iterations to be performed.
%    -  The verbose parameter is a positive integer indicating the figure number 
%       where the output will be plotted.
%    -  The drawstep parameter indicates the number of steps between updating the 
%       displayed image.
% 
%    y = NLDIF(..., 'imscale') uses imagesc instead of image to plot the diffusion.
%    y = NLDIF(..., 'grad') plots also the gradient image.
%    y = NLDIF(..., 'dfstep',n) only recalculates the diffusivities after n steps (increase speed)
%
%    See also: PMDIF, EEDIF, CEDIF, DIFFUSIVITY, GSDPLOT, FLUXPLOT
[verbose, drawstep, imscale, plotgrad, dfstep] = parse_inputs(varargin{:});

% Variable initialization
dif_time = 0;
if strcmp(class(u),'double')
    y = u;
else
    y = double(u);
end

% Verifying inputs
[lambda, sigma, stepsize] = verify_inputs(lambda, sigma, stepsize, nosteps);

% Initial Drawing
slicen = round(size(y,3)/2); % slice to be used in plots
if verbose
    figure(verbose);
    subplot(1,2,1); 
    if strcmp(imscale,'imscale')
        imagesc(y(:,:,slicen));
    else
        image(y(:,:,slicen));
    end
    colorbar
    title('Original Image Slice'); drawnow;
    difplot(y(:,:,slicen), 0, 0, 'Non Linear Diffusion', verbose, imscale);
end

% Calculate Cm constant
Cm = Cmcalc(m);
for i=1:nosteps
    
    if mod(i-1,dfstep) == 0 % diffusivity recalc step
        % Calculate gradient of smoothed image  
        smoothsize = floor(sigma(i));
        if mod(smoothsize,2) == 0;
            smoothsize = smoothsize + 1;
        end
        if smoothsize > 1; 
            ys = smooth3(y,'box',smoothsize);
        else
            ys = y;
            %disp('a')
        end
        [gradx, grady, gradz] = gradient(ys);
        grad = sqrt(gradx.^2 + grady.^2 + gradz.^2);
        
        if plotgrad
            figure(verbose+1)
            clf
            imagesc(grad(:,:,slicen));
            axis off
            colorbar
            Title('3D Gradient')
        end
        
        
        % Calculate difusivity
        g = 1 - exp(-Cm./( (grad+eps)./lambda(i)).^m);  % grad + eps to avoid division by zero (eps -> 0)
    end
    
    % Calculate dy/dt
    dy = isodifstep3(y, g);
    y = y + stepsize(i) * dy;  % updating
    
    % Plot actualization   
    dif_time = dif_time + stepsize(i);
    if verbose & ~mod(i,drawstep)
        difplot(y(:,:,slicen), dif_time, i, 'Non Linear Diffusion', verbose, imscale);
    end
    
end % for i = 1 : nosteps

% Last plot
if verbose
    difplot(y(:,:,slicen), dif_time, i, 'Non Linear Diffusion', verbose, imscale);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [verbose, drawstep, imscale, plotgrad, dfstep] = parse_inputs(varargin)
dfsteppos = -1;
dfstep = 1;
verbose = -1;
drawstep = -1;
imscale = 'null';
plotgrad = 0;

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
    elseif strcmp(varargin{i},'dfstep')
        dfstep = varargin{i+1};
        flag = 1;
        dfsteppos = i+1;
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