function y = cedif( u, lambda, sigma, rho, m, stepsize, nosteps, varargin)
%
% CEDIF  Coherence Enhancing Diffusion
%
%    y = CEDIF( u, lambda, sigma, rho, m, stepsize, steps, verbose, drawstep ) returns the 
%    image y as the result of the application of the Coherence Enhancing anisotropic 
%    diffusion to the image u.
%
%    dy/dt = div( D(S(grad(us))).grad(y) ),  y(0) = u,  y(t+T) = y(t) + T*dy/dt
%
%    The diffusion tensor (D) is given by:
%       D(S) =  / c1+c2+[(c2-c1)*(s11-s22)/a]         (c2-c1)*s12/a         \
%               \      (c2-c1)*s12/a            c1+c2-[(c2-c1)*(s11-s22)/a] /
%
%         where: a = sqrt{ (s11-s22)^2 + 4*(s12)^2 },
%                c1 and c2 are the diffusivities perpendicular and parallel to the structure orientation,
%                c1 = gamma                                             
%                c2 = / gamma+(1-gamma)*exp(-lambda/(mi1-mi2)^m ) , mi ~= mi2   
%                     \ gamma                                    , mi1 = mi2   
% 
%    The structure tensor (S) is given by:
%       S(G) =  / R(Gx*Gx)  R(Gx*Gy) \    , where G = [Gx Gy]' is the gradient of us
%               \ R(Gx*Gy)  R(Gy*Gy) /      R(.) is smoothing with a 2D rho std.dev. gaussian.
%
%    -  Lambda is the coherence parameter (if the coherence is inferior to lambda the 
%       flux is increasing with the coherence and if the coherence is larger then lambda
%       the flux decreases as the coherence grows.
%       For time-variable lambda use lambda as a row vector and length(lambda)=number of steps.
%    -  Sigma is the space regularization parameter (the stardard deviation of the 
%       gaussian which will be convolved with the image before the gradient is calculated.
%       For time-variable sigma use sigma as a row vector and length(sigma)=number of steps.
%    -  Rho is the space regularization parameter for structure S(G). A rho std.dev. gaussian will be 
%       convolved with the elements of the structure tensor before the diffusion tensor is calculated
%       For time-variable rho use rho as a row vector and length(rho)=number of steps.
%    -  The stepsize parameter can be a scalar (for constant stepsize) or a row vector for
%       variable stepsize. If stepsize is a row vector, length(stepsize) = number of steps.
%       For stability use stepsize < 0.25
%    -  The steps parameter indicates the number of iterations to be performed.
%    -  The verbose parameter is a positive integer indicating the figure number 
%       where the output will be plotted.
%    -  The drawstep parameter indicates the number of steps between updating the 
%       displayed image.
%    -  gamma is a small positive to keep D positive definite (gamma = 0.01) 
%       use y=CEDIF(...,'gamma', value) to change gamma
%    -  'm' defines the speed the difusivity parallel to the structure changes for a variation in 
%       the coherence. Big values of 'm' make the diffusivity change quickly. 'm' must be bigger than 1.
%
%    y = CEDIF(..., 'imscale') uses imagesc instead of image to plot the diffusion.
%    y = CEDIF(..., 'grad','struc','dif', 'dy') plots also the gradient, structure, diffusivity
%          perpendicular to gradient and dy/dt images.
%    y = CEDIF(..., 'dfstep',n) only recalculates the diffusion tensor after n steps (increase speed)
%
%    See also: PMDIF, NLDIF, EEDIF, DIFFUSIVITY, GSDPLOT, FLUXPLOT

[verbose, drawstep, imscale, plotgrad, plotstruc, plotdif, plotdy, dfstep, ori, gamma] = parse_inputs(varargin{:});

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
[lambda, sigma, rho, stepsize] = verify_inputs(lambda, sigma, rho, stepsize, nosteps);

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
   title('Original Image');
   difplot(y, 0, 0, 'Coherence Enhancing Diffusion', verbose, imscale);
   drawnow;
end

% Calculate Cm constant
Cm = Cmcalc(m);


for i=1:nosteps
   
   if mod(i-1,dfstep) == 0 % diffusivity recalc step
      % Calculate gradient of smoothed image
      if plotgrad
         figure(verbose+1)
         if ori
            [gradx, grady] = oridplot( gsderiv(y,sigma(i),0) );
         else
            [gradx, grady] = gsdplot(y, sigma(i), 1, 'modulus');
         end
      else
         if ori
            [gradx, grady] = orideriv( gsderiv(y,sigma(i),0) );
         else
            [gradx, grady] = gsderiv(y, sigma(i), 1);
         end
      end
      
      % Calculate structure tensor
      s11 = gsderiv(gradx.^2, rho(i), 0);
      s12 = gsderiv(gradx.*grady, rho(i), 0);
      s22 = gsderiv(grady.^2, rho(i), 0);
      
      %s1_p_s2 = s11 + s22;
      s1_m_s2 = s11 - s22;
      
      % Structure tensor autovalues     
      %alfa = sqrt( (s1 - s2).^2 + 4*s12.^2 );
      %mi1 = .5*(s11+s22-alfa); 
      %mi2 = .5*(s11+s22+alfa); 
      alfa = sqrt( (s1_m_s2).^2 + 4*s12.^2 ); % alfa = mi1 - mi2 = coherence
      
      % Plot structure tensor
      if plotstruc
         figure(verbose+plotgrad+1)
         imagesc(alfa);
         colorbar('vert');
         title('Structure Coherence (mi1-mi2)')
         figure(gcf+1)
         vx = 2*s12;
         vy = - s1_m_s2 + alfa;
         ang_st = 180/pi*atan(vx./(vy+eps));
         imagesc(ang_st)
         colormap(colormapc(4));
         colorbar('vert')
         title('Structure orientation (degrees)')
         drawnow;
      end

      
      % Diffusion tensor autovalues   
      c1 = gamma;
      c2 = gamma + (1-gamma)*exp(-Cm./ ( (alfa+eps)./lambda(i) ).^m );
      if plotdif
         figure(verbose+plotgrad+2*plotstruc+1)
         imagesc(c2)
         colorbar
         title('Diffusivity perp. to gradient')
         drawnow;
      end

      
      % Calculate diffusion tensor components
      c1_p_c2 = c1 + c2;
      c2_m_c1 = c2 - c1;
      
      dd = ( c2_m_c1 .* s1_m_s2 )./(alfa);
      
      d11 = .5 * ( c1_p_c2 + dd );
      d12 = -(c2_m_c1).*s12./(alfa);
      d22 = .5 * ( c1_p_c2 - dd);
   end
   
   
   % Calculate dy/dt
   if ori
      [dx dy] = orideriv(y);
      j1 = d22.*dx + d12.*dy;
      j2 = d12.*dx + d11.*dy;
      dy = orideriv(j1,1) + orideriv(j2,2);
   else
      dy = anidifstep(y, d11, d12, d22);
      %dy = anidifstep(y, d22, d12, d11);
   end
   
   % Calculate diffusion time
   dif_time = dif_time + stepsize(i);
   
   % Image updating
   y = y + stepsize(i) * dy;
   
   % Plot actualization   
   if verbose & ~mod(i,drawstep)
      difplot(y, dif_time, i, 'Coherence Enhancing Diffusion', verbose, imscale);
      if plotdy
         figure(verbose+plotgrad+2*plotstruc+plotdif+1)
         imagesc(abs(stepsize(i).*(dy)))
         colorbar('vert');
         title('dy/dt')
         drawnow;
      end
   end
end

% Last plot
if verbose
   difplot(y, dif_time, i, 'Coherence Enhancing Diffusion', verbose, imscale);
   if plotdy
      figure(verbose+plotgrad+2*plotstruc+plotdif+1)
      imagesc(abs(stepsize(i).*(dy)))
      colorbar('vert');
      drawnow;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [verbose, drawstep, imscale, plotgrad, plotstruc, plotdif, plotdy, dfstep, ori, gamma] = parse_inputs(varargin)
gamma = .01;
gammapos = -1;
ori = 0;
dfsteppos = -1;
dfstep = 1;
verbose = -1;
drawstep = -1;
imscale = 'null';
plotgrad = 0;
plotstruc = 0;
plotdif = 0;
plotdy = 0;

for i = 1 : length(varargin)
   flag = 0;
   if i==dfsteppos | i==gammapos
      flag = 1;
   end
   if strcmp(varargin{i},'imscale')
      imscale = 'imscale';
      flag = 1;
   elseif strcmp(varargin{i},'grad') 
      plotgrad = 1;
      flag = 1;
   elseif strcmp(varargin{i},'struc')
      plotstruc = 1;
      flag = 1;
   elseif strcmp(varargin{i},'dif')
      plotdif = 1;
      flag = 1;
   elseif strcmp(varargin{i},'dy')
      plotdy = 1;
      flag = 1;
   elseif strcmp(varargin{i},'dfstep')
      dfstep = varargin{i+1};
      flag = 1;
      dfsteppos = i+1;
   elseif strcmp(varargin{i},'ori')
      ori = 1;
      flag = 1;
   elseif strcmp(varargin{i},'gamma')
      gamma = varargin{i+1};
      flag = 1;
      gammapos = i+1;
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
function [nlambda, nsigma, nrho, nstepsize] = verify_inputs(lambda, sigma, rho, stepsize, nosteps)

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

% Verifying rho
if sum(size(rho)>1) == 0 % rho is constant
   nrho = linspace(rho,rho,nosteps);
else
   if sum(size(rho)>1) > 1
      error('rho must be a row vector')
      return
   end
   if length(rho)~=nosteps
      error('length(rho) must be equal to number of steps')
      return
   end
   rho = rho;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gradx, grady] = oridplot(u)
[gradx, grady] = orideriv(u);
imagesc(sqrt(gradx.^2+grady.^2));
title('Grad. Magnitude')
colorbar('vert')
figure(gcf+1);
ang = 180/pi*atan(grady./(gradx+eps));
imagesc(ang);
colormap(colormapc(4));
title('Grad. Angle')
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cm = Cmcalc(m)
if m <= 1
   error('Use m > 1')
   return
else
   Cm = fzero(strcat('1-exp(-x)-x*exp(-x)*',num2str(m)),[1e-10 1e100]);
end
