function y = nldifc( u, lambda, sigma, m, stepsize, nosteps, varargin)
%
% NLDIFC  Color Nonlinear Diffusion
%
%    y = NLDIFC( u, lambda, sigma, m, stepsize, steps, verbose, drawstep ) returns the 
%    image y as the result of the application of the Non Linear diffusion
%    to the image u (See NLDIF for details).
%
%    NLDIFC performs an independent diffusion in each color channel of u. If the used
%    parameters are scalars or row vectors, the same parameter will be used for every 
%    channel. To use different parameters for each channel, each parameter must be passed
%    as a 1xCH cell array. String parameters cannot be defined for each channel individually
%
%      Example: y = nldif(u, {3 4 5}, 1, 10, .2, 10) will diffuse ch 1 using lambda = 3
%         ch 2 with lambda = 4 and ch3 with lambda = 5. Every other parameter will be the
%         same for the three channels.
%
%    y = NLDIFC( ..., 'alt1') computes the diffusivity based on the p-norm of the gradient
%    in each channel and filters all channels at once instead of using each channel separately.
%    By default p = 2. For other norms use : y = NLDIFC(..., 'alt1','norm',p). p can be any
%    number or 'inf'.
%
%    See also: NLDIF, PMDIF, EEDIF, CEDIF, DIFFUSIVITY, GSDPLOT, FLUXPLOT
%

Nch = size(u,3);
[verbose, drawstep, imscale, plotgrad, plotflux, dfstep, aos, alt1, pnorm] = parse_inputs(varargin{:});
% Check Inputs
[lambda, sigma, m, stepsize, nosteps, verbose, drawstep, dfstep] = check_inputs(Nch,lambda, sigma, m, stepsize, nosteps, verbose, drawstep, dfstep, alt1);


if alt1 == 0
   
   y = zeros(size(u));
   
   for ch = 1 : size(u,3)
      % MAKEVAR
      v = makevar(verbose{ch},drawstep{ch},imscale,plotgrad,plotflux,dfstep{ch},aos,ch);
      disp( ['Processing channel : ', int2str(ch) , ' / ' , int2str(Nch) ] );
      eval([ 'y(:,:,' , int2str(ch) , ') = nldif(u(:,:, ' , int2str(ch) , ' ),lambda{' , int2str(ch) ,'},sigma{' , int2str(ch) ,'},m{' , int2str(ch) ,'},stepsize{' , int2str(ch) ,'},nosteps{' , int2str(ch) ,'}, ' , v , ');' ]) ;
      
   end
   y = scale(y,[0 1]); % Scaling

   %Plotting
   if any( cat(1,verbose{:}) )
      subplot(1,2,1)
      imagesc(u);
      title('Original Image')
      subplot(1,2,2)
      imagesc(y)
      title('Non Linear Diffusion')
   end
   
else  % alt1 = 1 ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
   % Variable initialization
   dif_time = 0;
   if strcmp(class(u),'double')
      y = u;
   else
      y = double(u);
   end
   g = cell(1,Nch);
   for j = 1 : Nch
      g{j} = zeros( size(y(:,:,1)) );
   end
   
   % Verifying inputs
   if alt1 == 0
      for j = 1 : Nch
         [lambda{j}, sigma{j}, stepsize{j}] = verify_inputs(lambda{j}, sigma{j}, stepsize{j}, nosteps{j});
      end
   else
      for j = 1 : Nch
         [lambda{j}, sigma{j}, stepsize] = verify_inputs(lambda{j}, sigma{j}, stepsize, nosteps);
      end
   end
   
   % Initial Drawing
   if verbose
      figure(verbose);
      subplot(1,2,1);
      if ndims(y) == 3
         image(scale(y,[0 1]));
      else
         if strcmp(imscale,'imscale')
            imagesc(y);
         else
            image(y);
         end
         colorbar
      end
      title('Original Image'); drawnow;
      difplot(y, 0, 0, 'Non Linear Diffusion', verbose, imscale);
   end
   
   % Calculate Cm constant
   for j = 1 : Nch
      Cm{j}= Cmcalc(m{j});
   end
   
   for i=1:nosteps
      
      if mod(i-1,dfstep) == 0 % diffusivity recalc step
         grad = zeros( size(y(:,:,1)) );
         % Calculate gradient of smoothed image
         for j = 1 : Nch    
            [gradx, grady] = gsderiv(y(:,:,j), sigma{j}(i), 1);
            if strcmp(pnorm,'inf')
               grad = max(grad, sqrt(gradx.^2 + grady.^2) );
            elseif pnorm == 2
               grad = grad + gradx.^2 + grady.^2;
            elseif pnorm == 1
               grad = grad + sqrt(gradx.^2 + grady.^2);
            else
               grad = grad + (sqrt(gradx.^2 + grady.^2)).^pnorm;
            end
         end
         
         if ~strcmp(pnorm,'inf') & pnorm~=1
            grad = grad.^(1/pnorm); % Not necessary to do anything if pnorm == 'inf' or 1.
         end
         
         if plotgrad
            figure(verbose+1)
            imagesc(grad);
            title('Grad. Magnitude')
            colorbar
            axis off
         end
         
         % Calculate difusivity
         if samedif(lambda,m,i) % If lambda and m are the same for all channels in this step
            g{1} = 1 - exp(-Cm{1}./( (grad+eps)./lambda{1}(i)).^m{1});
            for j = 2 : Nch
               g{j} = g{1};
            end
         else
            for j = 1 : Nch
               g{j} = 1 - exp(-Cm{j}./( (grad+eps)./lambda{j}(i)).^m{j});  % grad + eps to avoid division by zero (eps -> 0)
            end
         end
         
      end % diff. recalc
      
      % Calculate dy/dt
      if aos
         if plotflux
            yo = y;
         end
         for j = 1 : Nch
            y(:,:,j) = aosiso(y(:,:,j),g{j},stepsize(i)); % updating
         end
      else
         for i = 1 : Nch
            dy(:,:,1) = isodifstep(y(:,:,1), g{j});
         end
         y = y + stepsize(i) * dy;  % updating
      end
      
      % Calculate diffusion time
      dif_time = dif_time + stepsize(i);
      
      
      % Plot actualization   
      if verbose & ~mod(i,drawstep)
         difplot(y, dif_time, i, 'Non Linear Diffusion', verbose, imscale);
         if plotflux
            figure(verbose+plotgrad+1)
            if aos
               fluxplot(g,grad, norm3(abs(y-yo)), 'axisoff')
            else
               fluxplot(g,grad, norm3( abs(stepsize(i).*(dy)) ), 'axisoff')
            end
         end
      end % Plot actualization
      
   end % i = 1 : nosteps 
   % Last plot
   difplot(y, dif_time, i, 'Non Linear Diffusion', verbose, imscale);
   if plotflux
      figure(verbose+plotgrad+1)
      if aos
         fluxplot(g,grad, norm3(abs(y-yo)), 'axisoff')
      else
         fluxplot(g,grad, norm3( abs(stepsize(i).*(dy)) ), 'axisoff')
      end
   end
   
end  % alt1 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = cellcomp(x,n)
for i = 1 : n
   y{i} = x;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [verbose, drawstep, imscale, plotgrad, plotflux, dfstep, aos, alt1, pnorm] = parse_inputs(varargin)
pnorm = 2;
pnormpos = -1;
alt1 = 0;
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
   if i == pnormpos
      flag = 1;
   end
   if strcmp(varargin{i},'hsv')
      hsv = 1;
      flag = 1;
   elseif strcmp(varargin{i},'norm')
      pnorm = varargin{i+1};
      flag = 1;
      pnormpos = i+1;
   elseif strcmp(varargin{i},'alt1')
      alt1 = 1;
      flag = 1;
   elseif strcmp(varargin{i},'imscale')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = quostr(varargin)
y = strcat( char(39),varargin{1},char(39) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nlambda, nsigma, nm, nstepsize, nnosteps, nverbose, ndrawstep, ndfstep] = check_inputs(Nch, lambda, sigma, m, stepsize, nosteps, verbose, drawstep, dfstep, alt1);

if ( length(lambda)==1 & isa(lambda,'cell') ) | ~isa(lambda,'cell')
   nlambda = cellcomp(lambda,Nch);
end
if ( length(sigma)==1 & isa(sigma,'cell') ) | ~isa(sigma,'cell')
   nsigma = cellcomp(sigma,Nch);
end
if ( length(m)==1 & isa(m,'cell') ) | ~isa(m,'cell')
   nm = cellcomp(m,Nch);
end

if alt1 == 0
   if ( length(dfstep)==1 & isa(dfstep,'cell') ) | ~isa(dfstep,'cell')
      ndfstep = cellcomp(dfstep,Nch);
   end
   if ( length(verbose)==1 & isa(verbose,'cell') ) | ~isa(verbose,'cell')
      nverbose = cellcomp(verbose,Nch);
   end
   if ( length(drawstep)==1 & isa(drawstep,'cell') ) | ~isa(drawstep,'cell')
      ndrawstep = cellcomp(drawstep,Nch);
   end
   if ( length(nosteps)==1 & isa(nosteps,'cell') ) | ~isa(nosteps,'cell')
      nnosteps = cellcomp(nosteps,Nch);
   end
   if ( length(stepsize)==1 & isa(stepsize,'cell') ) | ~isa(stepsize,'cell')
      nstepsize = cellcomp(stepsize,Nch);
   end 
else % alt1 == 1
   ndfstep = alt1_check(dfstep, 'dfstep');
   nverbose = alt1_check(verbose ,'verbose');
   ndrawstep = alt1_check(drawstep, 'drawstep');
   nnosteps = alt1_check(nosteps, 'nosteps');
   nstepsize = alt1_check(stepsize, 'stepsize');
end

if length(nlambda) ~= Nch
   error('Incorrect lambda')
end
if length(nsigma) ~= Nch
   error('Incorrect sigma')
end
if length(nm) ~= Nch
   error('Incorrect m')
end
if (length(nstepsize) ~= Nch & alt1 == 0) 
   error('Incorrect stepsize')
end
if (length(nnosteps) ~= Nch & alt1 == 0)  
   error('Incorrect steps')
end
if (length(nverbose) ~= Nch & alt1 == 0)  
   error('Incorrect verbose')
end
if (length(ndrawstep) ~= Nch & alt1 == 0) 
   error('Incorrect drawstep')
end
if (length(ndfstep) ~= Nch & alt1 == 0) 
   error('Incorrect dfstep')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = alt1_check(x,name)
if ( length(x)~=1 & isa(x,'cell') ) 
   error(['When using mode alt1 ', name ,' must be the same for all channels.'])
elseif ( length(x)==1 & isa(x,'cell') )
   if length(x{1})~=1
      error(['When using mode alt1 ' , name  ,' must be the same for all channels.'])
   else
      y = x{1};
   end
%elseif length(x) ~= 1
%   error(['When using mode alt1 ' , name  ,' must be the same for all channels.'])
else
   y = x;
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
function v = makevar(verbose,drawstep,imscale,plotgrad,plotflux,dfstep,aos,ch);
v = ' ';
if verbose ~= 0
   v = strcat(v,',',int2str(verbose));
end
if drawstep ~= 0
   v = strcat(v,',',int2str(drawstep));
end
if strcmp('imscale',imscale)
   v = strcat(v,',imscale');
end
if plotgrad ~= 0
   v = strcat(v,',',quostr('grad'));
end
if plotflux ~= 0
   v = strcat(v,',',quostr('flux'));
end
if dfstep ~= 1
   v = strcat(v,',',quostr('dfstep'),',',int2str(dfstep));
end
if aos ~= 0
   v = strcat(v,',',quostr('aos'));
end
v = strcat(v, ',' , quostr('ch') , ',' , int2str(ch));
v = v(2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = norm3(x,varargin)
if nargin == 2
   p = varargin{1};
else
   p = 2;
end

for i = 1 : size(x,3)
   y = y + x(:,:,i).^p
end
y = y.^(1/p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cm = Cmcalc(m)
if m <= 1
   error('Use m > 1')
   return
else
   Cm = fzero(strcat('1-exp(-x)-x*exp(-x)*',num2str(m)),[1e-10 1e100]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = samedif(lambda,m,i)
y = 1;
for j = 1 : length(m)-1
   if m{j} ~= m{j+1} | lambda{j}(i) ~= lambda{j+1}(i)
      y = 0;
      return
   end
end            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
