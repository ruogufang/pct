function [varargout] = gsderiv(u, sigma, order, varargin)
%
% GSDERIV   Calculates the derivate of the gaussian smoothed image.
%
%    [out] = GSDERIV( u, sigma, order) calculates the "order"-th derivates of 
%    image u. Prior to derivating the image "u" is smoothed by a sigma 
%    standard deviation 2D gaussian.
%
%    The number of outputs depends on the "order". For "order"=0 there is only one
%    output corresponding to the smoothed input. For "order"=1 there are two outputs,
%    respectivelly the x and y derivates. For "order"=2 there are three outputs,
%    corresponding respectivelly to the 2º x derivate, the crossed x y derivate and
%    the 2º y derivate.
%
%    If order=1, Y = gsderiv(u,sigma,order,'modulus') returns the gradient modulus,
%                [Y1, Y2] = gsderiv(u,sigma,order,'mod_angle') returns the gradient modulus and angle,
%  
%    "u" must have two dimensions.
%
%    Examples:
%       [g] = GSDERIV(x,1,0)
%       [gx gy] = GSDERIV(x,1,1)
%       [gx2 gxy gy2] = GSDERIV(x,1,2)
%

% Verifying inputs
if ndims(u) ~= 2 
    error('Use 2D images.')
end

if all(order ~= [0 1 2])
   error('Order must be 0, 1 or 2.')
   return
end

modulus = 0;
if nargin > 3
    if strcmp(varargin{1},'modulus')
        modulus = 1;
    elseif strcmp(varargin{1},'mod_angle')
        modulus = 2;
    else
        error('Unknown parameter.')
    end
end

% Calculate smoothing kernel
if sigma~=0
   kernel_size = ceil(3*sigma);
   kernel_index = -kernel_size:kernel_size;
   kernel = exp(-.5*(kernel_index./sigma).^2);
   kernel = kernel./sum(kernel);
else % sigma == 0
   kernel = 1;
   kernel_size = 1;
   kernel_index  = 1;
end
   

% Convolution image x kernel (smoothing)
us = conv2br(kernel, kernel, u);

% Calculate image derivates
switch order
case 0
   y = us;
case 1
   y = zeros([size(us) 2]);
   [y(:,:,1), y(:,:,2)] = gradient(us);
   %y(:,:,1) = centdiff(us,2);
   %y(:,:,2) = centdiff(us,1);
case 2
   y = zeros([size(us) 3]);
   [y(:,:,1), y(:,:,3)] = gradient(us); %[dx dy]
   [y(:,:,1), y(:,:,2)] = gradient(y(:,:,1)); % [dx2 dxdy]
   y(:,:,3) = centdiff(y(:,:,3),1); % dydy
   
   %y(:,:,1) = centdiff(us,2);   % dx
   %y(:,:,2) = centdiff(y(:,:,1),1);   % dx dy
   %y(:,:,1) = centdiff(y(:,:,1),2);  % dx2
   %y(:,:,3) = centdiff(us,1);  % dy
   %y(:,:,3) = centdiff(y(:,:,3),1);  % dy2
end 

if order==1           
    if modulus==1
        varargout(1) = { sqrt(y(:,:,1).^2+y(:,:,2).^2)};
    elseif modulus==2
        varargout(1) = {sqrt(y(:,:,1).^2+y(:,:,2).^2)};
        varargout(2) = {180/pi*atan(y(:,:,2)./(y(:,:,1)+eps))};
    else
        varargout(1) = {y(:,:,1)};
        varargout(2) = {y(:,:,2)};
    end
else        
    for i = 1 : size(y,3)
        varargout(i) = {y(:,:,i)};       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = centdiff(x,dim)
% Approximates derivates by central differeces .5*(x(n+1) - x(n-1))
% Differences are calculated on dimension "dim".

p = zeros(1,ndims(x));
p(1,dim) = 1;

y = .5*( roll(x,-p) - roll(x,p) );
y = bordx2(y,dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = bordx2(x,dim)
% Correct Gradient in Borders (duplicates gradient values at image borders)
% dim is the direction in which the gradient was calculated
y = x;
switch dim
case 1
   y([1 end],:,:) = 2*y([1 end],:,:);
case 2
   y(:,[1 end],:) = 2*y(:,[1 end],:);
end
