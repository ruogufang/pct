function [varargout] = orieriv(varargin)
%
% ORIDERIV   Calculates the image gradient using Optimized Rotation Invariance kernels.
%
%    [gradx grady] = ORIDERIV(u) calculates the first order derivates of 
%    image u using Optimized Rotational Invariant kernels Fx and Fy.
%
%     Fx = 1/32 * [[ -3 0  3] = Fy'
%                  [-10 0 10]
%                  [ -3 0  3]]
%
%    "u" must have two dimensions.
%
%    Using an additional input parameter, the gradient can be computed only in the specified
%    direction (1 = x, 2 = y)
%
%    Examples: 
%     [gx] = orideriv(u,1)
%     [gy] = orideriv(u,2)

u = varargin{1};

% Verifying inputs
if ndims(u) ~= 2
   error('Input image must have two dimensions.')
   return
end

% ORI Kernel
kernel = 1/32*[3 0 -3; 10 0 -10; 3 0 -3];
%kernel = .5*[1 0 -1];

% ORI Derivation
if nargin == 1 
   varargout{1} = conv2br(u,kernel);  % gx
   varargout{2} = conv2br(u,kernel'); % gy
elseif nargin == 2
   if varargin{2} == 1
      varargout{1} = conv2br(u,kernel);  % gx
   elseif varargin{2} == 2
      varargout{1} = conv2br(u,kernel');  % gy
   else
      error('Unknown parameter.')
   end
else
   error('Too many inputs.')
end


