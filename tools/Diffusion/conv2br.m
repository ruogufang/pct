function y = conv2br(varargin)
% CONV2BR   2D Convolution with border repetition
%
%    y = CONV2BR(x, kernel) performs the 2D convolution of the image x and
%    the 2D (or 1D) kernel.
%
%    y = CONV2BR(kernel_x, kernel_y, x) performs the 1D convolution
%    of kernel_x and x on the ox direction and of kernel_y and x on 
%    the oy direction, resulting in a 2D convolution.
%    kernel_x and kernel_y must be vectors.
%
%    The image x has its borders repeated prior to the convolutions so
%    that the valid convolution area is the same size as x.
%    kernel_x and kernel_y must be odd vectors.
%
%    y = CONV2BR(kernel_x, kernel_y, x)
%      = conv2(kernel_x, kernel_y, grow(x, .5*[length(kernel_x)-1 length(kernel_y)-1, 'valid']
%

if nargin == 3
   kernel_x = varargin{1};
   kernel_y = varargin{2};
   x = varargin{3};
   
   if mod(length(kernel_x),2)==0 | mod(length(kernel_y),2)==0
      error('kernel_x and kernel_y must be odd sized vectors.');
   end
   
   kern_size = [length(kernel_x) length(kernel_y)];
   kern_size = .5*(kern_size - 1);
   y = conv2(kernel_x, kernel_y, grow(x, kern_size), 'valid');
   
elseif nargin == 2
   x = varargin{1};
   kernel = varargin{2};
   
   if ~any(mod(size(kernel),2))
      error('kernel size in both directoins must be odd.')
   end
   
   kern_size = .5* ( [size(kernel,1) size(kernel,2)] - 1);
   y = conv2(grow(x,kern_size),kernel,'valid');
   
else
   error('Incorrect number of inputs')
end