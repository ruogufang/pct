function y = colormapc(varargin)
% COLORMAPC   Circular colormap
%
%    y = COLORMAPC(t) creates a circular colormap type t. Circular colormaps
%    are build so that the first and last colors are the same. This kind of
%    colormap is useful for plotting angle images which have a circular behavior
%    (0 = 2*pi).
%    
%    y = COLORMAPC(t,n) creates the circular colormap t with size "n".
%
%    t = 1 -> blk, yel, red, gre, blk        t = 2 -> cya, blu, mag, cya (cool)
%    t = 3 -> blu, cya, gre, yel, blu        t = 4 -> blk, red, yel, blk (hot)
%

switch nargin
case 0
   t = 1;
   n = size(colormap,1);
case 1
   t = varargin{1};
   n = size(colormap,1);
case 2
   t = varargin{1};
   n = varargin{2};
otherwise      
      error('Too many parameters.')
end

y = zeros(n,3);

switch t
case 1 
   y = cmapc(n,[0 0 0],[1 0 0],[1 1 0],[0 1 0]);
case 2 
   y = cmapc(n,[0 1 1],[0 0 1],[1 0 1]);
case 3 
   y = cmapc(n,[0 0 1],[0 1 1],[0 1 0],[1 1 0]);
case 4 
   y = cmapc(n,[0 0 0],[1 0 0],[1 1 0]);
otherwise
   error('Unknown colormap type.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = cmapc(varargin)
if  nargin < 3
   error('Too few inputs.')
end
n = varargin{1};
c = zeros(nargin,3);
for i = 1 : nargin-1
   c(i,:) = varargin{i+1};
end
c(nargin,:) = c(1,:);

y = zeros(n,3);
ls = 1;

for i = 1 : size(c,1)-1
   li = ls;
   ls = round(i*n/(nargin-1));
   d = ls - li + 1;
   y(li:ls,:) = [linspace(c(i,1),c(i+1,1),d); linspace(c(i,2),c(i+1,2),d); linspace(c(i,3),c(i+1,3),d)]';
end
