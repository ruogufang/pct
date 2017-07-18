function y = grow(varargin)
%
% GROW    Expands a matrix.
%
%   Y = GROW(X,[M1 M2 ... Mn]) creates a matrix 'y' by the expansion
%   of the matrix 'x' of Mi elements on directions i e -i. The size 
%   of the new matrix 'y' will be:
%       size(y) = size(x) + 2*[M]
%
%   Example: x = [ 1 2 3; 4 5 6; 7 8 9]
%   y = GROW(x,[1 0])
%   y = [1 2 3; 1 2 3; 4 5 6; 7 8 9; 7 8 9]
%
%   y = GROW(x,[1 1])
%   y = [1 1 2 3 3; 1 1 2 3 3; 4 4 5 6 6; 7 7 8 9 9; 7 7 8 9 9]
%
%   y = GROW(x) uses as standard Mi = 1.
%
%    If DIM(x)>1, GROW(x,N) = GROW(x,N*ones(DIM(x)))
%
%   y = GROW(x,[M], 'a') performs assimetric expansion. If Mi is even
%   assimetric expansion is the same as simmetric expansion, execpt that
%   the expansion in the assimetric mode will be half of the one in the
%   simmetric mode. If Mi is odd, the matrix will be expanded ceil(Mi/2)
%   elements in the -i direction and floor(Mi/2) in the i. The size of
%   the expanded matrix will be:  
%       size(y) = size(x) + [M]
%
%    See also ROLL, PAD, PADC.


[X, M, assimetric] = parse_inputs(varargin{:});

ind = cell(length(M),1);

if assimetric == 0
   if length(M) == 1
      ind{1} = [repmat(1,[1 M]) 1:length(X) repmat(length(X),[1 M]) ];
   else
      for i = 1 : length(M)
         ind{i} = [repmat(1,[1 M(i)]) 1:size(X,i) repmat(size(X,i),[1 M(i)]) ];
      end
   end
else
   if length(M) == 1
      ind{1} = [repmat(1,[1 ceil(M/2)]) 1:length(X) repmat(length(X),[1 floor(M/2)]) ];
   else
      for i = 1 : length(M)
         ind{i} = [repmat(1,[1 ceil(M(i)/2)]) 1:size(X,i) repmat(size(X,i),[1 floor(M(i)/2)]) ];
      end
   end
end
   
y = X(ind{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,M,assimetric] = parse_inputs(varargin)
assimetric = 0;

switch nargin
case 0
   error('Too few inputs!')
   return
   
case 1
   X = varargin{1};
   M = ones(1,ndims(X));
   
case 2
   X = varargin{1};
   M = varargin{2};
   if strcmp(M,'a')
      M = ones(1,ndims(X));
      assimetric = 1;
   end
   
case 3
   X = varargin{1};
   M = varargin{2};
   if ~strcmp('a', lower(varargin{3}))
      error('Unknown parameter.');
   end
   assimetric = 1;
  
case 4
   error('Too many inputs!')
   return
end

if length(M)==1 & sum(size(X)>1)>1
   M = M*ones(1,ndims(X));
end

if length(M)~=sum(size(X)>1)>1
   error('Invalid dimensions')
end

