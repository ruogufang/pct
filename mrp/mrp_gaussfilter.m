function [ outmap ] = mrp_gaussfilter(inmap, l, stddev)
%MRP_GAUSSFILTER Filters each pixel in time in an image sequence
%
%   USAGE:  OUTMAP = MRP_GAUSSFILTER(INMAP);
%
%   PRE:
%       INMAP   - A MR / image sequence [T x X x Y x Z]
%       L       - Length of the filter [Scalar]
%       STDDEV  - The standard deviation of the filter [Scalar]
%
%   POST:
%       OUTMAP  - A filtered image sequence [T x X x Y x Z]
%
%   Filters each pixel individually by convolving it with a windowed Gaussian
%   function.
%
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University

% Get image dimensions
[T,X,Y,Z] = size(inmap);

% Pre-allocate the output variable
outmap = zeros(size(inmap));

if nargin >= 3
    % Create the Gaussian kernell
    alpha = 1/stddev;
    h = gausswin(l,alpha);
else
    h = gausswin(l);
end
h = h/sum(h); %normalize the filter

for i = 1:X
    for j = 1:Y
        for k = 1 : Z
            R = squeeze(inmap(:,i,j,k));
            R = conv(R,h,'same');
            outmap(:,i,j,k) = R;
        end
    end
end

end

