function out = pct_rotavg(in)
%PCT_ROTAVE Computes the rotational average of the input
%
%   USAGE:  OUT = PCT_ROTAVE(IN)
%
%   INPUT:
%       IN      - A [N x N x M] matrix
%
%   OUTPUT:
%       OUT     - Rotational average of the input [N x M]
%
%   This function transforms the input from Cartesian coordinates to Polar
%   coordinates and compute the average over all angles along the radius.
%
%   Ruogu Fang, 10/02/13
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

% Get size of the input
[N,N2,M] = size(in);
% If the size of dimension 1 and 2 differ, warning.
if N ~= N2
    error('The size of the first two dimensions must be same.');
end
% If N is not even, warning.
if mod(N,2)~=0
    error('The size of the first dimension must be even.');
end

[X,Y]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);

[theta,rho]=cart2pol(X,Y);

% Find the pixels with the same radius
rho=round(rho);
idx = cell(N/2+1,1);
for r = 0 : N/2
  idx{r+1} = find(rho==r);
end

out = zeros(N/2+1,M);

% Compute the average values of the pixels on the same radius
for m =1 : M
  image = in(:,:,m);
  for r = 0 : N/2
    out(r+1,m)=mean(image(idx{r+1}));
  end
end

end
