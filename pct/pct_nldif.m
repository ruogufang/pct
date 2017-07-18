function [ outmap ] = pct_nldif(ctmap)
%PCT_FILTER Spatially filters each frame of a CT series.
%
%   Ruogu Fang 06/21/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUTMAP = PCT_NLDIF(CTMAP, LAMBDA);
%
%   PRE:
%       CTMAP   - A CT Map in HU [T x Y x X]
%
%   POST:
%       OUTMAP  - A filtered CT map in HU [T x Y x X]
%
%   This function spatially filters each frame individually using nonlinear
%   diffusive filters of lambda
%

%Get the image dimensions
[time height width] = size(ctmap);

% Compute lambda
[gx gy] = gradient(squeeze(ctmap(1,:,:)));
g = sqrt(gx.^2 + gy.^2);
lambda = 0.5 * mean(g(:));

%Pre-allocate the ouput variable
outmap = zeros(size(ctmap));

%Filter each frame individually
for i = 1:time
    outmap(i,:,:) = nldif(squeeze(ctmap(i,:,:)),linspace(lambda,lambda,20),linspace(.1,.1,20),10,1000,20,1,2,'aos','dfstep',2,'imscale');
end

close;

end

