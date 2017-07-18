function out = pct_norm(in, lo, hi)
% Normalize value to between 0 and 1
%
%   USAGE:  OUT = PCT_NORM(CTMAP, LO, HI);
%
%   PRE:    
%       IN      - A  [Y x X]
%       LO      - The lower window threshold [Scalar]
%       HI      - The upper window threshold [Scalar]
%       
%   POST:
%       OUT     - An image sequence that has been windowed.
%
% Ruogu Fang 02/22/2013

%Set the output range
lout = 0.0;
hout = 1.0;

%Get series dimensions
[height width] = size(in);
total = height*width;

%Pre-allocate the output
out = zeros([height width], 'double');
in = double(in);

for i = 1:total
    if in(i) <= lo
        out(i) = lout;
    elseif in(i) >= hi
        out(i) = hout;
    else
        out(i) = (in(i) - lo) * (hout-lout) / (hi-lo) + lout;
    end
end

end
