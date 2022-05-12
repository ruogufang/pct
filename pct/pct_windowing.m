function [ out ] = pct_windowing(ctseries, lo, hi )
%PCT_WINDOWING Performs image contrast windowing on an image sequence
%
%   Windows the image contrast according to input paramters. The window 
%   width is (HI-LO) and the window center is (HI+LO)/2.
%
%   USAGE:  OUT = PCT_WINOWING(CTMAP, LO, HI);
%
%   PRE:    
%       CTMAP   - A  [Y x X x T]
%       LO      - The lower window threshold [Scalar]
%       HI      - The upper window threshold [Scalar]
%       
%   POST:
%       OUT     - An image sequence that has been windowed.
%
%   Note that this function takes input on the form [Y x X x T] and not 
%   [T x Y x X] like most other pct_ functions.
%
%   Kolbeinn Karlsson 06/05/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Set the output range
lout = 0.0;
hout = 1.0;

%Get series dimensions
[time height width] = size(ctseries);
total = time*height*width;

%Pre-allocate the output
out = zeros([time height width], 'double');
ctseries = double(ctseries);

for i = 1:total
    if ctseries(i) <= lo
        out(i) = lout;
    elseif ctseries(i) >= hi
        out(i) = hout;
    else
        out(i) = (ctseries(i) - lo) * (hout-lout) / (hi-lo) + lout;
    end
end

end

