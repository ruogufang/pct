function [ out ] = pct_contconv(in, frames)
%PCT_CONTCONV Converts CT values to contrast concentration
%   
%   Kolbeinn Karlsson 06/04/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUT = PCT_CONTCONV(IN, FRAMES);
%
%   PRE:    
%       IN     - A series of CT maps [T x Y x X]
%       FRAMES - No. of frames to average
%
%   POST:
%       OUT    - A series of contrast concentration maps [T x Y x X]
%
%   This function calcualtes the average of each pixel in the first
%   FRAMES frames and subtracts it from every frame
%

%Calculate the average - the CT values before the contrast is introduced
avg = mean(in(1:frames,:,:),1);

%Get input dimensions
[t, Y, X] = size(in);

%Pre-allocate output variable
out = zeros([t,Y,X]);

%Convert to contrast
for i = 1:t
    out(i,:,:) = in(i,:,:) - avg;
end

end

