function [ out ] = pct_interpolate(data,tv,mask)
%PCT_INTERPOLATE Interpolates a dataset with varying sampling intervals
%
%   USAGE:  OUT = PCT_INTERPOLATE(DATA,TV);
%
%   PRE:
%       DATA    - The data to be interpolated [T x Y x X]
%       TV      - The time vector [T x 1] specifying the sampling intervals.
%       MASK    - Optional parameter. A logical mask [Y x X] that indicates
%                 which pixels are to be processed (processes those that are
%                 TRUE).
%
%   POST:
%       OUT     - The interpolated data. [T' x Y x X]
%
%   This function performs simple interpolation for a dataset with two different
%   sampling rates.
%   
%   TODO: Make this function interpolate a dataset with any number of different
%   sampling rates that are integer multiples of the lowest one. 
%
%   WARNING: This function is not very robust (was written in a hurry) and
%   should only be used for datasets with two different contiguous sampling
%   rates where the latter one is twice the former one. I will finish this
%   fuction at a later time.
%
%   Kolbeinn Karlsson 07/07/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if nargin < 3
    mask = true(size(data,2),size(data,3));
end

%Get the time interval vector
sv = zeros(size(tv));
for i = 2:length(tv)
    sv(i) = tv(i)-tv(i-1);
end

%This is because the first element is always zero but will of course belong to
%the first group
sv(1) = sv(2);

%Get the low and high intervals
lowdt = find(sv == min(sv));
highdt = find(sv == max(sv));


%Get data dimension
[l h w] = size(data);

%Pre-allocate output variable
new_length = length(lowdt) + 2*length(highdt);
out = zeros(new_length, h, w);

for i = 1:h
    for j = 1:w
        if mask(i,j)
            temp_vector = squeeze(data(:,i,j));
            first_part = temp_vector(lowdt);
            interp_part = temp_vector(highdt);
            last_part = interp(interp_part,2);
            out(:,i,j) = [first_part; last_part];
        end
    end
end        



end


