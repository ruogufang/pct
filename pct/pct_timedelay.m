function [ aifout ] = pct_timedelay(aif, tau)
%PCT_TIMEDELAY Shifs a time-attenuation curve (TAC) to the right
%
%   USAGE:  AIFOUT = PCT_TIMEDELAY(AIF, TAU)
%
%   PRE:
%       AIF     - A time-attenuation curve [T x 1]
%       TAU     - A time delay constant [Scalar]
%
%   POST:
%       AIFOUT  - AIF shifted to the right by TAU [T x 1]. The first TAU-1
%                 values are padded with ones.
%
%   This function performs a simple right-shift with one-padding on a
%   time-attenuation curve. There is some room for improvement here.
%
%   NOTE: This function is still untested.
%
%   Kolbeinn Karlsson 06/07/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if tau <= 0
    aifout = aif;
else
    a1 = ones([tau 1]);
    a2 = aif(1:end-tau);
    aifout = [a1; a2];
end

end

