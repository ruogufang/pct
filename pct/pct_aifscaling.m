function [ aif_out ] = pct_aifscaling(aif, vof)
%PCT_AIFSCALING Scales an AIF to a reference pixel (the VOF)
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  AIF_OUT = PCT_AIFSCALING(AIF, VOF);
%
%   PRE:
%       AIF     - A time-attenuation curve of a AIF pixel [T x 1]
%       VOF     - A time-attenuation curve of a 100% blood pixel [T x 1]
%
%   POST:
%       AIF_OUT - AIF rescaled so that it's maximum is equal to the
%                 maximum of VOF. [T x 1]
%
%   This function works best if the AIF and the VOF have already
%   been converted from CT units to contrast concentration levels.
%

max_aif = max(aif);
max_vof = max(vof);

if max_aif >= max_vof
    aif_out = aif;
else
    scale_factor = max_vof / max_aif;
    aif_out = aif .* scale_factor;
end

%Truncate any values below zero
i = aif_out < 0;
aif_out(i) = 0;

end

