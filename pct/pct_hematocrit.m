function [ outmap ] = pct_hematocrit(ctmap, Ha, Hc)
%PCT_HEMATOCRIT Corrects a CT map for hematocrit differences.
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUTMAP = PCT_HEMATOCRIT(CTMAP, Ha, Hc);
%
%   PRE:    
%       CTMAP   - A contrast concentration map / CT map [T x Y x X]
%       Ha      - Hematocrit correction factor for artery
%       Hc      - Hematocrit correction factor for capillary
%
%   POST:
%       OUTMAP  - A contrast concentration map corrected for hematocrit
%                 differences
%
%   This function simply multiplies each pixel point-wise with k.

if nargin < 3
    outmap = Ha .* ctmap;
else
    outmap = (1-Ha)/(1-Hc) .* ctmap;
end


end

