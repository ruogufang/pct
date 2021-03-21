function [ tac ] = pct_tac(inmap, x, y, z)
%PCT_TAC Extracts a time-attenuation curve (TAC) from a CT map
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   For a given coordinate (x,y), this returns the values of that
%   coordinate for each time index.
%
%   USAGE:  TAC = PCT_TAC(CTMAP, X, Y);
%
%   PRE:
%       CTMAP   - A CT map on the form [T x X x Y x Z]
%       X       - The x coordinate of the AIF pixel
%       Y       - The y coordinate of the AIF pixel
%       Z       - The z coordinate of the AIF pixel (option, default=1)
%
%   POST:
%       TAC     - A time-attenuation curve [T x 1]
%

if nargin < 4
    z = 1;
end

tac = inmap(:,x,y,z);

end

