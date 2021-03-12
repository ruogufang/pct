function [ tac ] = pct_tac(inmap, x, y, z, n, r)
%PCT_TAC Extracts a time-attenuation curve (TAC) from a CT map
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   For a given coordinate (x,y), this returns the values of that
%   coordinate for each time index.
%
%   USAGE:  AIF = PCT_TAC(CTMAP, X, Y);
%
%   PRE:
%       CTMAP   - A CT map on the form [T x X x Y x Z]
%       X       - The x coordinate of the AIF pixel
%       Y       - The y coordinate of the AIF pixel
%       Z       - The z coordinate of the AIF pixel (option, default=1)
%       N       - Number of neighborning pixels to average (optional,default=1)
%       R       - Search radius [Scalar] (optional, default=1)
%
%   POST:
%       AIF     - A time-attenuation curve [T x 1]
%

if nargin < 4
    z = 1;
end

if nargin < 5
    n = 1;
end

if nargin < 6
    r = 1;
end

[T,X,Y,Z] = size(inmap);

curves = squeeze(inmap(:,x-r:x+r,y-r:y+r,max(z-r,1):min(z+r,Z)));

%Find the peaks of attenuation curves
peaks = squeeze(max(curves,[],1));

%Find the optimal AIF index with the highest peak
[maxPeak,aifIndex] = sort(peaks(:),'descend');

%Convert index to subscripts
[aif_x, aif_y, aif_z] = ind2sub(size(peaks),aifIndex(1:n));

%Find the optimal AIF
tac = mean(reshape(curves(:,aif_y,aif_x,aif_z),T,[]),2);

end

