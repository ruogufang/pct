function [out] = pct_tinterp(data,dt,method)
%PCT_TINTERP Interpolates a dataset with sampling rate dt in the temporal dimensional
%
% USAGE: OUT = PCT_TINTERP(DATA, DT);
%
% PRE:
%       DATA    - The data to be interpolated [T x Y x X]
%       DT      - The upsampling rate [Scaler]
%       METHOD  - Interpolation method. 'linear', 'nearest', 'cubic', or
%                 'spline'. The default method is 'spline'.
%
% POST:
%       OUT     - The interpolated data. [T*DT x Y x X]
%
% This function performs interpolates a dataset in the temporal dimensional
% with constant upsampling rate rt.
%
%
%   Ruogu Fang 5/20/2014
%   Advanced Multimedia Processing (AMP) Lab
%   Department of Electrical and Computer Engineering
%   Cornell University

if nargin < 3
    method = 'spline';
end

[l,h,w] = size(squeeze(data));

% One dimensional temporal data
if ismatrix(data)
    out = interp1(1:l,data,1:1/dt:l,method);
    out = out';
% 3-dimensional spatiao-temporal data
else
    % Interpolation in temporal dimension
    [T,X,Y] = meshgrid(1:h,1:l,1:w);
    [Tq,Xq,Yq] = meshgrid(1:h,1:1/dt:l,1:w);
    out = interp3(T,X,Y,data,Tq,Xq,Yq,'spline');
end

end