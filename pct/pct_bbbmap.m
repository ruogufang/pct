function [BBBP, X, YMAP] = pct_bbbmap(data, AIF, mask)
%PCT_BBBMAP computes the BBBP map from CTP data
%
%   Ruogu Fang Revised 03/013/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  [BBBP] = PCT_MAP(DATA, AIF, MASK);
%
%   PRE:
%       data    - CTP input data [T x Y x X]
%       AIF     - A time-attenuation curve of a AIF pixel [T x 1]
%       mask    - Mask of non-brain tissue [Y x X]
%
%   POST:
%       BBBP    - A map of blood-brain barrier permeability in
%       mL/100g/min [Y x X]
%       X       - The independent variable of the patlak plot [T x 1]
%       YMAP    - A map of dependent variables of the patlak plot [T x Y x X]
%       R *      - A map of coefficients of determination (R^2) 
%
%

% PCT Parameters
rho = 1.05;   %Average brain tissue density
PRE_bbbp = 30; %First frame in BBBP calculation
POST_bbbp = 90; %Last frame in BBBP calculation
POST_bbbp = min(POST_bbbp,size(data,1)); %Last frame for BBBP cannot exceed the actual last frame
dt = 1; % time interval between samples in seconds

% SVD parameters
lambda = 0.15; %Truncation parameter
m = 2;        %Extend the data matrix m time for block circulant

%Calculate the residue functions
K = pct_bsvd(data,AIF,1,lambda,m,mask);

%Calculate a CBV map
CBV = pct_cbv(K, rho);

% Calculate BBBP map
[BBBP, X, YMAP] = pct_bbbp(data,CBV,AIF,dt,rho,PRE_bbbp,POST_bbbp,mask);

