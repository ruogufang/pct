function [R, CBF, CBV, MTT, TTP, BBBP] = pct_map(data, AIF, mask, lambda, m)
%PCT_MAP computes the perfusion maps from CTP data
%
%   Ruogu Fang Revised 07/02/2012
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  [R, CBF, CBV, MTT, BBBP] = PCT_MAP(DATA, AIF, MASK);
%
%   PRE:
%       data    - CTP input data [T x Y x X]
%       AIF     - A time-attenuation curve of a AIF pixel [T x 1]
%       mask    - Mask of non-brain tissue [Y x X]
%       lambda  - Truncation parameter for cTSVD [Scalar 0<lambda<1]
%       m       - Extend the data matrix m time for block circulant
%
%   POST:
%       R       - A map of residue functions [T x Y x X]
%       CBF     - A map of regional cerebral blood flow in
%       mL/100g tissue/min [Y x X]
%       CBV     - A map of regional cerberal blood volume in
%                 mL/100g [Y x X]
%       MTT     - A map of Mean Transit Time in seconds [Y x X]
%       TTP     - A map of Time To Peak in seconds [Y x X]
%       BBBP    - A map of blood-brain barrier permeability in
%       mL/100g/min [Y x X]
%
%

% PCT Parameters
rho = 1.05;   %Average brain tissue density
PRE_bbbp = 30; %First frame in BBBP calculation
POST_bbbp = 90; %Last frame in BBBP calculation
POST_bbbp = min(POST_bbbp,size(data,1)); %Last frame for BBBP cannot exceed the actual last frame
dt = 1; % time interval between samples in seconds

% SVD parameters
if nargin < 4
    lambda = 0.15; %Truncation parameter
end
if nargin < 5
    m = 2;        %Extend the data matrix m time for block circulant
end

%Calculate the residue functions
R = pct_bsvd(data,AIF,1,lambda,m,mask);

%Calculate a CBV map
CBV = pct_cbv(R, rho);

%Calculate a CBF map
CBF = pct_cbf(R, rho);

%Multiply by a ratio, Why?
% CBF = CBF * 10;

%Calculate a MTT map
MTT = pct_mtt(R);

%Adjust MTT by dividing 4
% MTT = MTT / 4;

% Calculate TTP map
TTP = pct_ttp(data,dt);

% Calculate BBBP map
if nargout >= 6
    [BBBP, XMAP, YMAP, R2] = pct_bbbp(data,CBV,AIF,dt,rho,PRE_bbbp,POST_bbbp,mask);
end

%Filter non-brain tissue
CBF(~mask) = 0;
CBV(~mask) = 0;
MTT(~mask) = 0;
TTP(~mask) = 0;

%Filter negative tissue
CBF(CBF<0)=0;
CBV(CBV<0)=0;
MTT(MTT<0)=0;
TTP(TTP<0)=0;
