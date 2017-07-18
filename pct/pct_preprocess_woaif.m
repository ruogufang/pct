function data = pct_preprocess_woaif(data, PRE, POST)
%PCT_PREPROCESS preprocess the perfuison CTP data
%
%   Ruogu Fang 08/09/2014
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  [DATA, AIF, MASK, MA] = PCT_PREPROCESS(DATA, AIF, VOF);
%
%   PRE:
%       data    - CTP input data [T x Y x X]
%       PRE     - Pre-enhancement cutoff [Scalar]
%       POST    - Post-enhancement cutoff [Scalar]
%
%   POST:
%       data    - Pre-processed CTP data [T x Y x X]
%
%   This function works best if the AIF and the VOF have already
%   been converted from CT units to contrast concentration levels.
%

% PCT Parameters
Ha = 0.45; %Hematocrit correction factor for artery
Hc = 0.25; %Hematocrit correction factor for capillary
loth = 0;     %Lower segmentation threshold
hith = 120;   %Upper segmentation threshold
fsize = 5;    %Size of Gaussian filter
ftsize = 19;    %Size of temporal gaussian filter.
k = 1;          %Contrast conversion factor


% Get the data
data = squeeze(data);
data = double(data);

%Segmentation
data = pct_segment(data, loth, hith, PRE);

% Spatial filtering
% data = pct_filter(data, fsize);

% Time filtering
data = pct_gaussfilter(data,ftsize);

%Subtract base image
data = pct_subtractbase(data, PRE, k);

%Correct for hematocrit differences
data = pct_hematocrit(data,Ha,Hc);

%Truncate the data
data = pct_truncate(data,PRE,POST);

end
