function [data,AIF,VOF] = pct_preprocess_direct(data, PRE, POST, aif_x, aif_y, vof_x, vof_y)
%PCT_PREPROCESS preprocess the perfuison CTP data for direct estimation
%
%   Ruogu Fang 10/21/2015
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  [DATA, AIF, MASK, MA] = PCT_PREPROCESS(DATA, AIF, VOF);
%
%   PRE:
%       data    - CTP input data [T x Y x X]
%       PRE     - Pre-enhancement cutoff [Scalar]
%       POST    - Post-enhancement cutoff [Scalar]
%       AIF_X, AIF_Y - Location of Arterial Input Function [Scalar]
%       VOF_X, VOF_Y - Location of Venous Output Function [Scalar]
%       BASELINE - Baseline image (optional) [Y x X]
%
%   POST:
%       data    - Pre-processed CTP data [T x Y x X]
%       AIF     - Arterial input function [T x 1]
%       VOF     - Venous output function [T x 1]
%


% PCT Parameters
Ha = 0.45; %Hematocrit correction factor for artery
Hc = 0.25; %Hematocrit correction factor for capillary
loth = 0;     %Lower segmentation threshold
hith = 120;   %Upper segmentation threshold
fsize = 5;    %Size of Gaussian filter
ftsize = 19;    %Size of temporal gaussian filter.
k = 1;          %Contrast conversion factor
min_preprocessed = 0;
max_preprocessed = Inf;
min_aif = 1;
max_aif = Inf;

% Get the data
data = squeeze(data);
data = double(data);

% Segmentation
% data = pct_segment(data, loth, hith, PRE);

% if exist('tv','var')
%     %Interpolation
%     data = pct_interpolate(data,tv);
% end

% Spatial filtering
data = pct_filter(data, fsize);

% Time filtering
% data = pct_gaussfilter(data,ftsize);

%Subtract base image
data = pct_subtractbase(data, PRE, k);

%Correct for hematocrit differences
data = pct_hematocrit(data,Ha,Hc);

%Apply min/max values
data = pct_truncatevalues(data,min_preprocessed,max_preprocessed);

%Truncate the data
% data = pct_truncate(data,PRE,POST);
AIF = 0; VOF = 0;

if nargin > 3 && ~isempty(aif_x)
    %Get AIF and VOF
    AIF = pct_tac(data, aif_y, aif_x);
    VOF = pct_tac(data, vof_y, vof_x);
    
    %Correct AIF for partial volume effects (PVE)
    AIF = pct_aifscaling(AIF,VOF);
    AIF = pct_truncatevalues(AIF,min_aif,max_aif);
end

end
