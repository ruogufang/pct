%Kolbeinn Karlsson 06/11/12
%Advanced Multimedia Processing (AMP) Lab, Cornell University
clear; close all; clc;

%Time the algorithm
ticID = tic;

%%
% PCT Parameters
PRE = 10;       %Pre-enhancement cutoff (first frame included in calc's)
POST = 145; %89    %Post-enhancement cutoff (last frame included in calc's)
kappa = 0.73;   %Hematocrit correction factor
loth = 0;       %Lower segmentation threshold
hith = 120;     %Upper segmentation threshold
rho = 1.05;     %Average brain tissue density
fsize = 3;      %Size of Gaussian filter (spatial filtering)
sigma = 0.9;    %Standard deviation of Gaussian filter (spatial filtering)
k = 1;          %Contrast conversion factor
ftsize = 5;    %Size of temporal gaussian filter.
dt = 1; %0.5;       %The time interval between CT samples

% SVD parameters
lambda = 0.3;   %Truncation parameter
m = 2;          %Extend the data matrix m time for block circulant

% Min and max values
min_cbv = 0;
max_cbv = Inf;
min_cbf = 0;
max_cbf = Inf;
min_mtt = 0;
max_mtt = Inf;
min_ttp = 0;
max_ttp = Inf;
min_preprocessed = 0;
max_preprocessed = Inf;
min_bbbp = 0;
max_bbbp = Inf;
min_aif = 1;
max_aif = Inf;
max_delay = 15; % average circulation time 6-11 sec

first_bbbp = 73;
last_bbbp = 89;

I0 = 190; % High-dose tube current level
I = 50; % Simulated low-dose tube current level

% y1 = 200; y2 = 202; x1 = 200; x2 = 202;

% Get the data
% load 'Data/IRB_001_S3.mat'
load '../DICOM/IRB/IRB_017.mat'
%This file contains the following variables:
% V  - containing the image data [T x Y x X]
% aif_x - the x coordinate of the AIF
% aif_y - the y coordinate of the AIF
% vof_x - the x coordinate of the VOF
% vof_y - the y coordinate of the VOF

data = V;
aif_x = AIFx; aif_y = AIFy; vof_x = VOFx; vof_y = VOFy;

%Get length
[len,h,w] = size(data);

%Simulate low-dose
sigma = pct_mA2sigma(I,I0);
% data = pct_noise(data,[],sigma);

%Registration
% data = pct_registration(data);

%Segmentation
data = pct_segment(data, loth, hith, PRE);

%Create time vector
tv = [0:0.5:44 45:74];
% tv = [0:1.25:85 86.5:1.5:115]; % AK data
% tv = 0:1:88;

% %Interpolation
% data = pct_interpolate(data,tv);

% Spatial filtering
data = pct_filter(data, fsize, sigma);

% Time filtering
data = pct_gaussfilter(data,ftsize);

% Filtering
%data = pct_nldif(data);

%Subtract base image
data = pct_subtractbase(data, PRE, k);

%Get AIF and VOF
AIF = pct_tac(data, aif_x, aif_y);
VOF = pct_tac(data, vof_x, vof_y);

%Correct AIF for partial volume effects (PVE)
AIF = pct_aifscaling(AIF,VOF);

%Correct for hematocrit differences
data = pct_hematocrit(data, kappa);

%Apply min/max values
data = pct_truncatevalues(data,min_preprocessed,max_preprocessed);
AIF = pct_truncatevalues(AIF,min_aif,max_aif);

%Truncate the data
data_t = pct_truncate(data,PRE,POST);
AIF_t = AIF(PRE:POST);

%Save "data" to disk to use less memory
save data.mat data
clear data

% % Select ROI
% data_t = data_t(:,y1:y2,x1:x2);
% mask = mask(y1:y2,x1:x2);

%Calculate the residue functions
% R = cTSVD(AIF_t, data_t, lambda, m);
R = pct_bsvd(data_t, AIF_t, dt, lambda, m, mask);

%Now get rid of the truncated data and get the full data back from disk
clear data_t AIF_t
load data.mat

%%

%Apply the time correction
aifpeak = find(AIF == max(AIF));
R = R*(1/dt);
CBV = pct_cbv(R, rho, dt);
CBF = pct_cbf(R, rho, dt);
MTT = pct_mtt(R, dt);
TTP = pct_ttp(data, dt, aifpeak, max_delay);
CBV = pct_truncatevalues(CBV,min_cbv,max_cbv);
CBF = pct_truncatevalues(CBF,min_cbf,max_cbf);
MTT = pct_truncatevalues(MTT,min_mtt,max_mtt);


CBV = mask .* CBV;
CBF = mask .* CBF;
MTT = mask .* MTT;
TTP = mask .* TTP;


%%

%Compute the BBBP map without delay correction
% [BBBP X YMAP R] = pct_bbbp(data,CBV,AIF,dt,rho,first_bbbp,last_bbbp,mask);
 
% TTP = ones(size(CBV))*aifpeak*dt;

% Compute BBBP with delay correction
[BBBP XMAP YMAP R] = pct_bbbpdc(data,CBV,TTP,AIF,dt,rho,...
   first_bbbp,last_bbbp,mask);

%Apply min/max values to the BBBP
BBBP = pct_truncatevalues(BBBP, min_bbbp, max_bbbp);
BBBP = mask .* BBBP;

%%

%Clear intermediate variables
clear PRE POST kappa loth hith fsize sigma ...
    lambda m aif_x aif_y vof_x vof_y time min_aif ...
    max_aif min_cbv max_cbv min_cbf max_cbf ...
    min_preprocessed max_preprocessed min_mtt ...
    max_mtt min_ttp max_ttp min_bbbp max_bbbp ...
    k h w first_bbbp ftsize last_bbbp i j R2 len tv data ...
    aifpeak max_delay I I0

clear R rho VOF AIF dt XMAP YMAP
    

%Get time performance
toc(ticID);
clear ticID

