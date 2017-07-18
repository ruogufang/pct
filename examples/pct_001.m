%Kolbeinn Karlsson 06/11/12
%Advanced Multimedia Processing (AMP) Lab, Cornell University
clear; close all; clc;

addpath(genpath(cd));

% Set path
addpath('../DICOM/IRB');

% PCT Parameters
kappa = 0.73; %Hematocrit correction factor
loth = 0;     %Lower segmentation threshold
hith = 120;   %Upper segmentation threshold
rho = 1.05;   %Average brain tissue density
fsize = 5;    %Size of Gaussian filter
PRE_bbbp = 1; %First frame in BBBP calculation
POST_bbbp = 110; %Last frame in BBBP calculation
sigma = 0;   %Standard deviation of added Gaussian noise to CT data
dt = 0.5;       % Time step in time series
ftsize = 19;    %Size of temporal gaussian filter.

% SVD parameters
lambda = 0.15; %Truncation parameter
m = 2;        %Extend the data matrix m time for block circulant

% Get the data
load('IRB_001.mat');
data = squeeze(V);
data = double(data);

% Compute Brain Mask
B = squeeze(mean(V(1:10,:,:),1));
mask = pct_brainMask(B,0,120);

% Add Gaussian Noise
data = pct_noise(data,[],sigma,'g');
mA = pct_sigma2mA(sigma);

% Spatial filtering
data = pct_filter(data, fsize);

% Non-linear diffusive filtering
% data = pct_nldif(data);

% Time filtering
data = pct_gaussfilter(data,ftsize);

%Segmentation
data = pct_segment(data, loth, hith, PRE);

%Create time vector
% tv = [0:0.5:44 45:74];
% tv = [0:1.25:85 86.5:1.5:115]; % AK data

%Interpolation
% data = pct_interpolate(data,tv);

%Subtract base image
data = pct_contconv(data, PRE);

%Correct for hematocrit differences
data = pct_hematocrit(data, kappa);

% %Auto select AIF
% [AIF_auto,aif_x,aif_y]=pct_aifautoselect(data,mask);
% 
% [AIF_manual]=pct_aifsearch(data,5);

%Get AIF and VOF
AIF = pct_tac(data, AIFy, AIFx);
VOF = pct_tac(data, VOFy, VOFx);

%Correct AIF for partial volume effects (PVE)
AIF = pct_aifscaling(AIF,VOF);

%Truncate the data
data = pct_truncate(data,PRE,POST);
AIF = AIF(PRE:POST);
VOF = VOF(PRE:POST);

%Calculate the residue functions
% R = cTSVD(AIF,data,lambda,m);
R = pct_bsvd(data,AIF,dt,lambda,m,mask);

%Calculate a CBV map
CBV = pct_cbv(R, rho);

%Calculate a CBF map
CBF = pct_cbf(R, rho);

%Multiply by a ratio, Why?
% CBF = CBF * 10;

%Calculate a MTT map
MTT = CBV ./ (CBF + eps) * 60;
% MTT = pct_mtt(R, rho);

% Threshold MTT
MTT(MTT<0)=0;
MTT(MTT>50)=0;

%Remove non-brain tissue
mask = pct_brainMask(CBF);

% Show perfusion maps
figure; ctshow(CBF,mask,[0 80]);
figure; ctshow(CBV,mask,[0 15]);
figure; ctshow(MTT,mask,[0 30]);

save('aif_022.mat','AIF');

%Multiply by a ratio of 2. Why?
% MTT = MTT * 2;

% Calculate BBBP map
% BBBP = pct_bbbp(data,CBV,AIF,dt,rho,PRE_bbbp,POST_bbbp,mask);

% % Save perfusion maps
% dicomwrite(int16(CBF*4),'Output/CBF.dcm');
% dicomwrite(int16(CBV*10),'Output/CBV.dcm');
% dicomwrite(int16(MTT*100),'Output/MTT.dcm');
% dicomwrite(int16(BBBP*100),'Output/BBBP.dcm');
% 
% %Read GE maps
% CBV_GE=double(dicomread('CBV'));
% CBF_GE=double(dicomread('CBF'));
% MTT_GE=double(dicomread('MTT'));
% 
% %Correct GE maps by dividing a ratio
% CBV_GE(CBV_GE==-256) = 0;
% CBV_GE = CBV_GE / 10;
% CBF_GE(CBF_GE==-256) = 0;
% CBF_GE = CBF_GE / 4;
% MTT_GE(MTT_GE==-256) = 0;
% MTT_GE = MTT_GE / 100;