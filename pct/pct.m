%Kolbeinn Karlsson 06/11/12
%Ruogu Fang Revised 06/25/2012
%Advanced Multimedia Processing (AMP) Lab, Cornell University
clear; close all; clc;

% Set path
addpath('tprod');
addpath('Diffusion');
addpath(genpath('~/Dropbox/Research/DICOM'));
outputpath = '~/Dropbox/Research/DICOM/PerfusionMaps';

% PCT Parameters
PRE = 10;     %Pre-enhancement cutoff (first frame included in calc's)
POST = 119;   %Post-enhancement cutoff (last frame included in calc's)
kappa = 0.73; %Hematocrit correction factor
loth = 0;     %Lower segmentation threshold
hith = 120;   %Upper segmentation threshold
rho = 1.05;   %Average brain tissue density
fsize = 5;    %Size of Gaussian filter
PRE_bbbp = 30; %First frame in BBBP calculation
POST_bbbp = 90; %Last frame in BBBP calculation
sigma = 20;   %Standard deviation of added Gaussian noise to CT data

% SVD parameters
lambda = 0.3; %Truncation parameter
m = 2;        %Extend the data matrix m time for block circulant

% Reference artery and veins
AIFx = 246;
AIFy = 157;
VOFx = 233;
VOFy = 432;

% Get the data
load 'CTP_003.mat'
data = squeeze(V);
data = double(data);

% Add Gaussian Noise
% data = pct_noise(data,[],sigma);
% mA = pct_sigma2mA(sigma);

% Make the indexing right
data = permute(data, [3 1 2]);

% Spatial filtering
data = pct_filter(data, fsize);

% Find brain mask
I = squeeze(data(1,:,:));
mask = (I < loth | I > hith);

%Segmentation
data = pct_segment(data, loth, hith, PRE);

%Create time vector
tv = [0:0.5:44 45:74];
% tv = [0:1.25:85 86.5:1.5:115]; % AK data

%Interpolation
data = pct_interpolate(data,tv);

% Non-linear diffusive filtering
% data = pct_nldif(data);

%Subtract base image
data = pct_contconv(data, PRE);

%Get AIF and VOF
AIF = pct_tac(data, AIFx, AIFy);
VOF = pct_tac(data, VOFx, VOFy);

%Correct AIF for partial volume effects (PVE)
AIF = pct_aifscaling(AIF,VOF);

%Correct for hematocrit differences
data = pct_hematocrit(data, kappa);

%Truncate the data
data = pct_truncate(data,PRE,POST);
AIF = AIF(PRE:POST);

%Calculate the residue functions
R = cTSVD(AIF, data, lambda, m);

%Calculate a CBV map
CBV = pct_cbv(R, rho);

%Calculate a CBF map
CBF = pct_cbf(R, rho);

%Multiply by a ratio, Why?
CBF = CBF * 10;

%Calculate a MTT map
MTT = pct_mtt(R);

%Adjust MTT by dividing 4
MTT = MTT / 4;

%Filter non-brain tissue
MTT(mask) = 0;

% Calculate BBBP map
BBBP = pct_bbbp(data,CBV,AIF,rho,PRE_bbbp,POST_bbbp);

% Save perfusion maps
dicomwrite(int16(CBF),'Output/CBF.dcm');
dicomwrite(int16(CBV),'Output/CBV.dcm');
dicomwrite(int16(MTT),'Output/MTT.dcm');
dicomwrite(int16(BBBP),'Output/BBBP.dcm');

%Read GE maps
CBV_GE=double(dicomread('CBV'));
CBF_GE=double(dicomread('CBF'));
MTT_GE=double(dicomread('MTT'));

%Correct GE maps by dividing a ratio
CBV_GE(CBV_GE==-256) = 0;
CBV_GE = CBV_GE / 10;
CBF_GE(CBF_GE==-256) = 0;
CBF_GE = CBF_GE / 4;
MTT_GE(MTT_GE==-256) = 0;
MTT_GE = MTT_GE / 100;
