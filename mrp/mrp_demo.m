%Ruogu Fang 2014-11-5

clear; close all; clc;

% Set path
addpath('~/Dropbox/Research/MRP/data/MAT');

% Get the data
load('TCGA-06-0133-MRP.mat');
dim = ndims(V);
data = double(squeeze(shiftdim(V,dim-1)));

% MRP Parameters
k_H = 0.73;   %Hematocrit correction factor
loth = 100;     %Lower segmentation threshold
hith = 4500;   %Upper segmentation threshold
rho = 1.05;   %Average brain tissue density
fsize = 5;    %Size of Gaussian filter
ftsize = 10;    %Size of temporal gaussian filter.
TE = info.EchoTime;      % Echo time (ms)
PRE = 18;        %First frame in calculation after delay correction
POST = 95;      %Last frame in calculation after delay correction
dt = 1;         % Time step in time series

% SVD parameters
lambda = 0.15; %Truncation parameter
m = 2;        %Extend the data matrix m time for block circulant

% Compute Brain Mask
B = squeeze(mean(data(1:PRE,:,:,:),1));
mask = pct_brainMask(B,loth,hith);

%Segmentation
data = mrp_segment(data, loth, hith, PRE);

%Convert MR signal intensity to contrast concentration
data = mrp_convert(data);

% Add Gaussian Noise
% standard deviation of Gaussian noise in real and imaginary parts of the
% complex MR signal (assume not vary with time)
sigma = 15;
data_n = mrp_noise(data,sigma,mask);

% Spatial filtering
data = mrp_sfilter(data, fsize);

% Time filtering data = mrp_gaussfilter(data,ftsize);

%Correct for hematocrit differences
data = mrp_hematocrit(data, k_H);

% %Truncate the data data = mrp_truncate(data,PRE,POST);
[T,X,Y,Z]=size(data);

[AIF]=pct_aifsearch(squeeze(data(:,:,:,round(Z/2))),5,[0 10]);
[VOF]=pct_aifsearch(squeeze(data(:,:,:,round(Z/2))),5,[0 10],'Please select the vein');

% %Get AIF and VOF AIF = pct_tac(data, AIFx, AIFy); VOF = pct_tac(data, VOFx,
% VOFy);

%Correct AIF for partial volume effects (PVE)
AIF = pct_aifscaling(AIF,VOF);

%Calculate the residue functions
R = mrp_bsvd(data,AIF,dt,lambda,m,mask);

%Calculate a CBV map
CBV = pct_cbv(R, rho);

%Calculate a CBF map
CBF = pct_cbf(R, rho);

%Calculate a MTT map
% MTT = CBV ./ (CBF + eps) * 60;
MTT = pct_mtt(R, rho);

% Display parameter maps
figure;ctshow(CBF,mask);
figure;ctshow(CBV,mask);
figure;ctshow(MTT,mask);

