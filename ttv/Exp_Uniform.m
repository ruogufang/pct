% Simulation Experiment on Uniform CBF region estimation
% Figure 4
%   Ruogu Fang Revised 09/24/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

clear; clc; close all;

% Set path
addpath(genpath('toolbox')); % Include toolbox package

% Parameters
PSNR = 15;
reg = 100; % TV weight
lambda = 0.15; % bSVD truncation threshold
m = 2; % circulant number
CBV = 4; %ml/100g
CBF = 15; % 20 : 10 : 80 ml/100g/min
X = 40; Y = 40; Z = 1; T = 60; % volumn size

% Step 1: Generate Artery Input Function (AIF)
t = (0 : T-1)';
t0 = 0;
a = 1;
b = 3;
c = 1.5;
AIF = pct_aif(t,t0,a,b,c);

% Step 2: Generate Impulse Residue Function RIF using Gamma Family
CBF_gt = repmat(CBF,X,Y);
MTT = CBV / CBF * 60; %sec (3 ~ 12 s)
MTT_gt = repmat(MTT,[X Y]);
CBV_gt = repmat(CBV,[X Y]);
alpha = 10;
beta = MTT / alpha;
RIF = pct_irf(t, alpha, beta);
% RIF = exp(-(t-t0)/MTT); % exponential RIF

% Step 3: Generate Tissue Time Enhancement Curve
TEC = pct_tec(AIF, RIF, CBF);

% Step 4: Correlated Gaussian Noise Generation
load acf;
sigma = max(TEC(:))/(10^(PSNR/20));
Cv = repmat(TEC,[1 X Y Z]);
Cn = pct_noise(Cv, acf, sigma,'s');
TTP_gt = pct_ttp(Cv,1);

% Method 1: Standard Truncated Singular Value Decomposition (sSVD)
RIF_sSVD = pct_ssvd(Cn,AIF,1,lambda,ones(X,Y));
CBF_sSVD = squeeze(max(RIF_sSVD));
MTT_sSVD = squeeze(sum(RIF_sSVD))./CBF_sSVD;
CBV_sSVD = squeeze(sum(RIF_sSVD))/60;
TTP_sSVD = pct_ttp(pct_tec(AIF,RIF_sSVD),1);

% Method 2: block-circulant SVD (bSVD)
RIF_bSVD = pct_bsvd(Cn,AIF,1,lambda,m,ones(X,Y));
CBF_bSVD = squeeze(max(RIF_bSVD));
MTT_bSVD = squeeze(sum(RIF_bSVD))./CBF_bSVD;
CBV_bSVD = squeeze(sum(RIF_bSVD))/60;
TTP_bSVD = pct_ttp(pct_tec(AIF,RIF_bSVD),1);

% Method 2-1: block-circulant SVD (bSVD) no noise
RIF0_bSVD = pct_bsvd(Cv,AIF,1,lambda,m,ones(X,Y));
CBF0_bSVD = squeeze(max(RIF0_bSVD));
MTT0_bSVD = squeeze(sum(RIF0_bSVD))./CBF0_bSVD;
CBV0_bSVD = squeeze(sum(RIF0_bSVD))/60;
TTP0_bSVD = pct_ttp(pct_tec(AIF,RIF0_bSVD),1);

% Method 3: Tikhonov
RIF_tikh = pct_tikh(Cn,AIF,1,lambda,1,ones(X,Y));
CBF_tikh = squeeze(max(RIF_tikh));
MTT_tikh = squeeze(sum(RIF_tikh))./CBF_tikh;
CBV_tikh = squeeze(sum(RIF_tikh))/60;
TTP_tikh = pct_ttp(pct_tec(AIF,RIF_tikh),1);

% Method 4: Sparse Perfusion Deconvolution (SPD)
params.cbf0 = CBF_bSVD;
params.mtt0 = MTT_bSVD;
params.cbv0 = CBV_bSVD;
params.ttp0 = TTP_bSVD;
params.numThreads=2; % number of threads
params.K = 256; % dictionary size
params.sigma = sigma;
params.maxval = max(CBF_bSVD(:));
params.trainnum = 40000;
params.verbose = false;
params.lambda = 10; % weight of the sparsity term in SPAMS package
params.mixture = 1; % weight of input signal
params.J = []; % initialize energy function
params.round = 2; % MAP rounds
params.C = Cn;
params.AIF = AIF;
params.R = RIF_bSVD;
params.verbose = false;
params.rho = 1.05;
params.m = m;
params.beta = 10; % weight of prior
params.blocksize = 6;
overlap = 5; % patch overlap
params.stepsize = params.blocksize - overlap;

% %Load pre-trained dictionary
load D_ODL_6x6.mat
params.D = dict;

% Online deconvolution
[CBF_SPD, MTT_SPD, CBV_SPD, TTP_SPD] = spd(params);

% Method 5: Time-Intensity Profile Similarity (TIPS) bilateral filtering +bSVD 
w     = 5;       % bilateral filter half-width
sigma_tips = [5 sigma]; % bilateral filter standard deviations
Cn_TIPS = permute(pct_tips(permute(Cn,[2 3 1]),w,sigma_tips),[3 1 2]);
RIF_TIPS_bSVD = pct_bsvd(Cn_TIPS,AIF,1,lambda,m,ones(X,Y));
CBF_TIPS_bSVD = squeeze(max(RIF_TIPS_bSVD));
MTT_TIPS_bSVD = squeeze(sum(RIF_TIPS_bSVD))./CBF_TIPS_bSVD;
CBV_TIPS_bSVD = squeeze(sum(RIF_TIPS_bSVD))/60;
TTP_TIPS_bSVD = pct_ttp(pct_tec(AIF,RIF_TIPS_bSVD),1);

% Method 6: Tensor Total Variation on RIF (TTV)
% Parameters for TV-based Methods
input.reg = reg;
input.n1=X;input.n2=Y; input.nt = T;
input.maxitr=500;
input.L=1;input.num=1; 
input.l=0; input.u=80;
input.no = 50;
Ca = pct_circ(AIF,[],[],1); % block-circulant version
input.A = Ca;
input.b = reshape(Cn,T,[]);

out = TTV_FCSA(input);
RIF_TTV = reshape(out.y,[T,X,Y]);
CBF_TTV = squeeze(max(RIF_TTV));
MTT_TTV = squeeze(sum(RIF_TTV))./CBF_TTV;
CBV_TTV = squeeze(sum(RIF_TTV))/60;
TTP_TTV = pct_ttp(pct_tec(AIF,RIF_TTV),1);

% Method 7: TIPS+TTV
input.b = reshape(Cn_TIPS,T,[]);
out_tipsttv = TTV_FCSA(input);
RIF_TIPS_TTV = reshape(out_tipsttv.y,[T,X,Y]);
RIF_TIPS_TTV = RIF_TIPS_TTV(1:T,:,:); % keep first T data points
CBF_TIPS_TTV = squeeze(max(RIF_TIPS_TTV));
MTT_TIPS_TTV = squeeze(sum(RIF_TIPS_TTV))./CBF_TIPS_TTV;
CBV_TIPS_TTV = squeeze(sum(RIF_TIPS_TTV))/60;
TTP_TIPS_TTV = pct_ttp(pct_tec(AIF,RIF_TIPS_TTV),1);

% Method 8: Initialization with baseline
input.b = reshape(Cn,T,[]);
input.yr = reshape(RIF_bSVD,T,[]);
out_ttv_init = TTV_FCSA_init(input);
RIF_TTV_init = reshape(out_ttv_init.y,[T,X,Y]);
RIF_TTV_init(RIF_TTV_init<0) = 0;
CBF_TTV_init = squeeze(max(RIF_TTV_init));
MTT_TTV_init = squeeze(sum(RIF_TTV_init))./CBF_TTV_init;
CBV_TTV_init = squeeze(sum(RIF_TTV_init))/60;
TTP_TTV_init = pct_ttp(pct_tec(AIF,RIF_TTV_init),1);

%% Show results
fs = 20;
figure; 
subplot(491); imshow(CBF_gt,[0 50]); set(gca,'FontSize',fs); title('Reference'); 
text(-20,20,'CBF','FontSize',fs);
subplot(492); imshow(CBF_sSVD,[0 50]); set(gca,'FontSize',fs); title('sSVD'); 
subplot(493); imshow(CBF_bSVD,[0 50]); set(gca,'FontSize',fs); title('bSVD'); 
subplot(494); imshow(CBF_TIPS_bSVD,[0 50]); set(gca,'FontSize',fs); title('TIPS+bSVD');
subplot(495); imshow(CBF_tikh,[0 50]); set(gca,'FontSize',fs); title('Tikhonov'); 
subplot(496); imshow(CBF_SPD,[0 50]); set(gca,'FontSize',fs); title('SPD');
subplot(497); imshow(CBF_TTV,[0 50]); set(gca,'FontSize',fs); title('TTV');
subplot(498); imshow(CBF_TIPS_TTV,[0 50]); set(gca,'FontSize',fs); title('TIPS+TTV');
subplot(499); imshow([],[0 50]);
c=colorbar; set(gca,'FontSize',fs);
x1=get(gca,'position');
x=get(c,'Position');
x(1)=x1(1);x(2)=x1(2);x(3)=0.015;x(4)=0.15;
set(c,'Position',x)
set(gca,'position',x1)

subplot(4,9,10); imshow(CBV_gt,[0 6]); 
text(-20,20,'CBV','FontSize',20);
subplot(4,9,11); imshow(CBV_sSVD,[0 6]);
subplot(4,9,12); imshow(CBV_bSVD,[0 6]);
subplot(4,9,13); imshow(CBV_TIPS_bSVD,[0 6]); 
subplot(4,9,14); imshow(CBV_tikh,[0 6]); 
subplot(4,9,15); imshow(CBV_SPD,[0 6]); 
subplot(4,9,16); imshow(CBV_TIPS_TTV,[0 6]); 
subplot(4,9,17); imshow(CBV_TTV,[0 6]); 
subplot(4,9,18); imshow([],[0 6]);
c=colorbar; set(gca,'FontSize',fs);
x1=get(gca,'position');
x=get(c,'Position');
x(1)=x1(1);x(2)=x1(2);x(3)=0.015;x(4)=0.15;
set(c,'Position',x);
set(gca,'position',x1);

subplot(4,9,19); imshow(MTT_gt,[0 20]);
text(-20,20,'MTT','FontSize',20);
subplot(4,9,20); imshow(MTT_sSVD,[0 20]); 
subplot(4,9,21); imshow(MTT_bSVD,[0 20]);
subplot(4,9,22); imshow(MTT_TIPS_bSVD,[0 20]);
subplot(4,9,23); imshow(MTT_tikh,[0 20]); 
subplot(4,9,24); imshow(MTT_SPD,[0 20]);
subplot(4,9,25); imshow(MTT_TTV,[0 20]);
subplot(4,9,26); imshow(MTT_TIPS_TTV,[0 20]);
subplot(4,9,27); imshow([],[0 20]);
c=colorbar; set(gca,'FontSize',fs);
x1=get(gca,'position');
x=get(c,'Position');
x(1)=x1(1);x(2)=x1(2);x(3)=0.015;x(4)=0.15;
set(c,'Position',x);
set(gca,'position',x1);

ttp_range = [0 20];
subplot(4,9,28); imshow(TTP_gt,ttp_range);
text(-20,20,'TTP','FontSize',20);
subplot(4,9,29); imshow(TTP_sSVD,ttp_range); 
subplot(4,9,30); imshow(TTP_bSVD,ttp_range);
subplot(4,9,31); imshow(TTP_TIPS_bSVD,ttp_range);
subplot(4,9,32); imshow(TTP_tikh,ttp_range); 
subplot(4,9,33); imshow(TTP_SPD,ttp_range);
subplot(4,9,34); imshow(TTP_TTV,ttp_range);
subplot(4,9,35); imshow(TTP_TIPS_TTV,ttp_range);
subplot(4,9,36); imshow([],ttp_range);
c=colorbar; set(gca,'FontSize',fs);
x1=get(gca,'position');
x=get(c,'Position');
x(1)=x1(1);x(2)=x1(2);x(3)=0.015;x(4)=0.15;
set(c,'Position',x);
set(gca,'position',x1);

% RMSE
rmse_CBF_bSVD = pct_rmse(CBF_bSVD,CBF_gt);
rmse_CBF_TTV = pct_rmse(CBF_TTV,CBF_gt);

rmse_MTT_bSVD = pct_rmse(MTT_bSVD,MTT_gt);
rmse_MTT_TTV = pct_rmse(MTT_TTV,MTT_gt);

% Save figures
figpath = '../../../Figures/pdf';
w = 25; h = 10;
set(gcf, 'papersize', [w h]);
set(gcf, 'paperposition', [0 0 w h]);
print(gcf,fullfile(figpath,'uniform'),'-dpdf');
saveas(gcf,'../../../Figures/fig/uniform.fig');

%% Answer to R3 question on initialization
figure;plot(out.funv,'LineWidth',2);
hold on; plot(out_ttv_init.funv,'^r','LineWidth',2);
plot(out_tipsttv.funv,'-.g','LineWidth',2);
set(gca,'FontSize',20);
legend('TTV','TTV init with bSVD','TIPS+TTV');
