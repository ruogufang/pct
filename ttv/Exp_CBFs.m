% Simulation Experiment on CBF estimation at various CBF ground truth
% values
% Fig 3 (a)-(b) & Fig 5 (a)-(b)
%   Ruogu Fang Revised 09/24/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

clear; clc; close all; warning off;

% Set path
addpath(genpath('toolbox')); % Include toolbox package

% Parameters
PSNR = 15; % about 5mA
reg = 100; % TV weight
lambda = 0.15; % bSVD truncation threshold
m = 2; % circulant number
CBV = 4; %ml/100g
CBFs = (10 : 10 : 80)'; % ml/100g/min
N = length(CBFs);
MTTs = CBV./CBFs * 60;
CBVs = repmat(CBV,N,1);
X = 40; Y = 40; Z = 1; T = 60; % volumn size

% Step 0: Initilization
mean_sSVD = zeros(N,4);
std_sSVD = zeros(N,4);
mean_bSVD = zeros(N,4);
std_bSVD = zeros(N,4);
mean_TIPS_bSVD = zeros(N,4);
std_TIPS_bSVD = zeros(N,4);
mean_tikh = zeros(N,4);
std_tikh = zeros(N,4);
mean_SPD = zeros(N,4);
std_SPD = zeros(N,4);
mean_TTV = zeros(N,4);
std_TTV = zeros(N,4);
mean_TIPS_TTV = zeros(N,4);
std_TIPS_TTV = zeros(N,4);
TTPs = zeros(size(CBFs));

% Parameters for TV-based Methods
input.reg = reg;
input.n1=X;input.n2=Y; input.nt = T;
input.maxitr=500;
input.L=1;input.num=1; 
input.l=0; input.u=inf;
input.no = 50;
alpha = 10;
load acf;

% Step 1: Generate Artery Input Function (AIF)
t = (0 : T-1)';
t0 = 0;
a = 1;
b = 3;
c = 1.5;
AIF = pct_aif(t,t0,a,b,c);

for n = 1 : N
    CBF = CBFs(n);

% Step 2: Generate Impulse Residue Function RIF using Gamma Family
CBF_gt = repmat(CBF,X,Y);
MTT = CBV / CBF * 60; %sec (3 ~ 12 s)
beta = MTT / alpha;
RIF = pct_irf(t, alpha, beta);

% Step 3: Generate Tissue Time Enhancement Curve
TEC = pct_tec(AIF, RIF, CBF);

% Step 4: Correlated Gaussian Noise Generation
% sigma = max(TEC(:))/PSNR;
sigma = max(TEC(:))/(10^(PSNR/20));
Cv = repmat(TEC,[1 X Y Z]);
Cn = pct_noise(Cv, acf, sigma,'s');
TTPs(n) = pct_ttp(TEC,1);

% Method 1: Standard Truncated Singular Value Decomposition (sSVD)
RIF_sSVD = pct_ssvd(Cn,AIF,1,lambda,ones(X,Y));
CBF_sSVD = squeeze(max(RIF_sSVD));
MTT_sSVD = squeeze(sum(RIF_sSVD))./CBF_sSVD;
CBV_sSVD = squeeze(sum(RIF_sSVD))/60;
TTP_sSVD = pct_ttp(pct_tec(AIF,RIF_sSVD),1);
mean_sSVD(n,1) = mean(CBF_sSVD(:));
std_sSVD(n,1) = std(CBF_sSVD(:));
mean_sSVD(n,2) = mean(MTT_sSVD(:));
std_sSVD(n,2) = std(MTT_sSVD(:));
mean_sSVD(n,3) = mean(CBV_sSVD(:));
std_sSVD(n,3) = std(CBV_sSVD(:));
mean_sSVD(n,4) = mean(TTP_sSVD(:));
std_sSVD(n,4) = std(TTP_sSVD(:));

% Method 2: block-circulant SVD (bSVD)
RIF_bSVD = pct_bsvd(Cn,AIF,1,lambda,m,ones(X,Y));
CBF_bSVD = squeeze(max(RIF_bSVD));
MTT_bSVD = squeeze(sum(RIF_bSVD))./CBF_bSVD;
CBV_bSVD = squeeze(sum(RIF_bSVD))/60;
TTP_bSVD = pct_ttp(pct_tec(AIF,RIF_bSVD),1);
mean_bSVD(n,1) = mean(CBF_bSVD(:));
std_bSVD(n,1) = std(CBF_bSVD(:));
mean_bSVD(n,2) = mean(MTT_bSVD(:));
std_bSVD(n,2) = std(MTT_bSVD(:));
mean_bSVD(n,3) = mean(CBV_bSVD(:));
std_bSVD(n,3) = std(CBV_bSVD(:));
mean_bSVD(n,4) = mean(TTP_bSVD(:));
std_bSVD(n,4) = std(TTP_bSVD(:));

% % Method 2-1: block-circulant SVD (bSVD)no noise
% RIF0_bSVD = pct_bsvd(Cv,AIF,1,lambda,m,ones(X,Y));
% CBF0_bSVD = squeeze(max(RIF0_bSVD));
% MTT0_bSVD = squeeze(sum(RIF0_bSVD))./CBF_bSVD;
% CBV0_bSVD = squeeze(sum(RIF0_bSVD))/60;
% mean0_bSVD(n,1) = mean(CBF0_bSVD(:));
% std0_bSVD(n,1) = std(CBF0_bSVD(:));
% mean0_bSVD(n,2) = mean(MTT0_bSVD(:));
% std0_bSVD(n,2) = std(MTT0_bSVD(:));
% mean0_bSVD(n,3) = mean(CBV0_bSVD(:));
% std0_bSVD(n) = std(CBV0_bSVD(:));

% Method 3: Tikhonov
RIF_tikh = pct_tikh(Cn,AIF,1,lambda,1,ones(X,Y));
CBF_tikh = squeeze(max(RIF_tikh));
MTT_tikh = squeeze(sum(RIF_tikh))./CBF_tikh;
CBV_tikh = squeeze(sum(RIF_tikh))/60;
TTP_tikh = pct_ttp(pct_tec(AIF,RIF_tikh),1);
mean_tikh(n,1) = mean(CBF_tikh(:));
std_tikh(n,1) = std(CBF_tikh(:));
mean_tikh(n,2) = mean(MTT_tikh(:));
std_tikh(n,2) = std(MTT_tikh(:));
mean_tikh(n,3) = mean(CBV_tikh(:));
std_tikh(n,3) = std(CBV_tikh(:));
mean_tikh(n,4) = mean(TTP_tikh(:));
std_tikh(n,4) = std(TTP_tikh(:));

% Method 4: Sparse Perfusion Deconvolution (SPD)
% set parameters %
% params.xdic = CBF_gt;
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
mean_SPD(n,1) = mean(CBF_SPD(:));
std_SPD(n,1) = std(CBF_SPD(:));
mean_SPD(n,2) = mean(MTT_SPD(:));
std_SPD(n,2) = std(MTT_SPD(:));
mean_SPD(n,3) = mean(CBV_SPD(:));
std_SPD(n,3) = std(CBV_SPD(:));
mean_SPD(n,4) = mean(TTP_SPD(:));
std_SPD(n,4) = std(TTP_SPD(:));

% Method 5: Time-Intensity Profile Similarity (TIPS) bilateral filtering+bSVD 
w     = 5;       % bilateral filter half-width
sigma_tips = [5 sigma]; % bilateral filter standard deviations

Cn_TIPS = permute(pct_tips(permute(Cn,[2 3 1]),w,sigma_tips),[3 1 2]);
RIF_TIPS_bSVD = pct_bsvd(Cn_TIPS,AIF,1,lambda,m,ones(X,Y));
CBF_TIPS_bSVD = squeeze(max(RIF_TIPS_bSVD));
MTT_TIPS_bSVD = squeeze(sum(RIF_TIPS_bSVD))./CBF_TIPS_bSVD;
CBV_TIPS_bSVD = squeeze(sum(RIF_TIPS_bSVD))/60;
TTP_TIPS_bSVD = pct_ttp(pct_tec(AIF,RIF_TIPS_bSVD),1);
mean_TIPS_bSVD(n,1) = mean(CBF_TIPS_bSVD(:));
std_TIPS_bSVD(n,1) = std(CBF_TIPS_bSVD(:));
mean_TIPS_bSVD(n,2) = mean(MTT_TIPS_bSVD(:));
std_TIPS_bSVD(n,2) = std(MTT_TIPS_bSVD(:));
mean_TIPS_bSVD(n,3) = mean(CBV_TIPS_bSVD(:));
std_TIPS_bSVD(n,3) = std(CBV_TIPS_bSVD(:));
mean_TIPS_bSVD(n,4) = mean(TTP_TIPS_bSVD(:));
std_TIPS_bSVD(n,4) = std(TTP_TIPS_bSVD(:));

% Method 6: Tensor Total Variation (TTV)
% Parameters for TV-based Methods
input.reg = reg;
input.n1=X;input.n2=Y; input.nt = T;
input.maxitr=500;
input.L=1;input.num=1; 
input.l=0; input.u=100;
input.no = 50;
Ca = pct_circ(AIF,[],[],1); % block-circulant version
input.A = Ca;
input.b = reshape(Cn,T,[]);

% TTV
out_ttv = TTV_FCSA(input);
RIF_TTV = reshape(out_ttv.y,[T,X,Y]);
CBF_TTV = squeeze(max(RIF_TTV));
MTT_TTV = squeeze(sum(RIF_TTV))./CBF_TTV;
CBV_TTV = squeeze(sum(RIF_TTV))/60;
TTP_TTV = pct_ttp(pct_tec(AIF,RIF_TTV),1);
mean_TTV(n,1) = mean(CBF_TTV(:));
std_TTV(n,1) = std(CBF_TTV(:));
mean_TTV(n,2) = mean(MTT_TTV(:));
std_TTV(n,2) = std(MTT_TTV(:));
mean_TTV(n,3) = mean(CBV_TTV(:));
std_TTV(n,3) = std(CBV_TTV(:));
mean_TTV(n,4) = mean(TTP_TTV(:));
std_TTV(n,4) = std(TTP_TTV(:));

% Method 7: TIPS+TTV
input.b = reshape(Cn_TIPS,T,[]);
out_tips_ttv = TTV_FCSA(input);
RIF_TIPS_TTV = reshape(out_ttv.y,[T,X,Y]);
CBF_TIPS_TTV = squeeze(max(RIF_TIPS_TTV));
MTT_TIPS_TTV = squeeze(sum(RIF_TIPS_TTV))./CBF_TIPS_TTV;
CBV_TIPS_TTV = squeeze(sum(RIF_TIPS_TTV))/60;
TTP_TIPS_TTV = pct_ttp(pct_tec(AIF,RIF_TIPS_TTV),1);
mean_TIPS_TTV(n,1) = mean(CBF_TIPS_TTV(:));
std_TIPS_TTV(n,1) = std(CBF_TIPS_TTV(:));
mean_TIPS_TTV(n,2) = mean(MTT_TIPS_TTV(:));
std_TIPS_TTV(n,2) = std(MTT_TIPS_TTV(:));
mean_TIPS_TTV(n,3) = mean(CBV_TIPS_TTV(:));
std_TIPS_TTV(n,3) = std(CBV_TIPS_TTV(:));
mean_TIPS_TTV(n,4) = mean(TTP_TIPS_TTV(:));
std_TIPS_TTV(n,4) = std(TTP_TIPS_TTV(:));

end

%% Plot results
close all;
fs = 30; % font size for x and y axises
lw = 4;
ms = 10;
spec = {'+-b','d--k','vc','o-.g','.-.m','s-r','^y'};

% CBF
h1 = figure;
line(0:100,0:100,'Color',[0.5 0.5 0.5],'LineWidth',5);
set(gca,'XTick',0:20:100); xlim([10 85]);
hold on;
errorbar(CBFs,mean_sSVD(:,1), std_sSVD(:,1),spec{1},'MarkerSize',ms,'MarkerFaceColor',spec{1}(end),'LineWidth',lw);
errorbar(CBFs,mean_bSVD(:,1), std_bSVD(:,1),spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
errorbar(CBFs,mean_TIPS_bSVD(:,1), std_TIPS_bSVD(:,1),spec{3},'MarkerFaceColor',spec{3}(end),'MarkerSize',ms,'LineWidth',lw);
errorbar(CBFs,mean_tikh(:,1), std_bSVD(:,1),spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
errorbar(CBFs,mean_SPD(:,1), std_SPD(:,1),spec{5},'MarkerSize',ms,'MarkerFaceColor',spec{5}(end),'LineWidth',lw);
errorbar(CBFs,mean_TTV(:,1), std_TTV(:,1),spec{6},'MarkerSize',ms,'MarkerFaceColor',spec{6}(end),'LineWidth',lw);
errorbar(CBFs,mean_TIPS_TTV(:,1), std_TIPS_TTV(:,1),spec{7},'MarkerSize',ms,'MarkerFaceColor',spec{7}(end),'LineWidth',lw);
set(gca,'FontSize',fs);
xlabel('True CBF (ml/100g/min)');
ylabel('Estimated CBF (ml/100g/min)');
legend('Reference','sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');


% Variation (in STD)
h2 = figure;
plot(CBFs,std_sSVD(:,1),spec{1},'MarkerFaceColor',spec{1}(end),'MarkerSize',ms,'LineWidth',lw);
hold on; xlim([8 82]);
plot(CBFs,std_bSVD(:,1),spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
plot(CBFs,std_TIPS_bSVD(:,1),spec{3},'MarkerFaceColor',spec{3}(end),'MarkerSize',ms,'LineWidth',lw);
plot(CBFs,std_tikh(:,1),spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
plot(CBFs,std_SPD(:,1),spec{5},'MarkerFaceColor',spec{5}(end),'MarkerSize',ms,'LineWidth',lw);
plot(CBFs,std_TTV(:,1),spec{6},'MarkerFaceColor',spec{6}(end),'MarkerSize',ms,'LineWidth',lw);
plot(CBFs,std_TIPS_TTV(:,1),spec{7},'MarkerFaceColor',spec{7}(end),'MarkerSize',ms,'LineWidth',lw);
set(gca,'FontSize',fs);
xlabel('True CBF (ml/100g/min)');
ylabel('CBF Variation (ml/100g/min)');
legend('sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');


% MTT
h3 = figure;
line(0:24,0:24,'Color',[0.5 0.5 0.5],'LineWidth',5);
xlim([2 25]);
hold on;
errorbar(MTTs,mean_sSVD(:,2), std_sSVD(:,2),spec{1},'MarkerSize',ms,'MarkerFaceColor',spec{1}(end),'LineWidth',lw);
errorbar(MTTs,mean_bSVD(:,2), std_bSVD(:,2),spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
errorbar(MTTs,mean_TIPS_bSVD(:,2), std_TIPS_bSVD(:,2),spec{3},'MarkerFaceColor',spec{3}(end),'MarkerSize',ms,'LineWidth',lw);
errorbar(MTTs,mean_tikh(:,2), std_bSVD(:,2),spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
errorbar(MTTs,mean_SPD(:,2), std_SPD(:,2),spec{5},'MarkerSize',ms,'MarkerFaceColor',spec{5}(end),'LineWidth',lw);
errorbar(MTTs,mean_TTV(:,2), std_TTV(:,2),spec{6},'MarkerSize',ms,'MarkerFaceColor',spec{6}(end),'LineWidth',lw);
errorbar(MTTs,mean_TIPS_TTV(:,2), std_TIPS_TTV(:,2),spec{7},'MarkerSize',ms,'MarkerFaceColor',spec{7}(end),'LineWidth',lw);
set(gca,'FontSize',fs);
xlabel('True MTT (s)');
ylabel('Estimated MTT (s)');
legend('Reference','sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');

% Variation (in STD)
h4  = figure; xlim([2 25]); hold on;
plot(MTTs,std_sSVD(:,2),spec{1},'MarkerFaceColor',spec{1}(end),'MarkerSize',ms,'LineWidth',lw);
plot(MTTs,std_bSVD(:,2),spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
plot(MTTs,std_TIPS_bSVD(:,2),spec{3},'MarkerFaceColor',spec{3}(end),'MarkerSize',ms,'LineWidth',lw);
plot(MTTs,std_tikh(:,2),spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
plot(MTTs,std_SPD(:,2),spec{5},'MarkerFaceColor',spec{5}(end),'MarkerSize',ms,'LineWidth',lw);
plot(MTTs,std_TTV(:,2),spec{6},'MarkerFaceColor',spec{6}(end),'MarkerSize',ms,'LineWidth',lw);
plot(MTTs,std_TIPS_TTV(:,2),spec{7},'MarkerFaceColor',spec{7}(end),'MarkerSize',ms,'LineWidth',lw);
set(gca,'FontSize',fs);
xlabel('True MTT (s)');
ylabel('MTT Variation (s)');
legend('sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');


%% RMSE and Lin's Concordance
% CBF
rmse_sSVD(1) = pct_rmse(mean_sSVD(:,1),CBFs);
rmse_bSVD(1) = pct_rmse(mean_bSVD(:,1),CBFs);
rmse_TIPS_bSVD(1) = pct_rmse(mean_TIPS_bSVD(:,1),CBFs);
rmse_tikh(1) = pct_rmse(mean_tikh(:,1),CBFs);
rmse_SPD(1) = pct_rmse(mean_SPD(:,1),CBFs);
rmse_TTV(1) = pct_rmse(mean_TTV(:,1),CBFs);
rmse_TIPS_TTV(1) = pct_rmse(mean_TIPS_TTV(:,1),CBFs);

[lin_sSVD(1),ci_sSVD] = pct_lincon(mean_sSVD(:,1),CBFs);
[lin_bSVD(1),ci_bSVD] = pct_lincon(mean_bSVD(:,1),CBFs);
[lin_TIPS_bSVD(1),ci_TIPS_bSVD] = pct_lincon(mean_TIPS_bSVD(:,1),CBFs);
[lin_tikh(1),ci_tikh] = pct_lincon(mean_tikh(:,1),CBFs);
[lin_SPD(1),ci_SPD] = pct_lincon(mean_SPD(:,1),CBFs);
[lin_TTV(1),ci_TTV] = pct_lincon(mean_TTV(:,1),CBFs);
[lin_TIPS_TTV(1),ci_TIPS_TTV] = pct_lincon(mean_TIPS_TTV(:,1),CBFs);

% MTT
rmse_sSVD(2) = pct_rmse(mean_sSVD(:,2),MTTs);
rmse_bSVD(2) = pct_rmse(mean_bSVD(:,2),MTTs);
rmse_TIPS_bSVD(2) = pct_rmse(mean_TIPS_bSVD(:,2),MTTs);
rmse_tikh(2) = pct_rmse(mean_tikh(:,2),MTTs);
rmse_SPD(2) = pct_rmse(mean_SPD(:,2),MTTs);
rmse_TTV(2) = pct_rmse(mean_TTV(:,2),MTTs);
rmse_TIPS_TTV(2) = pct_rmse(mean_TIPS_TTV(:,2),MTTs);

[lin_sSVD(2),ci_sSVD] = pct_lincon(mean_sSVD(:,2),MTTs);
[lin_bSVD(2),ci_bSVD] = pct_lincon(mean_bSVD(:,2),MTTs);
[lin_TIPS_bSVD(2),ci_TIPS_bSVD] = pct_lincon(mean_TIPS_bSVD(:,2),MTTs);
[lin_tikh(2),ci_tikh] = pct_lincon(mean_tikh(:,2),MTTs);
[lin_SPD(2),ci_SPD] = pct_lincon(mean_SPD(:,2),MTTs);
[lin_TTV(2),ci_TTV] = pct_lincon(mean_TTV(:,2),MTTs);
[lin_TIPS_TTV(2),ci_TIPS_TTV] = pct_lincon(mean_TIPS_TTV(:,2),MTTs);

% CBV
rmse_sSVD(3) = pct_rmse(mean_sSVD(:,3),CBVs);
rmse_bSVD(3) = pct_rmse(mean_bSVD(:,3),CBVs);
rmse_TIPS_bSVD(3) = pct_rmse(mean_TIPS_bSVD(:,3),CBVs);
rmse_tikh(3) = pct_rmse(mean_tikh(:,3),CBVs);
rmse_SPD(3) = pct_rmse(mean_SPD(:,3),CBVs);
rmse_TTV(3) = pct_rmse(mean_TTV(:,3),CBVs);
rmse_TIPS_TTV(3) = pct_rmse(mean_TIPS_TTV(:,3),CBVs);

% TTP
rmse_sSVD(4) = pct_rmse(mean_sSVD(:,4),TTPs);
rmse_bSVD(4) = pct_rmse(mean_bSVD(:,4),TTPs);
rmse_TIPS_bSVD(4) = pct_rmse(mean_TIPS_bSVD(:,4),TTPs);
rmse_tikh(4) = pct_rmse(mean_tikh(:,4),TTPs);
rmse_SPD(4) = pct_rmse(mean_SPD(:,4),TTPs);
rmse_TTV(4) = pct_rmse(mean_TTV(:,4),TTPs);
rmse_TIPS_TTV(4) = pct_rmse(mean_TIPS_TTV(:,4),TTPs);

[lin_sSVD(4),ci_sSVD] = pct_lincon(mean_sSVD(:,4),TTPs);
[lin_bSVD(4),ci_bSVD] = pct_lincon(mean_bSVD(:,4),TTPs);
[lin_TIPS_bSVD(4),ci_TIPS_bSVD] = pct_lincon(mean_TIPS_bSVD(:,4),TTPs);
[lin_tikh(4),ci_tikh] = pct_lincon(mean_tikh(:,4),TTPs);
[lin_SPD(4),ci_SPD] = pct_lincon(mean_SPD(:,4),TTPs);
[lin_TTV(4),ci_TTV] = pct_lincon(mean_TTV(:,4),TTPs);
[lin_TIPS_TTV(4),ci_TIPS_TTV] = pct_lincon(mean_TIPS_TTV(:,4),TTPs);

results = zeros(7,8);
results(:,1:2:end) = [rmse_sSVD; rmse_bSVD; rmse_TIPS_bSVD; rmse_tikh; rmse_SPD; rmse_TTV; rmse_TIPS_TTV];
results(:,2:2:end) = [lin_sSVD; lin_bSVD; lin_TIPS_bSVD; lin_tikh; lin_SPD; lin_TTV; lin_TIPS_TTV];
results(:,6) = [];

%% Save figures
figpath = '../../../Figures/pdf';
w = 12; h = 10;

set(h1, 'papersize', [w h]);
set(h1, 'paperposition', [0 0 w h]);
print(h1,fullfile(figpath,'cbfs.pdf'),'-dpdf');

set(h2, 'papersize', [w h]);
set(h2, 'paperposition', [0 0 w h]);
print(h2,fullfile(figpath,'cbfs_var.pdf'),'-dpdf');

set(h3, 'papersize', [w h]);
set(h3, 'paperposition', [0 0 w h]);
print(h3,fullfile(figpath,'mtts.pdf'),'-dpdf');

set(h4, 'papersize', [w h]);
set(h4, 'paperposition', [0 0 w h]);
print(h4,fullfile(figpath,'mtts_var.pdf'),'-dpdf');

