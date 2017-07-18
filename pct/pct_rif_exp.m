% Simulation Experiment on Residue Function Recovery
%
%   Ruogu Fang Revised 09/24/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

clear; clc; close all;

% Parameters
mA = 15;
reg = 10; % TV regularization parameter
lambda = 0.06; % bSVD truncation threshold
m = 2; % circulant number
CBV = 4; %ml/100g
CBF = 20; % 20 : 10 : 80 ml/100g/min
X = 1; Y = 1; Z = 1; T = 60; % volumn size
Tn = 30;

% Step 1: Generate Artery Input Function (AIF)
t = (0 : T-1)';
t0 = 5;
a = 1;
b = 3;
c = 1.5;
AIF = pct_aif(t,t0,a,b,c);

% Step 2: Generate Impulse Residue Function RIF using Gamma Family
MTT = CBV / CBF * 60; %sec (3 ~ 12 s)
alpha = 10;
beta = MTT / alpha;
RIF = pct_irf(t, alpha, beta); % Gamma RIF
% RIF = exp(-(t-t0)/MTT); % exponential RIF

% Step 3: Generate Tissue Time Enhancement Curve
TEC = pct_tec(AIF, RIF, CBF);

% Step 4: Correlated Gaussian Noise Generation
load acf;
% sigma = max(TEC(:)) / (10^(PSNR/20));
% sigma = pct_mA2sigma(mA,190);
sigma=0;
PSNR = 20*log10(max(TEC(:))/sigma);
Cv = repmat(TEC,[1 X Y Z]);
Cn = pct_noise(Cv, acf, sigma,'s');

% Method 1: Standard Truncated Singular Value Decomposition (sSVD)
% High-dose
RIF0_sSVD = pct_ssvd(Cv,AIF,1,lambda,ones(X,Y));
CBF0_sSVD = squeeze(max(RIF0_sSVD));
% Low-dose
RIF_sSVD = pct_ssvd(Cn,AIF,1,lambda,ones(X,Y));
CBF_sSVD = squeeze(max(RIF_sSVD));


% Method 2: Block-circulant SVD (bSVD)
% High-dose
RIF0_bSVD = pct_bsvd(Cv,AIF,1,lambda,m,ones(X,Y));
CBF0_bSVD = squeeze(max(RIF0_bSVD));
% Low-dose
RIF_bSVD = pct_bsvd(Cn,AIF,1,lambda,m,ones(X,Y));
CBF_bSVD = squeeze(max(RIF_bSVD));

% Method 3: Tikhonov Regularization
RIF_tikh = pct_tikh(Cn,AIF,1,lambda,1,ones(X,Y));
CBF_tikh = squeeze(max(RIF_tikh));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot RIF for sSVD
% ideal residue function RIF
h = figure;
subplot(231);
set(gca,'FontSize',20);
plot(RIF * CBF,'k','LineWidth',6);
xlabel('t(s)','FontSize',20);
ylabel('BF.R(t)','FontSize',20);
grid on;
xlim([0 60]);
ylim([-10 CBF+10]);
title('Reference');

% sSVD
subplot(232);
set(gca,'FontSize',20);
plot(t, mean(reshape(RIF_sSVD,T,[]),2),'k','LineWidth',4);
xlabel('t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 60]);
ylim([-10 CBF+10]);
title('sSVD');

% bSVD
subplot(233);
set(gca,'FontSize',20);
plot(t, mean(reshape(RIF_bSVD,T,[]),2),'k','LineWidth',4);
xlabel('t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 60]);
ylim([-10 CBF+10]);
title('bSVD');

% Tikhonov
subplot(234);
set(gca,'FontSize',20);
plot(t, mean(reshape(RIF_tikh,T,[]),2),'k','LineWidth',4);
xlabel('t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 60]);
ylim([-10 CBF+10]);
title('Tikhonov');


% sSVD no noise
subplot(236);
set(gca,'FontSize',20);
plot(t, mean(reshape(RIF0_sSVD,T,[]),2),'k','LineWidth',6);
xlabel('t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 60]);
ylim([-10 CBF+10]);
title('sSVD no noise');

