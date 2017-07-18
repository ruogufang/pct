% Experiment: Noise Power Spectrum & RIF recovery with tracer delay in AIF
% Figure 2
% Ruogu Fang 06/16/2014

close all; clear; clc;

% Set path
addpath(genpath('toolbox')); % Include toolbox package
addpath(genpath('~/Dropbox/Research/pct')); % put pct package into toolbox finally

%% Set parameters for noise power spectrum
mA = 10; % low-dose tube current
mA0 = 190; % high-dose tube current
K = 103; % constant relating tube current to noise standard deviation
N = 256; % size of the image
sigma = pct_mA2sigma(mA,mA0,K); % noise standard variation
CF = 3.06; % factor to account for the smoothing effect in correlation

% Compute the power spectrum of the phantom
% [fname, pname] = uigetfile('../DICOM_noise_vs_mA/*.dcm','Select files:','MultiSelect','on');
fname = '400.dcm';
vol = dicomread(fname);
vol = vol(145:400,145:400);
vol = double(vol);
noise_phantom = vol-mean(vol(:));
ps_phantom = pct_powspec(noise_phantom);
ps_phantom = pct_rotavg(ps_phantom);

% Generated noise by computing the auto-correlation function directly from
% phantom data
noise_gen = randn(size(noise_phantom));
auto_corr = conv2(noise_phantom(end:-1:1,end:-1:1),noise_phantom,'same');
sac = auto_corr(N/2-5:N/2+5,N/2-5:N/2+5);
noise_gen = conv2(noise_gen,sac,'same');
noise_gen = (noise_gen-mean(noise_gen(:)))./std(noise_gen(:)) * sigma*CF;
ps = pct_powspec(noise_gen);
ps = pct_rotavg(ps);

%% Parameters
mA = 10;
reg = 100; % TV regularization parameter
lambda = 0.15; % bSVD truncation threshold
m = 2; % circulant number
CBV = 4; %ml/100g
CBF = 20; % 20 : 10 : 80 ml/100g/min
X = 1; Y = 1; Z = 1; T = 60; % volumn size
Tn = 30;

% Step 1: Generate Artery Input Function (AIF)
t = (0 : T-1)';
t0 = 0;
a = 1;
b = 3;
c = 1.5;
AIF0 = pct_aif(t,0,a,b,c);
AIF = pct_aif(t,t0,a,b,c);

% Step 2: Generate Impulse Residue Function RIF using Gamma Family
MTT = CBV / CBF * 60; %sec (3 ~ 12 s)
alpha = 10;
beta = MTT / alpha;
RIF = pct_irf(t, alpha, beta); % Gamma RIF
% RIF = exp(-(t-t0)/MTT); % exponential RIF

% Step 3: Generate Tissue Time Enhancement Curve
TEC = pct_tec(AIF0, RIF, CBF);

% Step 4: Correlated Gaussian Noise Generation
load acf;
% sigma = max(TEC(:)) / (10^(PSNR/20));
sigma = pct_mA2sigma(mA,190);
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

% Method 4: TIPS+bSVD
% Set bilateral filter parameters.
w     = 5;       % bilateral filter half-width
sigma_tips = [5 sigma]; % bilateral filter standard deviations

Cn_TIPS = pct_tips(Cn,w,sigma_tips);
RIF_TIPS_bSVD = pct_bsvd(Cn_TIPS,AIF,1,lambda,m,ones(X,Y));

% Method 5: Total Variation on RIF

% Parameters
m = 1;
input.reg = reg;
input.maxitr=500;
input.L=1;input.num=1; 
input.l=0; input.u=80;
input.no = 50;
[Ca,Cc] = pct_circ(AIF,Cn,[],m); % block-circulant version
Tc = T*m;
input.n1=X;input.n2=Y; input.nt = Tc;
input.A = Ca;
input.b = reshape(Cc,Tc,[]);
input.dim = [Tc X Y];

% Run TTV
out = TTV_FCSA(input);
RIF_TTV = reshape(out.y,[Tc,X,Y]);
RIF_TTV = RIF_TTV(1:T,:,:); % keep first T data points

% Method 6: TIPS+TTV
input.b = reshape(Cn_TIPS,Tc,[]);
out_tipsttv = TTV_FCSA(input);
RIF_TIPS_TTV = reshape(out_tipsttv.y,[Tc,X,Y]);
RIF_TIPS_TTV = RIF_TIPS_TTV(1:T,:,:); % keep first T data points

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot power spectrum
h = figure;
subplot(241);
hold on;
plot(1:length(ps_phantom),ps_phantom,'k','LineWidth',6);
plot(1:length(ps),ps,'b--','LineWidth',6);
set(gca,'FontSize',20);
xlabel('(a) Frequency (Hz)');
ylabel('Noise Power');
legend('Reference','Simulated');

fs = 20; % font size
ymax = CBF+10; % max value for y axis
ymin = -10; % min value for y axis

% subplot(242);
% hold on;
% set(gca,'FontSize',fs);
% plot(AIF,'r','LineWidth',6);
% xlabel('(a) t(s)','FontSize',fs);
% ylabel('AIF(t)','FontSize',fs);
% grid on;
% xlim([0 T]);
% ylim([-1 max(AIF)+1]);
% title('Arterial Input Function');

% Plot RIF for sSVD
% sSVD
subplot(242);
set(gca,'FontSize',fs);
plot(t, mean(reshape(RIF_sSVD,T,[]),2),'k','LineWidth',6);
xlabel('(b) t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 T]);
ylim([ymin ymax]);
title('sSVD');

% bSVD
subplot(243);
set(gca,'FontSize',fs);
plot(t, mean(reshape(RIF_bSVD,T,[]),2),'k','LineWidth',6);
xlabel('(c) t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 T]);
ylim([ymin ymax]);
title('bSVD');

% TIPS + bSVD
subplot(244);
set(gca,'FontSize',fs);
plot(t, mean(reshape(RIF_TIPS_bSVD,T,[]),2),'k','LineWidth',6);
xlabel('(d) t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 T]);
ylim([ymin ymax]);
title('TIPS+bSVD');

% Tikhonov
subplot(245);
set(gca,'FontSize',fs);
plot(t, mean(reshape(RIF_tikh,T,[]),2),'k','LineWidth',6);
xlabel('(e) t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 T]);
ylim([ymin ymax]);
title('Tikhonov');


% TTV
subplot(246);
set(gca,'FontSize',fs);
plot(t, mean(reshape(RIF_TTV,T,[]),2),'k','LineWidth',6);
xlabel('(f) t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 T]);
ylim([ymin ymax]);
title('TTV');

% TIPS + TTV
subplot(247);
set(gca,'FontSize',fs);
plot(t, mean(reshape(RIF_TIPS_TTV,T,[]),2),'k','LineWidth',6);
xlabel('(g) t(s)');
ylabel('CBF.R(t)');
grid on;
xlim([0 T]);
ylim([ymin ymax]);
title('TIPS+TTV');

% ideal residue function RIF
subplot(248);
set(gca,'FontSize',fs);
plot(RIF * CBF,'k','LineWidth',6);
xlabel('(h) t(s)','FontSize',fs);
ylabel('BF.R(t)','FontSize',fs);
grid on;
xlim([0 T]);
ylim([ymin ymax]);
title('Reference');

%% save figure
w = 22; h=10;
set(gcf, 'papersize', [w h]);
set(gcf, 'paperposition', [0 0 w h]);
print('../../../Figures/pdf/residue.pdf','-dpdf');
saveas(gcf,'../../../Figures/fig/residue.fig');
