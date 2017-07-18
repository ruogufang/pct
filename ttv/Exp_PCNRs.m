% Simulation Experiment on CBF estimation at various PCNR values
% Figure 6
%   Ruogu Fang Revised 09/24/2013
%
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

clear; clc; close all;

% Set path
addpath(genpath('~/Dropbox/Research/pct')); % Include PCT package
addpath(genpath('toolbox')); % SPD toolbox


% Parameters
% PCNRs = zeros(6,4);
% PCNRs(:,1) = [2 5 7 10 20 30]'; % 1 or 0.2
PCNRs = 0.2;
lambda = 0.15; % bSVD truncation threshold
m = 2; % circulant number
CBV1 = 4; %ml/100g
CBV2 = 2; %ml/100g
CBF1 = 70; % ml/100g/min
CBF2 = 20; % ml/100g/min
X = 40; Y = 40; Z = 1; T = 60; % volumn size
alpha = 10; % transpose function h parameter
load acf;


% Step 1: Generate Artery Input Function (AIF)
t = (0 : T-1)';
t0 = 0;
a = 1;
b = 3;
c = 1.5;
AIF = pct_aif(t,t0,a,b,c);

% Step 2: Generate Impulse Residue Function RIF using Gamma Family
CBF_gt = zeros(X,Y);
CBF_gt(:,1:Y/2) = CBF1;
CBF_gt(:,Y/2+1:end) = CBF2;
MTT1 = CBV1 / CBF1 * 60; %sec (3 ~ 12 s)
beta1 = MTT1 / alpha;
RIF1 = pct_irf(t, alpha, beta1);
MTT2 = CBV2 / CBF2 * 60; %sec (3 ~ 12 s)
beta2 = MTT2 / alpha;
RIF2 = pct_irf(t, alpha, beta2);
RIF = zeros(T,X,Y);
RIF(:,:,1:Y/2) = repmat(RIF1,[1 X Y/2]);
RIF(:,:,Y/2+1:end) = repmat(RIF2,[1 X Y/2]);
MTT_gt = zeros(X,Y);
MTT_gt(:,1:Y/2) = MTT1;
MTT_gt(:,Y/2+1:end) = MTT2;
CBV_gt = zeros(X,Y);
CBV_gt(:,1:Y/2) = CBV1;
CBV_gt(:,Y/2+1:end) = CBV2;

% Step 3: Generate Tissue Time Enhancement Curve
TEC = pct_tec(AIF, RIF, CBF_gt);
TEC1 = pct_tec(AIF, RIF1,CBF1);
TEC2 = pct_tec(AIF, RIF2, CBF2);
TTP_gt = pct_ttp(TEC,1);
TTP1 = mean(mean(TTP_gt(:,1:Y/2)));
TTP2 = mean(mean(TTP_gt(:,Y/2+1:end)));

PCNRs(:,2) = PCNRs(:,1)/max(abs(CBF1-CBF2))*max(abs(CBV1-CBV2));
PCNRs(:,3) = PCNRs(:,1)/max(abs(CBF1-CBF2))*max(abs(MTT1-MTT2));
PCNRs(:,4) = PCNRs(:,1)/max(abs(CBF1-CBF2))*max(abs(TTP1-TTP2));

% Initialization
PCNR_sSVD = zeros(length(PCNRs),4);
PCNR_bSVD = zeros(length(PCNRs),4);
PCNR_TIPS_bSVD = zeros(length(PCNRs),4);
PCNR_tikh = zeros(length(PCNRs),4);
PCNR_SPD = zeros(length(PCNRs),4);
PCNR_TTV = zeros(length(PCNRs),4);
PCNR_TIPS_TTV = zeros(length(PCNRs),4);

for i = 1 : size(PCNRs,1)
    PCNR = PCNRs(i,1);
    % Step 4: Correlated Gaussian Noise Generation
%     sigma = (max(TEC1(:)-TEC2(:))) / (10^(PCNR/20));
sigma = max(abs(CBF1-CBF2))/PCNR;
Cn = pct_noise(TEC, acf, sigma,'s');
    
% Method 1: Standard Truncated Singular Value Decomposition (sSVD)
RIF_sSVD = pct_ssvd(Cn,AIF,1,lambda,ones(X,Y));
CBF_sSVD = squeeze(max(RIF_sSVD));
MTT_sSVD = squeeze(sum(RIF_sSVD))./CBF_sSVD;
CBV_sSVD = squeeze(sum(RIF_sSVD))/60;
TTP_sSVD = pct_ttp(pct_tec(AIF,RIF_sSVD),1);
CBF1_sSVD = CBF_sSVD(:,1:Y/2);
CBF2_sSVD = CBF_sSVD(:,Y/2+1:end);
PCNR_sSVD(i,1) = max(abs(mean(mean(CBF_sSVD(:,1:Y/2)))-mean(mean(CBF_sSVD(:,Y/2+1:end)))))/sigma;
PCNR_sSVD(i,2) = max(abs(mean(mean(CBV_sSVD(:,1:Y/2)))-mean(mean(CBV_sSVD(:,Y/2+1:end)))))/sigma;
PCNR_sSVD(i,3) = max(abs(mean(mean(MTT_sSVD(:,1:Y/2)))-mean(mean(MTT_sSVD(:,Y/2+1:end)))))/sigma;
PCNR_sSVD(i,4) = max(abs(mean(mean(TTP_sSVD(:,1:Y/2)))-mean(mean(TTP_sSVD(:,Y/2+1:end)))))/sigma;

% Method 2: block-circulant SVD (bSVD)
RIF_bSVD = pct_bsvd(Cn,AIF,1,lambda,m,ones(X,Y));
CBF_bSVD = squeeze(max(RIF_bSVD));
MTT_bSVD = squeeze(sum(RIF_bSVD))./CBF_bSVD;
CBV_bSVD = squeeze(sum(RIF_bSVD))/60;
TTP_bSVD = pct_ttp(pct_tec(AIF,RIF_bSVD),1);
PCNR_bSVD(i,1) = max(abs(mean(mean(CBF_bSVD(:,1:Y/2)))-mean(mean(CBF_bSVD(:,Y/2+1:end)))))/sigma;
PCNR_bSVD(i,2) = max(abs(mean(mean(CBV_bSVD(:,1:Y/2)))-mean(mean(CBV_bSVD(:,Y/2+1:end)))))/sigma;
PCNR_bSVD(i,3) = max(abs(mean(mean(MTT_bSVD(:,1:Y/2)))-mean(mean(MTT_bSVD(:,Y/2+1:end)))))/sigma;
PCNR_bSVD(i,4) = max(abs(mean(mean(TTP_bSVD(:,1:Y/2)))-mean(mean(TTP_bSVD(:,Y/2+1:end)))))/sigma;


% Method 3: Tikhonov
RIF_tikh = pct_tikh(Cn,AIF,1,lambda,1,ones(X,Y));
CBF_tikh = squeeze(max(RIF_tikh));
MTT_tikh = squeeze(sum(RIF_tikh))./CBF_tikh;
CBV_tikh = squeeze(sum(RIF_tikh))/60;
TTP_tikh = pct_ttp(pct_tec(AIF,RIF_tikh),1);
PCNR_tikh(i,1) = max(abs(mean(mean(CBF_tikh(:,1:Y/2)))-mean(mean(CBF_tikh(:,Y/2+1:end)))))/sigma;
PCNR_tikh(i,2) = max(abs(mean(mean(CBV_tikh(:,1:Y/2)))-mean(mean(CBV_tikh(:,Y/2+1:end)))))/sigma;
PCNR_tikh(i,3) = max(abs(mean(mean(MTT_tikh(:,1:Y/2)))-mean(mean(MTT_tikh(:,Y/2+1:end)))))/sigma;
PCNR_tikh(i,4) = max(abs(mean(mean(TTP_tikh(:,1:Y/2)))-mean(mean(TTP_tikh(:,Y/2+1:end)))))/sigma;


% Parameters for TV-based Methods
Ca = pct_circ(AIF,[],[],1); % block-circulant version
input.A = Ca;
input.b = reshape(Cn,T,[]);

% Method 4: Sparse Perfusion Deconvolution (SPD)
% set parameters %
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
% params.lambda = 10; % weight of the sparsity term in SPAMS package
params.mixture = 1; % weight of input signal
params.J = []; % initialize energy function
params.round = 2; % MAP rounds
params.C = Cn;
params.AIF = AIF;
params.R = RIF_bSVD;
params.verbose = false;
params.rho = 1.05;
params.m = m;
params.beta = 1; % weight of prior
params.blocksize = 6;
overlap = 5; % patch overlap
params.stepsize = params.blocksize - overlap;

% %Load pre-trained dictionary
load D_ODL_6x6.mat
params.D = dict;

% Online deconvolution
[CBF_SPD, MTT_SPD, CBV_SPD, TTP_SPD] = spd(params);
PCNR_SPD(i,1) = max(abs(mean(mean(CBF_SPD(:,1:Y/2)))-mean(mean(CBF_SPD(:,Y/2+1:end)))))/sigma;
PCNR_SPD(i,2) = max(abs(mean(mean(CBV_SPD(:,1:Y/2)))-mean(mean(CBV_SPD(:,Y/2+1:end)))))/sigma;
PCNR_SPD(i,3) = max(abs(mean(mean(MTT_SPD(:,1:Y/2)))-mean(mean(MTT_SPD(:,Y/2+1:end)))))/sigma;
PCNR_SPD(i,4) = max(abs(mean(mean(TTP_SPD(:,1:Y/2)))-mean(mean(TTP_SPD(:,Y/2+1:end)))))/sigma;


% Method 5: Time-Intensity Profile Similarity (TIPS) bilateral filtering 
w     = 5;       % bilateral filter half-width
sigma_tips = [5 sigma]; % bilateral filter standard deviations

Cn_TIPS = permute(pct_tips(permute(Cn,[2 3 1]),w,sigma_tips),[3 1 2]);
RIF_TIPS_bSVD = pct_bsvd(Cn_TIPS,AIF,1,lambda,m,ones(X,Y));
CBF_TIPS_bSVD = squeeze(max(RIF_TIPS_bSVD));
MTT_TIPS_bSVD = squeeze(sum(RIF_TIPS_bSVD))./CBF_TIPS_bSVD;
CBV_TIPS_bSVD = squeeze(sum(RIF_TIPS_bSVD))/60;
TTP_TIPS_bSVD = pct_ttp(pct_tec(AIF,RIF_TIPS_bSVD),1);
PCNR_TIPS_bSVD(i,1) = max(abs(mean(mean(CBF_TIPS_bSVD(:,1:Y/2)))-mean(mean(CBF_TIPS_bSVD(:,Y/2+1:end)))))/sigma;
PCNR_TIPS_bSVD(i,2) = max(abs(mean(mean(CBV_TIPS_bSVD(:,1:Y/2)))-mean(mean(CBV_TIPS_bSVD(:,Y/2+1:end)))))/sigma;
PCNR_TIPS_bSVD(i,3) = max(abs(mean(mean(MTT_TIPS_bSVD(:,1:Y/2)))-mean(mean(MTT_TIPS_bSVD(:,Y/2+1:end)))))/sigma;
PCNR_TIPS_bSVD(i,4) = max(abs(mean(mean(TTP_TIPS_bSVD(:,1:Y/2)))-mean(mean(TTP_TIPS_bSVD(:,Y/2+1:end)))))/sigma;

% Method 6: Tensor Total Variation on RIF (TTV)

% Parameters for TV-based Methods
% input.reg = [2 0.5]; % TV weight for PCNR = 1 [2 0.5]
input.reg = [10 2.5]; % TV weight for PCNR=0.2 [10 1]
input.n1=X;input.n2=Y; input.nt = T;
input.maxitr=500;
input.L=1;input.num=1;
input.l=0; input.u=inf;
input.no = 10;
input.tt = T;
% input.dim = [T X Y];
Ca = pct_circ(AIF,[],[],1); % block-circulant version
input.A = Ca;
input.b = reshape(Cn,T,[]);

out = TTV_FCSA(input);
RIF_TTV = reshape(out.y,[T,X,Y]);
% [RIF_TTV,funv] = mrp_ttv(input);
CBF_TTV = squeeze(max(RIF_TTV));
MTT_TTV = squeeze(sum(RIF_TTV))./CBF_TTV;
CBV_TTV = squeeze(sum(RIF_TTV))/60;
TTP_TTV = pct_ttp(pct_tec(AIF,RIF_TTV),1);
PCNR_TTV(i,1) = max(abs(mean(CBF_TTV(:,1:Y/2))-mean(CBF_TTV(:,Y/2+1:end))))/sigma;
PCNR_TTV(i,2) = max(abs(mean(CBV_TTV(:,1:Y/2))-mean(CBV_TTV(:,Y/2+1:end))))/sigma;
PCNR_TTV(i,3) = max(abs(mean(MTT_TTV(:,1:Y/2))-mean(MTT_TTV(:,Y/2+1:end))))/sigma;
PCNR_TTV(i,4) = max(abs(mean(TTP_TTV(:,1:Y/2))-mean(TTP_TTV(:,Y/2+1:end))))/sigma;


% TIPS + TTV
input.b = reshape(Cn_TIPS,T,[]);

out = TTV_FCSA(input);
RIF_TIPS_TTV = reshape(out.y,[T,X,Y]);
% [RIF_TIPS_TTV,funv] = mrp_ttv(input);
CBF_TIPS_TTV = squeeze(max(RIF_TIPS_TTV));
MTT_TIPS_TTV = squeeze(sum(RIF_TIPS_TTV))./CBF_TIPS_TTV;
CBV_TIPS_TTV = squeeze(sum(RIF_TIPS_TTV))/60;
TTP_TIPS_TTV = pct_ttp(pct_tec(AIF,RIF_TIPS_TTV),1);
PCNR_TIPS_TTV(i,1) = max(abs(mean(mean(CBF_TIPS_TTV(:,1:Y/2)))-mean(mean(CBF_TIPS_TTV(:,Y/2+1:end)))))/sigma;
PCNR_TIPS_TTV(i,2) = max(abs(mean(mean(CBV_TIPS_TTV(:,1:Y/2)))-mean(mean(CBV_TIPS_TTV(:,Y/2+1:end)))))/sigma;
PCNR_TIPS_TTV(i,3) = max(abs(mean(mean(MTT_TIPS_TTV(:,1:Y/2)))-mean(mean(MTT_TIPS_TTV(:,Y/2+1:end)))))/sigma;
PCNR_TIPS_TTV(i,4) = max(abs(mean(mean(TTP_TIPS_TTV(:,1:Y/2)))-mean(mean(TTP_TIPS_TTV(:,Y/2+1:end)))))/sigma;

end

%% Show results
fs = 20;
figure; 
cbf_max = 100;
subplot(491); imshow(CBF_gt,[0 cbf_max]); set(gca,'FontSize',fs); title('Reference'); 
text(-fs,fs,'CBF','FontSize',fs);
subplot(492); imshow(CBF_sSVD,[0 cbf_max]); set(gca,'FontSize',fs); title('sSVD'); 
subplot(493); imshow(CBF_bSVD,[0 cbf_max]); set(gca,'FontSize',fs); title('bSVD'); 
subplot(494); imshow(CBF_TIPS_bSVD,[0 cbf_max]); set(gca,'FontSize',fs); title('TIPS+bSVD');
subplot(495); imshow(CBF_tikh,[0 cbf_max]); set(gca,'FontSize',fs); title('Tikhonov'); 
subplot(496); imshow(CBF_SPD,[0 cbf_max]); set(gca,'FontSize',fs); title('SPD');
subplot(497); imshow(CBF_TTV,[0 cbf_max]); set(gca,'FontSize',fs); title('TTV');
subplot(498); imshow(CBF_TIPS_TTV,[0 cbf_max]); set(gca,'FontSize',fs); title('TIPS+TTV');
subplot(499); imshow([],[0 cbf_max]);
c=colorbar; set(gca,'FontSize',fs);
x1=get(gca,'position');
x=get(c,'Position');
x(1)=x1(1);x(2)=x1(2);x(3)=0.015;x(4)=0.15;
set(c,'Position',x)
set(gca,'position',x1)

cbv_max = 10;
subplot(4,9,10); imshow(CBV_gt,[0 cbv_max]); 
text(-20,20,'CBV','FontSize',20);
subplot(4,9,11); imshow(CBV_sSVD,[0 cbv_max]);
subplot(4,9,12); imshow(CBV_bSVD,[0 cbv_max]);
subplot(4,9,13); imshow(CBV_TIPS_bSVD,[0 cbv_max]); 
subplot(4,9,14); imshow(CBV_tikh,[0 cbv_max]); 
subplot(4,9,15); imshow(CBV_SPD,[0 cbv_max]); 
subplot(4,9,16); imshow(CBV_TTV,[0 cbv_max]); 
subplot(4,9,17); imshow(CBV_TIPS_TTV,[0 cbv_max]); 
subplot(4,9,18); imshow([],[0 cbv_max]);
c=colorbar; set(gca,'FontSize',fs);
x1=get(gca,'position');
x=get(c,'Position');
x(1)=x1(1);x(2)=x1(2);x(3)=0.015;x(4)=0.15;
set(c,'Position',x);
set(gca,'position',x1);

mtt_max = 8;
subplot(4,9,19); imshow(MTT_gt,[0 mtt_max]);
text(-20,20,'MTT','FontSize',20);
subplot(4,9,20); imshow(MTT_sSVD,[0 mtt_max]); 
subplot(4,9,21); imshow(MTT_bSVD,[0 mtt_max]);
subplot(4,9,22); imshow(MTT_TIPS_bSVD,[0 mtt_max]);
subplot(4,9,23); imshow(MTT_tikh,[0 mtt_max]); 
subplot(4,9,24); imshow(MTT_SPD,[0 mtt_max]);
subplot(4,9,25); imshow(MTT_TTV,[0 mtt_max]);
subplot(4,9,26); imshow(MTT_TIPS_TTV,[0 mtt_max]);
subplot(4,9,27); imshow([],[0 mtt_max]);
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

% Save figures
figpath = '../../../Figures/';
w = 25; h = 10;
set(gcf, 'papersize', [w h]);
set(gcf, 'paperposition', [0 0 w h]);
print(gcf,fullfile(figpath,'pdf',['pcnr_' num2str(PCNR) '.pdf']),'-dpdf');
saveas(gcf,fullfile(figpath,'fig',['pcnr_' num2str(PCNR) '.fig']));

% %% Plot Estiamted PCNR vs True PCNR
% fs = 30; % font size for x and y axises
% lw = 4;
% ms = 10;
% spec = {'+-b','d--k','vc','o-.g','.-.m','s-r','^y'};
% 
% % CBF
% h1 = figure;
% line(0:PCNRs(end),0:PCNRs(end),'Color',[0.5 0.5 0.5],'LineWidth',5);
% set(gca,'XTick',0:5:30); xlim([0 max(PCNRs(:,1))]);
% hold on;
% plot(PCNRs(:,1),PCNR_sSVD(:,1), spec{1},'MarkerSize',ms,'MarkerFaceColor',spec{1}(end),'LineWidth',lw);
% plot(PCNRs(:,1),PCNR_bSVD(:,1), spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
% plot(PCNRs(:,1),PCNR_TIPS_bSVD(:,1), spec{3},'MarkerSize',ms,'MarkerFaceColor',spec{3}(end),'LineWidth',lw);
% plot(PCNRs(:,1),PCNR_tikh(:,1), spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
% plot(PCNRs(:,1),PCNR_SPD(:,1), spec{5},'MarkerSize',ms,'MarkerFaceColor',spec{5}(end),'LineWidth',lw);
% plot(PCNRs(:,1),PCNR_TTV(:,1), spec{6},'MarkerSize',ms,'MarkerFaceColor',spec{6}(end),'LineWidth',lw);
% plot(PCNRs(:,1),PCNR_TIPS_TTV(:,1), spec{7},'MarkerSize',ms,'MarkerFaceColor',spec{7}(end),'LineWidth',lw);
% set(gca,'FontSize',fs);
% xlabel('True PCNR');
% ylabel('Estimated PCNR');
% legend('Reference','sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');
% 
% % CBV
% h2 = figure;
% line(0:PCNRs(end),0:PCNRs(end),'Color',[0.5 0.5 0.5],'LineWidth',5);
% set(gca,'XTick',0:5:30); xlim([0 max(PCNRs(:,2))]);
% hold on;
% plot(PCNRs(:,2),PCNR_sSVD(:,2), spec{1},'MarkerSize',ms,'MarkerFaceColor',spec{1}(end),'LineWidth',lw);
% plot(PCNRs(:,2),PCNR_bSVD(:,2), spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
% plot(PCNRs(:,2),PCNR_TIPS_bSVD(:,2), spec{3},'MarkerSize',ms,'MarkerFaceColor',spec{3}(end),'LineWidth',lw);
% plot(PCNRs(:,2),PCNR_tikh(:,2), spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
% plot(PCNRs(:,2),PCNR_SPD(:,2), spec{5},'MarkerSize',ms,'MarkerFaceColor',spec{5}(end),'LineWidth',lw);
% plot(PCNRs(:,2),PCNR_TTV(:,2), spec{6},'MarkerSize',ms,'MarkerFaceColor',spec{6}(end),'LineWidth',lw);
% plot(PCNRs(:,2),PCNR_TIPS_TTV(:,2), spec{7},'MarkerSize',ms,'MarkerFaceColor',spec{7}(end),'LineWidth',lw);
% set(gca,'FontSize',fs);
% xlabel('True PCNR');
% ylabel('Estimated PCNR');
% legend('Reference','sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');
% 
% 
% h3 = figure;
% line(0:PCNRs(end),0:PCNRs(end),'Color',[0.5 0.5 0.5],'LineWidth',5);
% set(gca,'XTick',0:5:30); xlim([0 max(PCNRs(:,3))]);
% hold on;
% plot(PCNRs(:,3),PCNR_sSVD(:,3), spec{1},'MarkerSize',ms,'MarkerFaceColor',spec{1}(end),'LineWidth',lw);
% plot(PCNRs(:,3),PCNR_bSVD(:,3), spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
% plot(PCNRs(:,3),PCNR_TIPS_bSVD(:,3), spec{3},'MarkerSize',ms,'MarkerFaceColor',spec{3}(end),'LineWidth',lw);
% plot(PCNRs(:,3),PCNR_tikh(:,3), spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
% plot(PCNRs(:,3),PCNR_SPD(:,3), spec{5},'MarkerSize',ms,'MarkerFaceColor',spec{5}(end),'LineWidth',lw);
% plot(PCNRs(:,3),PCNR_TTV(:,3), spec{6},'MarkerSize',ms,'MarkerFaceColor',spec{6}(end),'LineWidth',lw);
% plot(PCNRs(:,3),PCNR_TIPS_TTV(:,3), spec{7},'MarkerSize',ms,'MarkerFaceColor',spec{7}(end),'LineWidth',lw);
% set(gca,'FontSize',fs);
% xlabel('True PCNR');
% ylabel('Estimated PCNR');
% legend('Reference','sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');
% 
% 
% h1 = figure;
% line(0:PCNRs(end),0:PCNRs(end),'Color',[0.5 0.5 0.5],'LineWidth',5);
% set(gca,'XTick',0:5:30); xlim([0 max(PCNRs(:,4))]);
% hold on;
% plot(PCNRs(:,4),PCNR_sSVD(:,4), spec{1},'MarkerSize',ms,'MarkerFaceColor',spec{1}(end),'LineWidth',lw);
% plot(PCNRs(:,4),PCNR_bSVD(:,4), spec{2},'MarkerSize',ms,'MarkerFaceColor',spec{2}(end),'LineWidth',lw);
% plot(PCNRs(:,4),PCNR_TIPS_bSVD(:,4), spec{3},'MarkerSize',ms,'MarkerFaceColor',spec{3}(end),'LineWidth',lw);
% plot(PCNRs(:,4),PCNR_tikh(:,4), spec{4},'MarkerSize',ms,'MarkerFaceColor',spec{4}(end),'LineWidth',lw);
% plot(PCNRs(:,4),PCNR_SPD(:,4), spec{5},'MarkerSize',ms,'MarkerFaceColor',spec{5}(end),'LineWidth',lw);
% plot(PCNRs(:,4),PCNR_TTV(:,4), spec{6},'MarkerSize',ms,'MarkerFaceColor',spec{6}(end),'LineWidth',lw);
% plot(PCNRs(:,4),PCNR_TIPS_TTV(:,4), spec{7},'MarkerSize',ms,'MarkerFaceColor',spec{7}(end),'LineWidth',lw);
% set(gca,'FontSize',fs);
% xlabel('True PCNR');
% ylabel('Estimated PCNR');
% legend('Reference','sSVD','bSVD','TIPS+bSVD','Tikhonov','SPD','TTV','TIPS+TTV','Location','NorthWest');
% 
