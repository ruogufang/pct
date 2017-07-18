% Experiment on Digital Head Phantom
%
%   Ruogu Fang 12/01/2014
%   Smadt Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University


close all; clear; clc;

% Set paths
addpath(genpath('toolbox')); % Include toolbox package

%%%%%%%%%% Setting parameters %%%%%%%%%%%%%%%%%%
mA = 190; % tube current-exposure time product
dt = 1; % downsampling rate in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tid_all = tic;

% classes of tissue
GM = 1; WM = 2; GMR = 3; GMSR = 4; WMR = 5; WMSR = 6; AI = 7; VO = 8; CSF = 9; HEGM = 10; HEWM = 11;

% Parameters
m = 2; % block-circulant Ca extended for m times
lambda = 0.15; % cTSVD truncation parameter for singular values
Tn = 30; % truncated time in TTV
rho = 1; % tissue density g/mL

% Load data
imgname = 'phantom';
load(imgname);
% Data format:
%   V: CTP data [T x X x Y]
[T,X,Y]=size(ctp);
PRE = 1; POST = T;
V = ctp - permute(repmat(baseline,[1 1 T]),[3 1 2]);
AIF = aif(PRE:POST); % T x 1 vector

load acf;

% noisy level
mA0 = 190;
sigma = pct_mA2sigma(mA,mA0);

% Remove negative values
V(V<0) = 0;

% Find Brain Mask
Mask = cbf>0;

% Add correlated Gaussian (spectral) noise to simulate low-dose
Vn = pct_noise(permute(V,[2 3 1]),acf,sigma,'s',Mask);
Vn = permute(Vn,[3 1 2]);

% PCT preprocess
% Whole Brain CBF
C = pct_preprocess(V, PRE, POST); % noiseless data
Cn = pct_preprocess(Vn, PRE, POST); % noisy data
% Cn = C;

% Crop data to Region of Interest
if exist('x0','var')
    C = C(:,y0:y0+h-1,x0:x0+w-1);
    Cn = Cn(:,y0:y0+h-1,x0:x0+w-1);
    Mask = Mask(y0:y0+h-1,x0:x0+w-1);
    brain = brain(y0:y0+h-1,x0:x0+w-1);
end

% Find the minimum bounding box
bb = pct_minBoundingBox(Mask);
C = C(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
Cn = Cn(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
Mask = Mask(bb(1,1):bb(1,2),bb(2,1):bb(2,3));
brain = brain(bb(1,1):bb(1,2),bb(2,1):bb(2,3));

cbf = cbf(bb(1,1):bb(1,2),bb(2,1):bb(2,3));
cbv = cbv(bb(1,1):bb(1,2),bb(2,1):bb(2,3));
mtt = mtt(bb(1,1):bb(1,2),bb(2,1):bb(2,3));
ttp = ttp(bb(1,1):bb(1,2),bb(2,1):bb(2,3));

clear V Vn VOF PRE POST B acf

% Method 1: Standard Truncated Singular Value Decomposition (sSVD)
tic;
RIF_sSVD = pct_ssvd(Cn,AIF,dt,lambda,Mask);
RIF_sSVD(RIF_sSVD<0) = 0;
CBF_sSVD = pct_cbf(RIF_sSVD,rho);
CBV_sSVD = pct_cbv(RIF_sSVD,rho);
MTT_sSVD = pct_mtt(RIF_sSVD);
TTP_sSVD = pct_ttp(pct_tec(AIF,RIF_sSVD),1);
t_sSVD = toc;

% Method 2: block-circulant SVD (bSVD)
tic;
RIF_bSVD = pct_bsvd(Cn,AIF,dt,lambda,m,Mask);
RIF_bSVD(RIF_bSVD<0) = 0;
CBF_bSVD = pct_cbf(RIF_bSVD,rho);
CBV_bSVD = pct_cbv(RIF_bSVD,rho);
MTT_bSVD = pct_mtt(RIF_bSVD);
TTP_bSVD = pct_ttp(pct_tec(AIF,RIF_bSVD),1);
t_bSVD = toc;

% Method 3: Tikhonov
tic;
RIF_tikh = pct_tikh(Cn,AIF,dt,lambda,1,Mask);
RIF_tikh(RIF_tikh<0) = 0;
CBF_tikh = pct_cbf(RIF_tikh,rho);
CBV_tikh = pct_cbv(RIF_tikh,rho);
MTT_tikh = pct_mtt(RIF_tikh);
TTP_tikh = pct_ttp(pct_tec(AIF,RIF_tikh),1);
t_tikh = toc;

% Method 4: TIPS+bSVD
tic;

% Set bilateral filter parameters.
w     = 5;       % bilateral filter half-width
sigma = [3 0.1]; % bilateral filter standard deviations

Cn_TIPS = pct_tips(Cn,w,sigma);
RIF_TIPS_bSVD = pct_ssvd(Cn_TIPS,AIF,dt,lambda,Mask);
RIF_TIPS_bSVD(RIF_TIPS_bSVD<0) = 0;
CBF_TIPS_bSVD = pct_cbf(RIF_TIPS_bSVD,rho);
CBV_TIPS_bSVD = pct_cbv(RIF_TIPS_bSVD,rho);
MTT_TIPS_bSVD = pct_mtt(RIF_TIPS_bSVD);
TTP_TIPS_bSVD = pct_ttp(pct_tec(AIF,RIF_TIPS_bSVD),1);
t_TIPS_bSVD = toc;

%% Method 5: Sparse Perfusion Deconvolution (SPD)
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
params.lambda = 1e-4; % weight of the sparsity term in SPAMS package
params.mixture = 10; % weight of input signal
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


%% Method 6: Total Variation Regularization (TTV)
% parameters
scale = 6000;
m = 1;
[T,h,w]=size(Cn);
% input.reg = 1e-4; % when dt=1
input.reg = [1e-4 1e-2] * 1e-4; % no noise
% input.reg = 1e-10;
input.maxitr=500;
input.l=-inf; input.u=inf;
input.no = 30; % max number of iterations allowed
[Ca,Cc] = pct_circ(AIF,Cn,[],m); % block-circulant version
Tc = T*m;
input.A=Ca;
input.b=reshape(Cc,Tc,[]);
input.tt = Tn;
input.n1=h; input.n2=w; input.nt = Tc;
input.dim = [Tc h w];
input.dt = dt;
input.T = T;
input.tv = 'l1';

% TTV
tic;
out_ttv = TTV_FCSA(input);
RIF_TTV = reshape(out_ttv.y,[Tc h w]);
% RIF_TTV = pct_ttv(input);
RIF_TTV(RIF_TTV<0) = 0;
t_TTV = toc;
CBF_TTV = pct_cbf(RIF_TTV,rho,Mask);
CBV_TTV = pct_cbv(RIF_TTV,rho,Mask);
MTT_TTV = pct_mtt(RIF_TTV,Mask);
TTP_TTV = pct_ttp(pct_tec(AIF,RIF_TTV),1,Mask);
fprintf(1,'TTV: time = %.2fsec, Iteration Num = %d\n',t_TTV, input.no);

% Method 7: TIPS+TTV
input.b = reshape(Cn_TIPS,T,[]);
out_tips_ttv = TTV_FCSA(input);
RIF_TIPS_TTV = reshape(out_tips_ttv.y,[Tc h w]);
% [RIF_TIPS_TTV,funv] = pct_ttv(input);
RIF_TIPS_TTV(RIF_TIPS_TTV<0) = 0;
CBF_TIPS_TTV = pct_cbf(RIF_TIPS_TTV,rho,Mask);
CBV_TIPS_TTV = pct_cbv(RIF_TIPS_TTV,rho,Mask);
MTT_TIPS_TTV = pct_mtt(RIF_TIPS_TTV,Mask);
TTP_TIPS_TTV = pct_ttp(pct_tec(AIF,RIF_TIPS_TTV),1,Mask);

% Method 8: Initialization with baseline
input.b = reshape(Cn,T,[]);
input.yr = reshape(RIF_bSVD,T,[]);
out_ttv_init = TTV_FCSA_init(input);
RIF_TTV_init = reshape(out_ttv_init.y,[Tc h w]);
RIF_TTV_init(RIF_TTV_init<0) = 0;
CBF_TTV_init = squeeze(max(RIF_TTV_init));
MTT_TTV_init = squeeze(sum(RIF_TTV_init))./CBF_TTV_init;
CBV_TTV_init = squeeze(sum(RIF_TTV_init))/60;
TTP_TTV_init = pct_ttp(pct_tec(AIF,RIF_TTV_init),1);

%% Plot results
close all;

% % RIFs
% x = 106; y = 40;
% figure;
% plot(RIF_sSVD(:,x,y),'b');
% hold on;
% plot(RIF_bSVD(:,x,y),'g');
% plot(RIF_TTV(:,x,y),'r');


fs = 25;

CBF = [cbf CBF_sSVD CBF_bSVD CBF_TIPS_bSVD CBF_tikh CBF_SPD CBF_TTV CBF_TIPS_TTV];
CBV = [cbv CBV_sSVD CBV_bSVD CBV_TIPS_bSVD CBV_tikh CBV_SPD CBV_TTV CBV_TIPS_TTV];
MTT = [mtt MTT_sSVD MTT_bSVD MTT_TIPS_bSVD MTT_tikh MTT_SPD MTT_TTV MTT_TIPS_TTV];
TTP = [ttp TTP_sSVD TTP_bSVD TTP_TIPS_bSVD TTP_tikh TTP_SPD TTP_TTV TTP_TIPS_TTV];
% 
% figure; ctshow(CBF,[],[0 80]); colorbar; set(gca,'FontSize',fs);
% figure; ctshow(CBV,[],[0 4]); colorbar; set(gca,'FontSize',fs);
% figure; ctshow(MTT,[],[0 15]); colorbar; set(gca,'FontSize',fs);
% figure; ctshow(TTP,[],[0 30]); colorbar; set(gca,'FontSize',fs);

% figure;
% subplot(411); ctshow(CBF,[],[0 80]); h=text(-20,120,'CBF','FontSize',20); set(h,'rotation',90); c=colorbar; set(gca,'FontSize',fs);
% subplot(412); ctshow(CBV,[],[0 4]); h=text(-20,120,'CBV','FontSize',20); set(h,'rotation',90);colorbar; set(gca,'FontSize',fs);
% subplot(413); ctshow(MTT,[],[0 15]); h=text(-20,120,'MTT','FontSize',20); set(h,'rotation',90);colorbar; set(gca,'FontSize',fs);
% subplot(414); ctshow(TTP,[],[0 30]); h=text(-20,120,'TTP','FontSize',20); set(h,'rotation',90);colorbar; set(gca,'FontSize',fs);

[X,Y] = size(cbf);
figure;
ha = tight_subplot(4,1,[.02 .03],[.01 .05],[.01 .01]);
axes(ha(1)); ctshow(CBF,[],[0 80]); h=text(-15,120,'CBF','FontSize',fs); set(h,'rotation',90); c=colorbar; set(gca,'FontSize',fs);
text(25,-20,'Reference','FontSize',fs); text(25+Y,-20,'sSVD','FontSize',fs); text(25+2*Y,-20,'bSVD','FontSize',fs); text(25+3*Y,-20,'TIPS+bSVD','FontSize',fs); text(25+4*Y,-20,'Tikhonov','FontSize',fs); text(25+5*Y,-20,'SPD','FontSize',fs); text(25+6*Y,-20,'TTV','FontSize',fs);text(25+7*Y,-20,'TIPS+TTV','FontSize',fs);
axes(ha(2)); ctshow(CBV,[],[0 4]); h=text(-15,120,'CBV','FontSize',fs); set(h,'rotation',90);colorbar; set(gca,'FontSize',fs);
axes(ha(3)); ctshow(MTT,[],[0 15]); h=text(-15,120,'MTT','FontSize',fs); set(h,'rotation',90);colorbar; set(gca,'FontSize',fs);
axes(ha(4)); ctshow(TTP,[],[0 30]); h=text(-15,120,'TTP','FontSize',fs); set(h,'rotation',90);colorbar; set(gca,'FontSize',fs);

% Save figures
figpath = '../../../Figures/';
w = 25; h = 15;
set(gcf, 'papersize', [w h]);
set(gcf, 'paperposition', [0 0 w h]);
print(gcf,fullfile(figpath,'pdf','phantom.pdf'),'-dpdf');
saveas(gcf,fullfile(figpath,'fig','phantom.fig'));

% Convergence with different initialization
figure;
plot(out_ttv.funv,'b','LineWidth',2);
hold on;
plot(out_ttv_init.funv,'-^k','LineWidth',2);
plot(out_tips_ttv.funv,'--r','LineWidth',2);
legend('TTV init with zero','TTV init with bSVD','TIPS+TTV');
set(gca,'FontSize',fs);

figure;
plot(out_ttv.funv(1:5),'b','LineWidth',2);
hold on;
plot(out_ttv_init.funv(1:5),'-^k','LineWidth',2);
plot(out_tips_ttv.funv(1:5),'--r','LineWidth',2);
legend('TTV init with zero','TTV init with bSVD','TIPS+TTV');
set(gca,'FontSize',fs);

%% RMSE and Lin's Concordance

rmse_CBF_sSVD = pct_rmse(CBF_sSVD,cbf,Mask);
rmse_CBF_bSVD = pct_rmse(CBF_bSVD,cbf,Mask);
rmse_CBF_TIPS_bSVD = pct_rmse(CBF_TIPS_bSVD,cbf,Mask);
rmse_CBF_tikh = pct_rmse(CBF_tikh,cbf,Mask);
rmse_CBF_SPD = pct_rmse(CBF_SPD,cbf,Mask);
rmse_CBF_TTV = pct_rmse(CBF_TTV,cbf,Mask);
rmse_CBF_TIPS_TTV = pct_rmse(CBF_TIPS_TTV,cbf,Mask);

rmse_CBV_sSVD = pct_rmse(CBV_sSVD,cbv,Mask);
rmse_CBV_bSVD = pct_rmse(CBV_bSVD,cbv,Mask);
rmse_CBV_TIPS_bSVD = pct_rmse(CBV_TIPS_bSVD,cbv,Mask);
rmse_CBV_tikh = pct_rmse(CBV_tikh,cbv,Mask);
rmse_CBV_SPD = pct_rmse(CBV_SPD,cbv,Mask);
rmse_CBV_TTV = pct_rmse(CBV_TTV,cbv,Mask);
rmse_CBV_TIPS_TTV = pct_rmse(CBV_TIPS_TTV,cbv,Mask);

rmse_MTT_sSVD = pct_rmse(MTT_sSVD,mtt,Mask);
rmse_MTT_bSVD = pct_rmse(MTT_bSVD,mtt,Mask);
rmse_MTT_TIPS_bSVD = pct_rmse(MTT_TIPS_bSVD,mtt,Mask);
rmse_MTT_tikh = pct_rmse(MTT_tikh,mtt,Mask);
rmse_MTT_SPD = pct_rmse(MTT_SPD,mtt,Mask);
rmse_MTT_TTV = pct_rmse(MTT_TTV,mtt,Mask);
rmse_MTT_TIPS_TTV = pct_rmse(MTT_TIPS_TTV,mtt,Mask);

rmse_TTP_sSVD = pct_rmse(TTP_sSVD,ttp,Mask);
rmse_TTP_bSVD = pct_rmse(TTP_bSVD,ttp,Mask);
rmse_TTP_TIPS_bSVD = pct_rmse(TTP_TIPS_bSVD,ttp,Mask);
rmse_TTP_tikh = pct_rmse(TTP_tikh,ttp,Mask);
rmse_TTP_SPD = pct_rmse(TTP_SPD,ttp,Mask);
rmse_TTP_TTV = pct_rmse(TTP_TTV,ttp,Mask);
rmse_TTP_TIPS_TTV = pct_rmse(TTP_TIPS_TTV,ttp,Mask);

% Tissue
rmse_CBF_sSVD_tissue = zeros(1,6);
rmse_CBF_bSVD_tissue = zeros(1,6);
rmse_CBF_TIPS_bSVD_tissue = zeros(1,6);
rmse_CBF_tikh_tissue = zeros(1,6);
rmse_CBF_SPD_tissue = zeros(1,6);
rmse_CBF_TTV_tissue = zeros(1,6);
rmse_CBF_TIPS_TTV_tissue = zeros(1,6);

rmse_CBV_sSVD_tissue = zeros(1,6);
rmse_CBV_bSVD_tissue = zeros(1,6);
rmse_CBV_TIPS_bSVD_tissue = zeros(1,6);
rmse_CBV_tikh_tissue = zeros(1,6);
rmse_CBV_SPD_tissue = zeros(1,6);
rmse_CBV_TTV_tissue = zeros(1,6);
rmse_CBV_TIPS_TTV_tissue = zeros(1,6);

rmse_MTT_sSVD_tissue = zeros(1,6);
rmse_MTT_bSVD_tissue = zeros(1,6);
rmse_MTT_TIPS_bSVD_tissue = zeros(1,6);
rmse_MTT_tikh_tissue = zeros(1,6);
rmse_MTT_SPD_tissue = zeros(1,6);
rmse_MTT_TTV_tissue = zeros(1,6);
rmse_MTT_TIPS_TTV_tissue = zeros(1,6);

rmse_TTP_sSVD_tissue = zeros(1,6);
rmse_TTP_bSVD_tissue = zeros(1,6);
rmse_TTP_TIPS_bSVD_tissue = zeros(1,6);
rmse_TTP_tikh_tissue = zeros(1,6);
rmse_TTP_SPD_tissue = zeros(1,6);
rmse_TTP_TTV_tissue = zeros(1,6);
rmse_TTP_TIPS_TTV_tissue = zeros(1,6);

for i = 1:6
    mask = (brain==i);
    rmse_CBF_sSVD_tissue(i) = pct_rmse(CBF_sSVD,cbf,mask);
    rmse_CBF_bSVD_tissue(i) = pct_rmse(CBF_bSVD,cbf,mask);
    rmse_CBF_TIPS_bSVD_tissue(i) = pct_rmse(CBF_TIPS_bSVD,cbf,mask);
    rmse_CBF_tikh_tissue(i) = pct_rmse(CBF_tikh,cbf,mask);
    rmse_CBF_SPD_tissue(i) = pct_rmse(CBF_tikh,cbf,mask);
    rmse_CBF_TTV_tissue(i) = pct_rmse(CBF_TTV,cbf,mask);
    rmse_CBF_TIPS_TTV_tissue(i) = pct_rmse(CBF_TIPS_TTV,cbf,mask);
    
    rmse_CBV_sSVD_tissue(i) = pct_rmse(CBV_sSVD,cbv,mask);
    rmse_CBV_bSVD_tissue(i) = pct_rmse(CBV_bSVD,cbv,mask);
    rmse_CBV_TIPS_bSVD_tissue(i) = pct_rmse(CBV_TIPS_bSVD,cbv,mask);
    rmse_CBV_tikh_tissue(i) = pct_rmse(CBV_tikh,cbv,mask);
    rmse_CBV_SPD_tissue(i) = pct_rmse(CBV_tikh,cbv,mask);
    rmse_CBV_TTV_tissue(i) = pct_rmse(CBV_TTV,cbv,mask);
    rmse_CBV_TIPS_TTV_tissue(i) = pct_rmse(CBV_TIPS_TTV,cbv,mask);
   
    rmse_MTT_sSVD_tissue(i) = pct_rmse(MTT_sSVD,mtt,mask);
    rmse_MTT_bSVD_tissue(i) = pct_rmse(MTT_bSVD,mtt,mask);
    rmse_MTT_TIPS_bSVD_tissue(i) = pct_rmse(MTT_TIPS_bSVD,mtt,mask);
    rmse_MTT_tikh_tissue(i) = pct_rmse(MTT_tikh,mtt,mask);
    rmse_MTT_SPD_tissue(i) = pct_rmse(MTT_tikh,mtt,mask);
    rmse_MTT_TTV_tissue(i) = pct_rmse(MTT_TTV,mtt,mask);
    rmse_MTT_TIPS_TTV_tissue(i) = pct_rmse(MTT_TIPS_TTV,mtt,mask);
    
    rmse_TTP_sSVD_tissue(i) = pct_rmse(TTP_sSVD,ttp,mask);
    rmse_TTP_bSVD_tissue(i) = pct_rmse(TTP_bSVD,ttp,mask);
    rmse_TTP_TIPS_bSVD_tissue(i) = pct_rmse(TTP_TIPS_bSVD,ttp,mask);
    rmse_TTP_tikh_tissue(i) = pct_rmse(TTP_tikh,ttp,mask);
    rmse_TTP_SPD_tissue(i) = pct_rmse(TTP_tikh,ttp,mask);
    rmse_TTP_TTV_tissue(i) = pct_rmse(TTP_TTV,ttp,mask);
    rmse_TTP_TIPS_TTV_tissue(i) = pct_rmse(TTP_TIPS_TTV,ttp,mask);
    
end

rmse = [rmse_CBF_sSVD_tissue rmse_CBF_sSVD;rmse_CBF_bSVD_tissue rmse_CBF_bSVD; rmse_CBF_TIPS_bSVD_tissue rmse_CBF_TIPS_bSVD; rmse_CBF_tikh_tissue rmse_CBF_tikh; rmse_CBF_SPD_tissue rmse_CBF_SPD; rmse_CBF_TTV_tissue rmse_CBF_TTV; rmse_CBF_TIPS_TTV_tissue rmse_CBF_TIPS_TTV;
    rmse_CBV_sSVD_tissue rmse_CBV_sSVD;rmse_CBV_bSVD_tissue rmse_CBV_bSVD; rmse_CBV_TIPS_bSVD_tissue rmse_CBV_TIPS_bSVD; rmse_CBV_tikh_tissue rmse_CBV_tikh; rmse_CBV_SPD_tissue rmse_CBV_SPD; rmse_CBV_TTV_tissue rmse_CBV_TTV; rmse_CBV_TIPS_TTV_tissue rmse_CBV_TIPS_TTV;
    rmse_MTT_sSVD_tissue rmse_MTT_sSVD;rmse_MTT_bSVD_tissue rmse_MTT_bSVD; rmse_MTT_TIPS_bSVD_tissue rmse_MTT_TIPS_bSVD; rmse_MTT_tikh_tissue rmse_MTT_tikh; rmse_MTT_SPD_tissue rmse_MTT_SPD; rmse_MTT_TTV_tissue rmse_MTT_TTV; rmse_MTT_TIPS_TTV_tissue rmse_MTT_TIPS_TTV;
    rmse_TTP_sSVD_tissue rmse_TTP_sSVD;rmse_TTP_bSVD_tissue rmse_TTP_bSVD; rmse_TTP_TIPS_bSVD_tissue rmse_TTP_TIPS_bSVD; rmse_TTP_tikh_tissue rmse_TTP_tikh; rmse_TTP_SPD_tissue rmse_TTP_SPD; rmse_TTP_TTV_tissue rmse_TTP_TTV; rmse_TTP_TIPS_TTV_tissue rmse_TTP_TIPS_TTV;
    ];

rmse = round(rmse/0.01)*0.01;


% Lin's Coefficient
lin_CBF_sSVD = pct_lincon(CBF_sSVD,cbf,Mask);
lin_CBF_bSVD = pct_lincon(CBF_bSVD,cbf,Mask);
lin_CBF_TIPS_bSVD = pct_lincon(CBF_TIPS_bSVD,cbf,Mask);
lin_CBF_tikh = pct_lincon(CBF_tikh,cbf,Mask);
lin_CBF_SPD = pct_lincon(CBF_SPD,cbf,Mask);
lin_CBF_TTV = pct_lincon(CBF_TTV,cbf,Mask);
lin_CBF_TIPS_TTV = pct_lincon(CBF_TIPS_TTV,cbf,Mask);

lin_CBV_sSVD = pct_lincon(CBV_sSVD,cbv,Mask);
lin_CBV_bSVD = pct_lincon(CBV_bSVD,cbv,Mask);
lin_CBV_TIPS_bSVD = pct_lincon(CBV_TIPS_bSVD,cbv,Mask);
lin_CBV_tikh = pct_lincon(CBV_tikh,cbv,Mask);
lin_CBV_SPD = pct_lincon(CBV_SPD,cbv,Mask);
lin_CBV_TTV = pct_lincon(CBV_TTV,cbv,Mask);
lin_CBV_TIPS_TTV = pct_lincon(CBV_TIPS_TTV,cbv,Mask);

lin_MTT_sSVD = pct_lincon(MTT_sSVD,mtt,Mask);
lin_MTT_bSVD = pct_lincon(MTT_bSVD,mtt,Mask);
lin_MTT_TIPS_bSVD = pct_lincon(MTT_TIPS_bSVD,mtt,Mask);
lin_MTT_tikh = pct_lincon(MTT_tikh,mtt,Mask);
lin_MTT_SPD = pct_lincon(MTT_SPD,mtt,Mask);
lin_MTT_TTV = pct_lincon(MTT_TTV,mtt,Mask);
lin_MTT_TIPS_TTV = pct_lincon(MTT_TIPS_TTV,mtt,Mask);

lin_TTP_sSVD = pct_lincon(TTP_sSVD,ttp,Mask);
lin_TTP_bSVD = pct_lincon(TTP_bSVD,ttp,Mask);
lin_TTP_TIPS_bSVD = pct_lincon(TTP_TIPS_bSVD,ttp,Mask);
lin_TTP_tikh = pct_lincon(TTP_tikh,ttp,Mask);
lin_TTP_SPD = pct_lincon(TTP_SPD,ttp,Mask);
lin_TTP_TTV = pct_lincon(TTP_TTV,ttp,Mask);
lin_TTP_TIPS_TTV = pct_lincon(TTP_TIPS_TTV,ttp,Mask);

% Tissue
lin_CBF_sSVD_tissue = zeros(1,6);
lin_CBF_bSVD_tissue = zeros(1,6);
lin_CBF_TIPS_bSVD_tissue = zeros(1,6);
lin_CBF_tikh_tissue = zeros(1,6);
lin_CBF_SPD_tissue = zeros(1,6);
lin_CBF_TTV_tissue = zeros(1,6);
lin_CBF_TIPS_TTV_tissue = zeros(1,6);

lin_CBV_sSVD_tissue = zeros(1,6);
lin_CBV_bSVD_tissue = zeros(1,6);
lin_CBV_TIPS_bSVD_tissue = zeros(1,6);
lin_CBV_tikh_tissue = zeros(1,6);
lin_CBV_SPD_tissue = zeros(1,6);
lin_CBV_TTV_tissue = zeros(1,6);
lin_CBV_TIPS_TTV_tissue = zeros(1,6);

lin_MTT_sSVD_tissue = zeros(1,6);
lin_MTT_bSVD_tissue = zeros(1,6);
lin_MTT_TIPS_bSVD_tissue = zeros(1,6);
lin_MTT_tikh_tissue = zeros(1,6);
lin_MTT_SPD_tissue = zeros(1,6);
lin_MTT_TTV_tissue = zeros(1,6);
lin_MTT_TIPS_TTV_tissue = zeros(1,6);

lin_TTP_sSVD_tissue = zeros(1,6);
lin_TTP_bSVD_tissue = zeros(1,6);
lin_TTP_TIPS_bSVD_tissue = zeros(1,6);
lin_TTP_tikh_tissue = zeros(1,6);
lin_TTP_SPD_tissue = zeros(1,6);
lin_TTP_TTV_tissue = zeros(1,6);
lin_TTP_TIPS_TTV_tissue = zeros(1,6);

for i = 1:6
    mask = (brain==i);
    lin_CBF_sSVD_tissue(i) = pct_lincon(CBF_sSVD,cbf,mask);
    lin_CBF_bSVD_tissue(i) = pct_lincon(CBF_bSVD,cbf,mask);
    lin_CBF_TIPS_bSVD_tissue(i) = pct_lincon(CBF_TIPS_bSVD,cbf,mask);
    lin_CBF_tikh_tissue(i) = pct_lincon(CBF_tikh,cbf,mask);
    lin_CBF_SPD_tissue(i) = pct_lincon(CBF_tikh,cbf,mask);
    lin_CBF_TTV_tissue(i) = pct_lincon(CBF_TTV,cbf,mask);
    lin_CBF_TIPS_TTV_tissue(i) = pct_lincon(CBF_TIPS_TTV,cbf,mask);
    
    lin_CBV_sSVD_tissue(i) = pct_lincon(CBV_sSVD,cbv,mask);
    lin_CBV_bSVD_tissue(i) = pct_lincon(CBV_bSVD,cbv,mask);
    lin_CBV_TIPS_bSVD_tissue(i) = pct_lincon(CBV_TIPS_bSVD,cbv,mask);
    lin_CBV_tikh_tissue(i) = pct_lincon(CBV_tikh,cbv,mask);
    lin_CBV_SPD_tissue(i) = pct_lincon(CBV_tikh,cbv,mask);
    lin_CBV_TTV_tissue(i) = pct_lincon(CBV_TTV,cbv,mask);
    lin_CBV_TIPS_TTV_tissue(i) = pct_lincon(CBV_TIPS_TTV,cbv,mask);
   
    lin_MTT_sSVD_tissue(i) = pct_lincon(MTT_sSVD,mtt,mask);
    lin_MTT_bSVD_tissue(i) = pct_lincon(MTT_bSVD,mtt,mask);
    lin_MTT_TIPS_bSVD_tissue(i) = pct_lincon(MTT_TIPS_bSVD,mtt,mask);
    lin_MTT_tikh_tissue(i) = pct_lincon(MTT_tikh,mtt,mask);
    lin_MTT_SPD_tissue(i) = pct_lincon(MTT_tikh,mtt,mask);
    lin_MTT_TTV_tissue(i) = pct_lincon(MTT_TTV,mtt,mask);
    lin_MTT_TIPS_TTV_tissue(i) = pct_lincon(MTT_TIPS_TTV,mtt,mask);
    
    lin_TTP_sSVD_tissue(i) = pct_lincon(TTP_sSVD,ttp,mask);
    lin_TTP_bSVD_tissue(i) = pct_lincon(TTP_bSVD,ttp,mask);
    lin_TTP_TIPS_bSVD_tissue(i) = pct_lincon(TTP_TIPS_bSVD,ttp,mask);
    lin_TTP_tikh_tissue(i) = pct_lincon(TTP_tikh,ttp,mask);
    lin_TTP_SPD_tissue(i) = pct_lincon(TTP_tikh,ttp,mask);
    lin_TTP_TTV_tissue(i) = pct_lincon(TTP_TTV,ttp,mask);
    lin_TTP_TIPS_TTV_tissue(i) = pct_lincon(TTP_TIPS_TTV,ttp,mask);
    
end

lin = [lin_CBF_sSVD_tissue lin_CBF_sSVD;lin_CBF_bSVD_tissue lin_CBF_bSVD; lin_CBF_TIPS_bSVD_tissue lin_CBF_TIPS_bSVD; lin_CBF_tikh_tissue lin_CBF_tikh; lin_CBF_SPD_tissue lin_CBF_SPD; lin_CBF_TTV_tissue lin_CBF_TTV; lin_CBF_TIPS_TTV_tissue lin_CBF_TIPS_TTV;
    lin_CBV_sSVD_tissue lin_CBV_sSVD;lin_CBV_bSVD_tissue lin_CBV_bSVD; lin_CBV_TIPS_bSVD_tissue lin_CBV_TIPS_bSVD; lin_CBV_tikh_tissue lin_CBV_tikh; lin_CBV_SPD_tissue lin_CBV_SPD; lin_CBV_TTV_tissue lin_CBV_TTV; lin_CBV_TIPS_TTV_tissue lin_CBV_TIPS_TTV;
    lin_MTT_sSVD_tissue lin_MTT_sSVD;lin_MTT_bSVD_tissue lin_MTT_bSVD; lin_MTT_TIPS_bSVD_tissue lin_MTT_TIPS_bSVD; lin_MTT_tikh_tissue lin_MTT_tikh; lin_MTT_SPD_tissue lin_MTT_SPD; lin_MTT_TTV_tissue lin_MTT_TTV; lin_MTT_TIPS_TTV_tissue lin_MTT_TIPS_TTV;
    lin_TTP_sSVD_tissue lin_TTP_sSVD;lin_TTP_bSVD_tissue lin_TTP_bSVD; lin_TTP_TIPS_bSVD_tissue lin_TTP_TIPS_bSVD; lin_TTP_tikh_tissue lin_TTP_tikh; lin_TTP_SPD_tissue lin_TTP_SPD; lin_TTP_TTV_tissue lin_TTP_TTV; lin_TTP_TIPS_TTV_tissue lin_TTP_TIPS_TTV;
    ];

lin = round(lin/0.01)*0.01;

t_all = toc(tid_all);
