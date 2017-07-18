% TTV-CTP CT perfusion deconvolution using total variation regularization on
% the residue functions with non-iso regularization parameter lambda
% 
% Figure 7 & 8
%
% Ruogu Fang 9/21/2013 Advanced Multimedia Laboratory

close all; 
% clear; clc;
warning('off');

% Set paths
addpath(genpath('toolbox')); % Include toolbox package

%%%%%%%%%% Setting parameters %%%%%%%%%%%%%%%%%%
mA = 15; % tube current-exposure time product
rt = 1; % downsampling rate in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tid_all = tic;

% Parameters
m = 2; % block-circulant Ca extended for m times
lambda = 0.1; % cTSVD truncation parameter for singular values
rho = 1.05; %Average brain tissue density
Tn = 30; % truncated time in TTV
% dt = 0.25; % sampling frequency for AK
dt = 0.5; % sampling frequency

% Load data
% imgname = 'CTP';
imgname = 'IRB_AK';
load(imgname);
% Data format:
%   V: CTP data [T x X x Y] AIFx, AIFy, VOFx, VOFy: aif and vof
%   coordinations PRE: Pre-enhancement cutoff (first frame included in
%   calc's) POST: Post-enhancement cutoff (last frame included in calc's)

load acf;

% Region of Interest
% x0 = 200; y0 = 200; w = 80; h = 80; % RACA
% x0 = 150; y0 = 150; w = 200; h = 200; % RMCA
% x0 = 1; y0 = 1; w = 512; h = 512; % Whole image

% noisy level
mA0 = 190;
sigma = pct_mA2sigma(mA,mA0);

% Remove negative values
V(V<0) = 0;

% Find Brain Mask
B = squeeze(mean(V(1:10,:,:),1));
Mask = pct_brainMask(B,0,120);

% Add correlated Gaussian (spectral) noise to simulate low-dose
Vn = pct_noise(permute(V,[2 3 1]),acf,sigma,'s',Mask);
Vn = permute(Vn,[3 1 2]);

% PCT preprocess
% Whole Brain CBF
[C,AIF0] = pct_preprocess(V, PRE, POST, AIFx,AIFy,VOFx,VOFy); % noiseless data
[Cn,AIF] = pct_preprocess(Vn, PRE, POST, AIFx,AIFy,VOFx,VOFy); % noisy data

% Crop data to Region of Interest
if exist('x0','var')
    C = C(:,y0:y0+h-1,x0:x0+w-1);
    Cn = Cn(:,y0:y0+h-1,x0:x0+w-1);
    Mask = Mask(y0:y0+h-1,x0:x0+w-1);
end

% Find the minimum bounding box
bb = pct_minBoundingBox(Mask);
C = C(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
Cn = Cn(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
Mask = Mask(bb(1,1):bb(1,2),bb(2,1):bb(2,3));

% Downsampling
Cn = Cn(1:rt:end,:,:);
AIF = AIF(1:rt:end);
Tn = round(Tn/rt);

clear V Vn VOF PRE POST B acf

% Method 0: block-circulant SVD (bSVD) no noise as reference
RIF0_bSVD = pct_bsvd(C,AIF0,dt,lambda,m,Mask);
CBF0_bSVD = pct_cbf(RIF0_bSVD,rho);
CBV0_bSVD = pct_cbv(RIF0_bSVD,rho);
MTT0_bSVD = pct_mtt(RIF0_bSVD);
TTP0_bSVD = pct_ttp(pct_tec(AIF0,RIF0_bSVD),dt,Mask);

% Mask with vasuclar elimination
MaskV = pct_velim(CBV0_bSVD,Mask);
% CBF0_bSVD(~Mask)=0; CBV0_bSVD(~Mask) = 0; MTT0_bSVD(~Mask)=0;

% Method 1: Standard Truncated Singular Value Decomposition (sSVD)
tic;
RIF_sSVD = pct_ssvd(Cn,AIF,dt*rt,lambda,Mask);
CBF_sSVD = pct_cbf(RIF_sSVD,rho);
CBV_sSVD = pct_cbv(RIF_sSVD,rho);
MTT_sSVD = pct_mtt(RIF_sSVD);
TTP_sSVD = pct_ttp(pct_tec(AIF,RIF_sSVD),dt*rt,Mask);
t_sSVD = toc;

% Method 2: block-circulant SVD (bSVD)
tic;
RIF_bSVD = pct_bsvd(Cn,AIF,dt*rt,lambda,m,Mask);
CBF_bSVD = pct_cbf(RIF_bSVD,rho);
CBV_bSVD = pct_cbv(RIF_bSVD,rho);
MTT_bSVD = pct_mtt(RIF_bSVD);
TTP_bSVD = pct_ttp(pct_tec(AIF,RIF_bSVD),dt*rt,Mask);
t_bSVD = toc;

% Method 3: Tikhonov
tic;
RIF_tikh = pct_tikh(Cn,AIF,dt*rt,lambda,1,Mask);
CBF_tikh = pct_cbf(RIF_tikh,rho);
CBV_tikh = pct_cbv(RIF_tikh,rho);
MTT_tikh = pct_mtt(RIF_tikh);
TTP_tikh = pct_ttp(pct_tec(AIF,RIF_tikh),dt*rt,Mask);
t_tikh = toc;

% Method 4: TIPS+bSVD
tic;

% Set bilateral filter parameters.
w_tips     = 5;       % bilateral filter half-width
sigma_tips = [3 0.1]; % bilateral filter standard deviations

tid_tips = tic;
Cn_TIPS = pct_tips(Cn,w_tips,sigma_tips);
t_tips = toc(tid_tips);
RIF_TIPS_bSVD = pct_bsvd(Cn_TIPS,AIF,dt*rt,lambda,m,Mask);
RIF_TIPS_bSVD(RIF_TIPS_bSVD<0) = 0;
CBF_TIPS_bSVD = pct_cbf(RIF_TIPS_bSVD,rho);
CBV_TIPS_bSVD = pct_cbv(RIF_TIPS_bSVD,rho);
MTT_TIPS_bSVD = pct_mtt(RIF_TIPS_bSVD);
TTP_TIPS_bSVD = pct_ttp(pct_tec(AIF,RIF_TIPS_bSVD),dt*rt,Mask);
t_TIPS_bSVD = toc;

%% Method 4: Global Sparse Perfusion Deconvolution (GSPD)

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
params.lambda = 0.1; % weight of the sparsity term in SPAMS package
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
params.rt = dt*rt;

% %Load pre-trained dictionary
load D_ODL_6x6.mat
params.D = dict;

% Online deconvolution
tic;
[CBF_SPD, MTT_SPD, CBV_SPD, TTP_SPD] = spd(params);
t_spd = toc;

%% Method 5: Total Variation Regularization (TTV)
% parameters
[T,h,w]=size(Cn);
% input.reg = 1e-4; % when rt=1
input.reg = [1e-4 1] * 1e-4;
input.maxitr=500;
input.l=-inf; input.u=inf;
input.no = 20; % max number of iterations allowed
[Ca,Cc] = pct_circ(AIF,Cn,[],m); % block-circulant version
Tc = T*m;
input.A=Ca;
input.b=reshape(Cc,Tc,[]);
input.tt = Tn;
input.n1=h; input.n2=w; input.nt = Tc;
input.dt = dt*rt;

% TTV
tic;
out = TTV_FCSA(input);
t_TTV = toc;
RIF_TTV = reshape(out.y,[Tc,h,w]);
RIF_TTV = RIF_TTV(1:T,:,:); % keep only first T time points
CBF_TTV = pct_cbf(RIF_TTV,rho); CBF_TTV(~Mask) = 0;
CBV_TTV = pct_cbv(RIF_TTV,rho); CBV_TTV(~Mask) = 0;
MTT_TTV = pct_mtt(RIF_TTV); MTT_TTV(~Mask) = 0;
TTP_TTV = pct_ttp(pct_tec(AIF,RIF_TTV),dt*rt,Mask);
fprintf(1,'TTV: Iter_time = %.2fsec, Iteration Num = %d\n',out.xtime(end), input.no);

% Method 7: TIPS+TTV
tid_tips_ttv = tic;
[Ca,Cc_TIPS] = pct_circ(AIF,Cn_TIPS,[],m); % block-circulant version
input.b = reshape(Cc_TIPS,Tc,[]);
out = TTV_FCSA(input);
RIF_TIPS_TTV = reshape(out.y,[Tc,h,w]);
% [RIF_TIPS_TTV,funv] = pct_ttv(input);
RIF_TIPS_TTV = RIF_TIPS_TTV(1:T,:,:); % keep only first T time points
CBF_TIPS_TTV = pct_cbf(RIF_TIPS_TTV,rho,Mask);
CBV_TIPS_TTV = pct_cbv(RIF_TIPS_TTV,rho,Mask);
MTT_TIPS_TTV = pct_mtt(RIF_TIPS_TTV,Mask);
TTP_TIPS_TTV = pct_ttp(pct_tec(AIF,RIF_TIPS_TTV),dt*rt,Mask);
t_tips_ttv = toc(tid_tips_ttv);

%% Plot results
[X,Y] = size(CBF0_bSVD);
pad = zeros(X,10);
MaskV_all = repmat([MaskV pad],[1 8]); MaskV_all=MaskV_all(:,1:end-10);
Mask_all = repmat([Mask pad],[1 8]); Mask_all=Mask_all(:,1:end-10);

CBF = [CBF0_bSVD pad CBF_sSVD pad CBF_bSVD pad CBF_TIPS_bSVD pad CBF_tikh pad CBF_SPD pad CBF_TTV pad CBF_TIPS_TTV];
% CBF = [CBF0_bSVD pad CBF_sSVD pad CBF_bSVD pad CBF_TIPS_bSVD pad CBF_tikh pad CBF_TTV pad CBF_TIPS_TTV];

% CBV = [CBV0_bSVD CBV_sSVD CBV_bSVD CBV_TIPS_bSVD CBV_tikh CBV_SPD CBV_TTV CBV_TIPS_TTV];
% MTT = [MTT0_bSVD MTT_sSVD MTT_bSVD MTT_TIPS_bSVD MTT_tikh MTT_SPD MTT_TTV MTT_TIPS_TTV];
% TTP = [TTP0_bSVD TTP_sSVD TTP_bSVD TTP_TIPS_bSVD TTP_tikh TTP_SPD TTP_TTV TTP_TIPS_TTV];

% y1 = 60; x1 = 60; w = 200; h = 160; % ROI for tube current 018
y1 = 80; x1 = 50; w = 180; h = 140; % ROI for tube current 022 & AK
% y1 = 140; x1 = 170; w = 100; h = 100; % ROI for sampling rate

h1=figure; ctshow(CBF,MaskV_all,[0 80]);
colorbar; set(gca,'FontSize',20);
hold on;
rectangle('Position',[x1,y1,w,h],'EdgeColor',[1 1 1],'LineWidth',2);

% Save figures
figpath = '../../../Figures/';
% W = 10; H = 2;
% set(h1, 'papersize', [W H]);
% set(h1, 'paperposition', [0 0 W H]);
% print(h1,fullfile(figpath,'pdf',[imgname '_',num2str(mA) 'mA_rt' num2str(rt) '.pdf']),'-dpdf');
% saveas(h1,fullfile(figpath,'fig',[imgname '_',num2str(mA) 'mA_rt' num2str(rt) '.fig']));

% %% ROI
% 
% s = 6; % resize magnitude
% 
% CBF0_bSVD_ROI = CBF0_bSVD(y1:y1+h-1,x1:x1+w-1);
% CBF0_bSVD_ROI = imresize(CBF0_bSVD_ROI,s,'nearest');
% 
% CBF_sSVD_ROI = CBF_sSVD(y1:y1+h-1,x1:x1+w-1);
% CBF_sSVD_ROI = imresize(CBF_sSVD_ROI,s,'nearest');
% 
% CBF_bSVD_ROI = CBF_bSVD(y1:y1+h-1,x1:x1+w-1);
% CBF_bSVD_ROI = imresize(CBF_bSVD_ROI,s,'nearest');
% 
% CBF_TIPS_bSVD_ROI = CBF_bSVD(y1:y1+h-1,x1:x1+w-1);
% CBF_TIPS_bSVD_ROI = imresize(CBF_TIPS_bSVD_ROI,s,'nearest');
% 
% CBF_tikh_ROI = CBF_tikh(y1:y1+h-1,x1:x1+w-1);
% CBF_tikh_ROI = imresize(CBF_tikh_ROI,s,'nearest');
% 
% CBF_SPD_ROI = CBF_SPD(y1:y1+h-1,x1:x1+w-1);
% CBF_SPD_ROI = imresize(CBF_SPD_ROI,s,'nearest');
% 
% CBF_TTV_ROI = CBF_TTV(y1:y1+h-1,x1:x1+w-1);
% CBF_TTV_ROI = imresize(CBF_TTV_ROI,s,'nearest');
% 
% CBF_TIPS_TTV_ROI = CBF_TTV(y1:y1+h-1,x1:x1+w-1);
% CBF_TIPS_TTV_ROI = imresize(CBF_TIPS_TTV_ROI,s,'nearest');
% 
% pad = zeros(h*s,10);
% 
% MaskV_ROI = MaskV(y1:y1+h-1,x1:x1+w-1);
% MaskV_ROI = imresize(MaskV_ROI,s,'nearest');
% MaskV_ROI_all = repmat([MaskV_ROI pad],[1 8]); MaskV_ROI_all=MaskV_ROI_all(:,1:end-10);
% 
% CBF_ROI = [CBF0_bSVD_ROI pad CBF_sSVD_ROI pad CBF_bSVD_ROI pad CBF_TIPS_bSVD_ROI pad CBF_tikh_ROI pad CBF_SPD_ROI pad CBF_TTV_ROI pad CBF_TIPS_TTV_ROI];
% 
% h2=figure; ctshow(CBF_ROI,MaskV_ROI_all,[0 80]);
% % colorbar; set(gca,'FontSize',20);
% 
% w = 25; h = 5;
% set(h2, 'papersize', [w h]);
% set(h2, 'paperposition', [0 0 w h]);
% print(h2,fullfile(figpath,'pdf',[imgname '_ROI_',num2str(mA) 'mA.png']),'-dpng');
% saveas(h2,fullfile(figpath,'fig',[imgname '_ROI_',num2str(mA) 'mA.fig']));

% figure;
% set(gca,'FontSize',20);
% subplot(311); ctshow(CBF,[],[0 80]); colorbar;
% subplot(312); ctshow(CBV,[],[0 20]); colorbar;
% subplot(313); ctshow(MTT,[],[0 15]); colorbar;
% subplot(313); ctshow(TTP,[],[0 50]); colorbar;
% 
figure; plot(out.funv); title('TTV funv');

%% RMSE, PSNR, Lin's Concordance
rmse_CBF_sSVD = pct_rmse(CBF_sSVD,CBF0_bSVD,MaskV);
rmse_CBF_bSVD = pct_rmse(CBF_bSVD,CBF0_bSVD,MaskV);
rmse_CBF_tikh = pct_rmse(CBF_tikh,CBF0_bSVD,MaskV);
rmse_CBF_TIPS_bSVD = pct_rmse(CBF_TIPS_bSVD,CBF0_bSVD,MaskV);
rmse_CBF_SPD = pct_rmse(CBF_SPD,CBF0_bSVD,MaskV);
rmse_CBF_TTV = pct_rmse(CBF_TTV,CBF0_bSVD,MaskV);
rmse_CBF_TIPS_TTV = pct_rmse(CBF_TIPS_TTV,CBF0_bSVD,MaskV);

rmse_CBV_sSVD = pct_rmse(CBV_sSVD,CBV0_bSVD,MaskV);
rmse_CBV_bSVD = pct_rmse(CBV_bSVD,CBV0_bSVD,MaskV);
rmse_CBV_tikh = pct_rmse(CBV_tikh,CBV0_bSVD,MaskV);
rmse_CBV_TIPS_bSVD = pct_rmse(CBV_TIPS_bSVD,CBV0_bSVD,MaskV);
rmse_CBV_SPD = pct_rmse(CBV_SPD,CBV0_bSVD,MaskV);
rmse_CBV_TTV = pct_rmse(CBV_TTV,CBV0_bSVD,MaskV);
rmse_CBV_TIPS_TTV = pct_rmse(CBV_TIPS_TTV,CBV0_bSVD,MaskV);

rmse_MTT_sSVD = pct_rmse(MTT_sSVD,MTT0_bSVD,MaskV);
rmse_MTT_bSVD = pct_rmse(MTT_bSVD,MTT0_bSVD,MaskV);
rmse_MTT_tikh = pct_rmse(MTT_tikh,MTT0_bSVD,MaskV);
rmse_MTT_TIPS_bSVD = pct_rmse(MTT_TIPS_bSVD,MTT0_bSVD,MaskV);
rmse_MTT_SPD = pct_rmse(MTT_SPD,MTT0_bSVD,MaskV);
rmse_MTT_TTV = pct_rmse(MTT_TTV,MTT0_bSVD,MaskV);
rmse_MTT_TIPS_TTV = pct_rmse(MTT_TIPS_TTV,MTT0_bSVD,MaskV);

rmse_TTP_sSVD = pct_rmse(TTP_sSVD,TTP0_bSVD,MaskV);
rmse_TTP_bSVD = pct_rmse(TTP_bSVD,TTP0_bSVD,MaskV);
rmse_TTP_tikh = pct_rmse(TTP_tikh,TTP0_bSVD,MaskV);
rmse_TTP_TIPS_bSVD = pct_rmse(TTP_TIPS_bSVD,TTP0_bSVD,MaskV);
rmse_TTP_SPD = pct_rmse(TTP_SPD,TTP0_bSVD,MaskV);
rmse_TTP_TTV = pct_rmse(TTP_TTV,TTP0_bSVD,MaskV);
rmse_TTP_TIPS_TTV = pct_rmse(TTP_TIPS_TTV,TTP0_bSVD,MaskV);

psnr_CBF_sSVD = pct_psnr(CBF_sSVD,CBF0_bSVD,MaskV);
psnr_CBF_bSVD = pct_psnr(CBF_bSVD,CBF0_bSVD,MaskV);
psnr_CBF_tikh = pct_psnr(CBF_tikh,CBF0_bSVD,MaskV);
psnr_CBF_TIPS_bSVD = pct_psnr(CBF_TIPS_bSVD,CBF0_bSVD,MaskV);
psnr_CBF_SPD = pct_psnr(CBF_SPD,CBF0_bSVD,MaskV);
psnr_CBF_TTV = pct_psnr(CBF_TTV,CBF0_bSVD,MaskV);
psnr_CBF_TIPS_TTV = pct_psnr(CBF_TIPS_TTV,CBF0_bSVD,MaskV);

psnr_CBV_sSVD = pct_psnr(CBV_sSVD,CBV0_bSVD,MaskV);
psnr_CBV_bSVD = pct_psnr(CBV_bSVD,CBV0_bSVD,MaskV);
psnr_CBV_tikh = pct_psnr(CBV_tikh,CBV0_bSVD,MaskV);
psnr_CBV_TIPS_bSVD = pct_psnr(CBV_TIPS_bSVD,CBV0_bSVD,MaskV);
psnr_CBV_SPD = pct_psnr(CBV_SPD,CBV0_bSVD,MaskV);
psnr_CBV_TTV = pct_psnr(CBV_TTV,CBV0_bSVD,MaskV);
psnr_CBV_TIPS_TTV = pct_psnr(CBV_TIPS_TTV,CBV0_bSVD,MaskV);

psnr_MTT_sSVD = pct_psnr(MTT_sSVD,MTT0_bSVD,MaskV);
psnr_MTT_bSVD = pct_psnr(MTT_bSVD,MTT0_bSVD,MaskV);
psnr_MTT_tikh = pct_psnr(MTT_tikh,MTT0_bSVD,MaskV);
psnr_MTT_TIPS_bSVD = pct_psnr(MTT_TIPS_bSVD,MTT0_bSVD,MaskV);
psnr_MTT_SPD = pct_psnr(MTT_SPD,MTT0_bSVD,MaskV);
psnr_MTT_TTV = pct_psnr(MTT_TTV,MTT0_bSVD,MaskV);
psnr_MTT_TIPS_TTV = pct_psnr(MTT_TIPS_TTV,MTT0_bSVD,MaskV);

psnr_TTP_sSVD = pct_psnr(TTP_sSVD,TTP0_bSVD,MaskV);
psnr_TTP_bSVD = pct_psnr(TTP_bSVD,TTP0_bSVD,MaskV);
psnr_TTP_tikh = pct_psnr(TTP_tikh,TTP0_bSVD,MaskV);
psnr_TTP_TIPS_bSVD = pct_psnr(TTP_TIPS_bSVD,TTP0_bSVD,MaskV);
psnr_TTP_SPD = pct_psnr(TTP_SPD,TTP0_bSVD,MaskV);
psnr_TTP_TTV = pct_psnr(TTP_TTV,TTP0_bSVD,MaskV);
psnr_TTP_TIPS_TTV = pct_psnr(TTP_TIPS_TTV,TTP0_bSVD,MaskV);

[lin_CBF_sSVD,ci_CBF_sSVD] = pct_lincon(CBF_sSVD(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_bSVD,ci_CBF_bSVD] = pct_lincon(CBF_bSVD(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_tikh,ci_CBF_tikh] = pct_lincon(CBF_tikh(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_TIPS_bSVD,ci_CBF_TIPS_bSVD] = pct_lincon(CBF_TIPS_bSVD(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_SPD,ci_CBF_SPD] = pct_lincon(CBF_SPD(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_TTV,ci_CBF_TTV] = pct_lincon(CBF_TTV(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_TIPS_TTV,ci_CBF_TIPS_TTV] = pct_lincon(CBF_TIPS_TTV(MaskV),CBF0_bSVD(MaskV));

[lin_CBV_sSVD,ci_CBV_sSVD] = pct_lincon(CBV_sSVD(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_bSVD,ci_CBV_bSVD] = pct_lincon(CBV_bSVD(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_tikh,ci_CBV_tikh] = pct_lincon(CBV_tikh(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_TIPS_bSVD,ci_CBV_TIPS_bSVD] = pct_lincon(CBV_TIPS_bSVD(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_SPD,ci_CBV_SPD] = pct_lincon(CBV_SPD(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_TTV,ci_CBV_TTV] = pct_lincon(CBV_TTV(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_TIPS_TTV,ci_CBV_TIPS_TTV] = pct_lincon(CBV_TIPS_TTV(MaskV),CBV0_bSVD(MaskV));

[lin_MTT_sSVD,ci_MTT_sSVD] = pct_lincon(MTT_sSVD(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_bSVD,ci_MTT_bSVD] = pct_lincon(MTT_bSVD(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_tikh,ci_MTT_tikh] = pct_lincon(MTT_tikh(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_TIPS_bSVD,ci_MTT_TIPS_bSVD] = pct_lincon(MTT_TIPS_bSVD(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_SPD,ci_MTT_SPD] = pct_lincon(MTT_SPD(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_TTV,ci_MTT_TTV] = pct_lincon(MTT_TTV(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_TIPS_TTV,ci_MTT_TIPS_TTV] = pct_lincon(MTT_TIPS_TTV(MaskV),MTT0_bSVD(MaskV));

[lin_TTP_sSVD,ci_TTP_sSVD] = pct_lincon(TTP_sSVD(MaskV),TTP0_bSVD(MaskV));
[lin_TTP_bSVD,ci_TTP_bSVD] = pct_lincon(TTP_bSVD(MaskV),TTP0_bSVD(MaskV));
[lin_TTP_tikh,ci_TTP_tikh] = pct_lincon(TTP_tikh(MaskV),TTP0_bSVD(MaskV));
[lin_TTP_TIPS_bSVD,ci_TTP_TIPS_bSVD] = pct_lincon(TTP_TIPS_bSVD(MaskV),TTP0_bSVD(MaskV));
[lin_TTP_SPD,ci_TTP_SPD] = pct_lincon(TTP_SPD(MaskV),TTP0_bSVD(MaskV));
[lin_TTP_TTV,ci_TTP_TTV] = pct_lincon(TTP_TTV(MaskV),TTP0_bSVD(MaskV));
[lin_TTP_TIPS_TTV,ci_TTP_TIPS_TTV] = pct_lincon(TTP_TIPS_TTV(MaskV),TTP0_bSVD(MaskV));

P_sSVD = polyfit(CBF_sSVD,CBF0_bSVD,1);
P_bSVD = polyfit(CBF_bSVD,CBF0_bSVD,1);
P_tikh = polyfit(CBF_tikh,CBF0_bSVD,1);
P_TIPS_bSVD = polyfit(CBF_TIPS_bSVD,CBF0_bSVD,1);
P_SPD = polyfit(CBF_SPD,CBF0_bSVD,1);
P_TTV = polyfit(CBF_TTV,CBF0_bSVD,1);
P_TIPS_TTV = polyfit(CBF_TIPS_TTV,CBF0_bSVD,1);

t_all = toc(tid_all);

% Save Results
% save(['data/',imgname,'_',num2str(mA),'mA_rt',num2str(rt),'.mat'],'imgname','CBF*','CBV*','MTT*','TTP*','rmse*','psnr*','lin*','ci*','P*','Mask*','t*','mA','rt');
