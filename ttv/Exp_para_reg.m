% Experiment: Tuning TV weighting parameter reg
% Figure 10(a)
% Ruogu Fang 6/19/2014
% Advanced Multimedia Laboratory
% Cornell University

close all; clear; clc;
warning('off');

% Set paths
addpath(genpath('toolbox')); % Include toolbox package

tid_all = tic;

%%%%%%% Tuning Parater %%%%%%%%%%
regs = 10.^(-6:6); % all TV weight to be evaluated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
m = 2; % block-circulant Ca extended for m times
lambda = 0.1; % cTSVD truncation parameter for singular values
rho = 1.05; %Average brain tissue density
rt = 1; % downsampling rate in time
Tn = 30; % truncated time in TTV

% Initialization
N = length(regs);
rmse_CBF_TTV = zeros(N,1);
lin_CBF_TTV = zeros(N,1);

% Load data
imgname = 'CTP.mat';
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
mA = 15;
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

clear V Vn VOF PRE POST B acf

% Method 0: block-circulant SVD (bSVD) no noise as reference
RIF0_bSVD = pct_bsvd(C,AIF0,1,lambda,m,Mask);
CBF0_bSVD = pct_cbf(RIF0_bSVD,rho);
CBV0_bSVD = pct_cbv(RIF0_bSVD,rho);
MTT0_bSVD = pct_mtt(RIF0_bSVD);

% Mask with vasuclar elimination
MaskV = pct_velim(CBV0_bSVD,Mask);
% CBF0_bSVD(~Mask)=0; CBV0_bSVD(~Mask) = 0; MTT0_bSVD(~Mask)=0;

%% Total Variation Regularization (TTV)
% parameters
[T,h,w]=size(Cn);
input.maxitr=500;
input.l=-inf; input.u=inf;
input.no = 5; % max number of iterations allowed
[Ca,Cc] = pct_circ(AIF,Cn,[],m); % block-circulant version
Tc = T*m;
input.A=Ca;
input.b=reshape(Cc,Tc,[]);
input.tt = Tn;
input.n1=h; input.n2=w; input.nt = Tc;

% clear C Ca

% TTV
for i = 1 : length(regs)
    input.reg = regs(i);
    out = TTV_FCSA(input);
    t_TTV = toc;
    RIF_TTV = reshape(out.y,[Tc,h,w]);
    RIF_TTV = RIF_TTV(1:T,:,:); % keep only first T time points
    CBF_TTV = pct_cbf(RIF_TTV,rho); CBF_TTV(~Mask) = 0;
    CBV_TTV = pct_cbv(RIF_TTV,rho); CBV_TTV(~Mask) = 0;
    MTT_TTV = pct_mtt(RIF_TTV); MTT_TTV(~Mask) = 0;


% RMSE and Lin's Concordance
rmse_CBF_TTV(i) = pct_rmse(CBF_TTV(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_TTV(i),ci_CBF_TTV] = pct_lincon(CBF_TTV(MaskV),CBF0_bSVD(MaskV));

fprintf('reg = %f: rmse = %.2f\n',regs(i),rmse_CBF_TTV(i));

end

%% Plot
N = length(regs);
figure;
[haxes,hline1,hline2] = plotyy(1:N,rmse_CBF_TTV,1:N,lin_CBF_TTV,'plot','plot');
xlabel('\gamma','FontSize',20);
ylabel('Performance','FontSize',20,'Color','b');
set(haxes(1),'FontSize',20,'XTick',1:3:N,'XTickLabel',regs(1:3:N));
set(haxes(2),'FontSize',20,'XTick',1:3:N,'XTickLabel',regs(1:3:N));
legend('RMSE','CCC');
set(hline1,'LineStyle','--','Marker','s','LineWidth',2,'Color','k','MarkerEdgeColor','b')
set(hline2,'LineStyle','--','Marker','o','LineWidth',2,'Color','k','MarkerEdgeColor','g')

t_all = toc;
