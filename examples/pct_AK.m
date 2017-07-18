%Kolbeinn Karlsson 06/11/12
%Advanced Multimedia Processing (AMP) Lab, Cornell University

addpath('tprod');
addpath('~/Dropbox/Research/DICOM/CT_MAT');

% PCT Parameters
PRE = 19;     %Pre-enhancement cutoff (first frame included in calc's)
POST = 89;   %Post-enhancement cutoff (last frame included in calc's)
kappa = 0.73; %Hematocrit correction factor
loth = 0;     %Lower segmentation threshold
hith = 120;   %Upper segmentation threshold
rho = 1.05;   %Average brain tissue density
fsize = 5;    %Size of Gaussian filter
tsize = 5;    %Size of temporal average filter
slice = 1;      % Slice number

% SVD parameters
lambda = 0.3; %Truncation parameter
m = 2;        %Extend the data matrix m time for block circulant

% Reference artery and veins
AIFx = 244;
AIFy = 154;
VOFx = 234;
VOFy = 434;

% Get the data
load 'CTP_AK.mat'
V = data; % 4 slices
data = squeeze(V(:,:,slice,:));
data = double(data);

%Low-dose simulation
sigma = 0;
mA = 190*103^2/(190*sigma^2+103^2);
data = pct_noise(data, [], sigma);

% Make the indexing right
data = permute(data, [3 1 2]);

% Spatial filtering
data = pct_filter(data, fsize);

% Temporal filtering
data = pct_tfilter(data, tsize);

%Segmentation
data = pct_segment(data, loth, hith, PRE);

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
MTT = CBV ./ (CBF + eps) * 60;

%Thresho values
MTT(MTT<0)=0;
MTT(MTT>50)=0;

%Multiply by a ratio of 2. Why?
MTT = MTT * 2;

%Remove non-brain tissue
mask = pct_brainMask(CBF);
CBF = pct_mask(CBF,mask);
CBV = pct_mask(CBV,mask);
MTT = pct_mask(MTT,mask);

% % Scale for saving
% CBV = CBV * 10;
% CBF = CBF * 4;
% MTT = MTT * 100;

% %Write to DICOM
% dicomwrite(int16(CBF),'OUTPUT/CBF_TSVD_11mA_002.dcm');
% dicomwrite(int16(CBV),'OUTPUT/CBV_TSVD_11mA_002.dcm');
% dicomwrite(int16(MTT),'OUTPUT/MTT_TSVD_11mA_002.dcm');