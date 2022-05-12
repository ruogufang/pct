%% gen_PMs_ori_demo.m
%
% Load CTP vol from .mat format
%
% Project: MMCTIQE
%
% Yao Xiao
% xiaoyao@ufl.edu
%  @ SMILE BME | UF
%
%% settings
close all; clear; clc;

% Set path
p_ctp = genpath('/Users/yaoxiao/Dropbox (UFL)/SMILE_Laptop_Doc/Code_Others/pct_201612/');
addpath(p_ctp);
p_pma = genpath('/Users/yaoxiao/Dropbox (UFL)/SMILE_Laptop_Doc/Manuals_Tutorials/PerfusionMaps/Display_PMA_Colormaps/');
addpath(p_pma);
color_lut = 'PMA_lut.csv';

path_data = '/Users/yaoxiao/Documents/SMILE_lab_files/CTP_TemporalSR/Data/';
path_data = strcat(path_data,'CTP_vol_Ori');

patient_id = '00000101'; %'00001801';
slice = 181;
preset_location = 'manual'; % select the method to obtain AIF, VOF locations: manual, default, ori, gt
show_fig = 1;

path_save = './perfusion_maps_demo_20200730/';
save_CBF = strcat(path_save,['PM_',patient_id,'_Ori_S',num2str(slice),'_CBF.pdf']);
save_CBV = strcat(path_save,['PM_',patient_id,'_Ori_S',num2str(slice),'_CBV.pdf']);
save_MTT = strcat(path_save,['PM_',patient_id,'_Ori_S',num2str(slice),'_MTT.pdf']);
save_perfusion_maps_mat = strcat(path_save,['PM_',patient_id,'_Ori_S',num2str(slice),'.mat']);
save_perfusion_maps_AIFVOF_Curves = strcat(path_save,['AIFVOF_',patient_id,'_Ori_S',num2str(slice),'_Curves.pdf']);
save_perfusion_maps_AIFVOF = strcat(path_save,['AIFVOF_',patient_id,'_Ori_S',num2str(slice),'.mat']);
path_perfusion_maps_AIFVOF = strcat(path_save,['AIFVOF_',patient_id,'_Ori_S',num2str(slice),'.mat']);

if ~exist(path_save,'dir'), mkdir(path_save); end

%% load CTP vol [YXTS] 512x512x21x10 and select slice for calculation
volume = load(fullfile(path_data,['CTP_vol_',patient_id,'_Ori_S',num2str(slice),'.mat']));

ctp_vol = volume.CTP_vol_ori;
%ctp_vol = squeeze(CTP_vol_ori(:,:,:,slice));

[Y, X, T] = size(ctp_vol);

if show_fig == 1
    figure; imshow(ctp_vol(:,:,11),[0 160]);
    title('CTP slice');
end

data = squeeze(ctp_vol);%V to data
data = permute(data,[3 1 2]); % to T Y X
% data = squeeze(V);
data = double(data);

% PCT Parameters
kappa = 0.73;  %Hematocrit correction factor
loth = 0;      %Lower segmentation threshold
hith = 120;    %Upper segmentation threshold
rho = 1.05;    %Average brain tissue density
fsize = 5;     %Size of Gaussian filter
PRE_bbbp = 1;  %First frame in BBBP calculation
POST_bbbp = T; %21; %Last frame in BBBP calculation
sigma = 20;    %Standard deviation of added Gaussian noise to CT data
dt = 0.5;      % Time step in time series
ftsize = 3;    %Size of temporal gaussian filter.
PRE = 1;
POST = T;

% SVD parameters
lambda = 0.15;    %Truncation parameter
m = 3; %2;        %Extend the data matrix m time for block circulant


% Compute Brain Mask
% B = squeeze(mean(data(1:10,:,:),1));%V to data
B = squeeze(mean(data(1:3,:,:),1));%V to data
% mask = pct_brainMask(B,0,120);
mask = pct_brainMask(B,0,120,15);

% Add Gaussian Noise
% data = pct_noise(data,[],sigma,'g');
% mA = pct_sigma2mA(sigma);

% Spatial filtering
data = pct_filter(data, fsize);

% Non-linear diffusive filtering
% data = pct_nldif(data);

% Time filtering
data = pct_gaussfilter(data,ftsize);

%Segmentation
data = pct_segment(data, loth, hith, PRE);

%Create time vector
% tv = [0:0.5:44 45:74];
% tv = [0:1.25:85 86.5:1.5:115]; % AK data

%Interpolation
% data = pct_interpolate(data,tv);

%Subtract base image
data = pct_contconv(data, PRE);

%Correct for hematocrit differences
data = pct_hematocrit(data, kappa);

% preset location
if strcmp(preset_location,'default')
    AIFx = 298;
    AIFy = 188;
    VOFx = 293;
    VOFy = 416;
end

% preset location
if strcmp(preset_location,'gt')
    load(fullfile(path_perfusion_maps_AIFVOF));
    AIFx = AIFx_gt; AIFy = AIFy_gt; 
    VOFx = VOFx_gt; VOFy = VOFy_gt;
end

% preset location
if strcmp(preset_location,'ori')
    load(fullfile(path_perfusion_maps_AIFVOF));
    AIFx = AIFx_ori; AIFy = AIFy_ori; 
    VOFx = VOFx_ori; VOFy = VOFy_ori;
end

% manually select location
if strcmp(preset_location,'manual')
    
    % Manual select AIF VOF
    [AIF_manual,AIFx,AIFy]=pct_aifsearch(data,5);
    [VOF_manual,VOFx,VOFy]=pct_aifsearch(data,5);
    
    AIFx_ori = AIFx; AIFy_ori = AIFy;
    VOFx_ori = VOFx; VOFy_ori = VOFy;
    save(fullfile(save_perfusion_maps_AIFVOF),'AIFx_ori','AIFy_ori','VOFx_ori','VOFy_ori');
end

% Auto select AIF VOF
if strcmp(preset_location,'auto')
	[AIF_auto,aif_x,aif_y]=pct_aifautoselect(data,mask);
end

%Get AIF and VOF
AIF = pct_tac(data, AIFy, AIFx);
VOF = pct_tac(data, VOFy, VOFx);

if show_fig == 1
    x_axis = 1:T;
    figure; plot(x_axis,AIF,'r',x_axis,VOF,'b');
    title('AIF VOF');
end

%Correct AIF for partial volume effects (PVE)
AIF = pct_aifscaling(AIF,VOF);

%Truncate the data
data = pct_truncate(data,PRE,POST);
AIF = AIF(PRE:POST);
VOF = VOF(PRE:POST);

if show_fig == 1
    x_axis = 1:T;
    h0=figure; plot(x_axis,AIF,'r',x_axis,VOF,'b');
    title('AIF VOF Correct');
    set(gcf,'color','white');
    set(gcf, 'InvertHardCopy', 'off');
    print(h0,fullfile(save_perfusion_maps_AIFVOF_Curves),'-dpdf');
end

%Calculate the residue functions
% R = cTSVD(AIF,data,lambda,m);
R = pct_bsvd(data,AIF,dt,lambda,m,mask);

%Calculate a CBV map
CBV = pct_cbv(R, rho);

%Calculate a CBF map
CBF = pct_cbf(R, rho);

%Multiply by a ratio, Why?
% CBF = CBF * 10;

%Calculate a MTT map
MTT = CBV ./ (CBF + eps) * 60;
% MTT = pct_mtt(R, rho);

% Threshold MTT
MTT(MTT<0)=0;
MTT(MTT>50)=0;


% PMA maps
% load the full PMA color lookup table
clt_pma = readtable(color_lut); 

% select colormaps from the full table
CLT_ASIST = select_colormap(clt_pma,'ASIST');

h1=figure;
ctshow_pma(CBF,mask,[],'pma',CLT_ASIST);
c = colorbar; c.Label.FontSize = 100; c.Color = 'white';
set(gcf,'color','black');
set(gcf, 'InvertHardCopy', 'off');
print(h1,fullfile(save_CBF),'-dpdf');

h2=figure;
ctshow_pma(CBV,mask,[],'pma',CLT_ASIST);
c = colorbar; c.Label.FontSize = 100; c.Color = 'white';
set(gcf,'color','black');
set(gcf, 'InvertHardCopy', 'off');
print(h2,fullfile(save_CBV),'-dpdf');

h3=figure;
ctshow_pma(MTT,mask,[],'pma',CLT_ASIST);
c = colorbar; c.Label.FontSize = 100; c.Color = 'white';
set(gcf,'color','black');
set(gcf, 'InvertHardCopy', 'off');
print(h3,fullfile(save_MTT),'-dpdf');


PMs = CBF;
PMs(:,:,2) = CBV;
PMs(:,:,3) = MTT;
PMs(:,:,4) = mask;

save(fullfile(save_perfusion_maps_mat),'PMs');

fprintf('--------- done --------\n');
