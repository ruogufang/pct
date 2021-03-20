%% example_ctshow_pma_colormap_full.m
%
% Display perfusion maps by using PMA colormaps
% 1) Load the full PMA color lookup table
% 2) Select the Colormap that want to use
% 3) Load precalcualted perfusion maps, CBF, CBV
% 4) Display images with different color maps
%
% Yao Xiao
% Feb 14, 2019 @ SMILE | UF
%
%% settings
close all; clear; clc;

addpath(genpath(cd));

path_data = '.';
path_lut = './data';
path_save = './results_me';

if ~exist(path_save,'dir'), mkdir(path_save); end

%% load the full PMA color lookup table
clt_pma = readtable(fullfile(path_lut,'PMA_lut.csv'));

%% select colormaps from the full table
CLT_ASIST = select_colormap(clt_pma,'ASIST');
CLT_Siemens_CT = select_colormap(clt_pma,'Siemens_CT');
CLT_Philips_CBF = select_colormap(clt_pma,'Philips_CBF');
CLT_Philips_CBV = select_colormap(clt_pma,'Philips_CBV');
CLT_GE_3_Colors = select_colormap(clt_pma,'GE_3_Colors');

%% load precalculated perfusion maps
load(fullfile(path_data,'rCBF_test4_image.mat')); % precalculated CBF
load(fullfile(path_data,'rCBV_test4_image.mat')); % precalculated CBV
load(fullfile(path_data,'mask.mat')); % precalculated mask

CBF = rCBF_test4_image; CBV = rCBV_test4_image;
%% display images
figure; 

% default colormap
p1 = subplot(5,2,1);
[~,cm1,~] = ctshow_pma(CBF,mask,[],'default');
title('CBF default');

p2 = subplot(5,2,2);
[~,cm2,~] = ctshow_pma(CBV,mask,[],'default');
title('CBV default');

% PMA colormap ASIST
p3 = subplot(5,2,3);
[~,cm3,~] = ctshow_pma(CBF,mask,[],'pma',CLT_ASIST);
title('CBF PMA colormap ASIST');

p4 = subplot(5,2,4);
[~,cm4,~] = ctshow_pma(CBV,mask,[],'pma',CLT_ASIST);
title('CBV PMA colormap ASIST');

% PMA colormap Siemens_CT
p5 = subplot(5,2,5);
[~,cm5,~] = ctshow_pma(CBF,mask,[],'pma',CLT_Siemens_CT);
title('CBF PMA colormap Siemens CT');

p6 = subplot(5,2,6);
[~,cm6,~] = ctshow_pma(CBV,mask,[],'pma',CLT_Siemens_CT);
title('CBV PMA colormap Siemens CT');

% PMA colormap Philips_CBF & Philips_CBV
p7 = subplot(5,2,7);
[~,cm7,~] = ctshow_pma(CBF,mask,[],'pma',CLT_Philips_CBF);
title('CBF PMA colormap Philips CBF');

p8 = subplot(5,2,8);
[~,cm8,~] = ctshow_pma(CBV,mask,[],'pma',CLT_Philips_CBV);
title('CBV PMA colormap Philips CBV');

% PMA colormap PGE_3_Colors
p9 = subplot(5,2,9);
[~,cm9,~] = ctshow_pma(CBF,mask,[],'pma',CLT_GE_3_Colors);
title('CBF PMA colormap GE 3 Colors');

p10 = subplot(5,2,10);
[~,cm10,~] = ctshow_pma(CBV,mask,[],'pma',CLT_GE_3_Colors);
title('CBV PMA colormap GE 3 Colors');


% display color for different maps
colormap(p1,cm1); % default
colormap(p2,cm2); % default
colormap(p3,cm3); % ASIST
colormap(p4,cm4); % ASIST
colormap(p5,cm5); % Siemens_CT
colormap(p6,cm6); % Siemens_CT
colormap(p7,cm7); % Philips_CBF
colormap(p8,cm8); % Philips_CBV
colormap(p9,cm9); % PGE_3_Colors
colormap(p10,cm10); % PGE_3_Colors

% save image to file
saveas(gcf,fullfile(path_save,'result.jpg'));
saveas(gcf,fullfile(path_save,'result.fig'));

