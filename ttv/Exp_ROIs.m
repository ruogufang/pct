% Experiment: Clinical Perfusion Maps ROIs
% ROI in Figure 7 and 8
%   Ruogu Fang Revised 10/07/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

clear; close all; clc;

mA = 30;
rt = 1;

load(['data/IRB_AK_',num2str(mA),'mA_rt',num2str(rt),'.mat']);

y1 = 100; x1 = 10; w = 100; h = 100; % ROI for tube current
% y1 = 140; x1 = 170; w = 100; h = 100; % ROI for sampling rate
m = 6; % resize magnitude

CBF0_bSVD = CBF0_bSVD(y1:y1+h-1,x1:x1+w-1);
CBF0_bSVD = imresize(CBF0_bSVD,m,'nearest');

CBF_sSVD = CBF_sSVD(y1:y1+h-1,x1:x1+w-1);
CBF_sSVD = imresize(CBF_sSVD,m,'nearest');

CBF_bSVD = CBF_bSVD(y1:y1+h-1,x1:x1+w-1);
CBF_bSVD = imresize(CBF_bSVD,m,'nearest');

CBF_tikh = CBF_tikh(y1:y1+h-1,x1:x1+w-1);
CBF_tikh = imresize(CBF_tikh,m,'nearest');

CBF_SPD = CBF_SPD(y1:y1+h-1,x1:x1+w-1);
CBF_SPD = imresize(CBF_SPD,m,'nearest');

CBF_TTV = CBF_TTV(y1:y1+h-1,x1:x1+w-1);
CBF_TTV = imresize(CBF_TTV,m,'nearest');

pad = zeros(w*m,10);

CBF = [CBF0_bSVD pad CBF_sSVD pad CBF_bSVD pad CBF_tikh pad CBF_SPD pad CBF_TTV];
h = figure;
ctshow(CBF,[],[0 50]); colorbar;
