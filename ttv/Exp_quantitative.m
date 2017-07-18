% Quantitative Evaluation on Clinical Data in Table form
% Ruogu Fang
% 12/17/2014

close all; clear; clc;
addpath('data');
filenames = dir(fullfile('data','*_15mA_rt1.mat'));

PSNR = []; RMSE = []; Lin = []; P = {};

for i = 1 : length(filenames)
    i;
    filename = filenames(i).name;
    load(filename);
    %     imgname = filename(1:end-13);
    %     rt = 1;
    %     quantify;
    rmse = [rmse_CBF_sSVD; rmse_CBF_bSVD; rmse_CBF_tikh; rmse_CBF_TIPS_bSVD; rmse_CBF_SPD; rmse_CBF_TTV; rmse_CBF_TIPS_TTV];
    RMSE = [RMSE rmse];
    psnr = [psnr_CBF_sSVD; psnr_CBF_bSVD; psnr_CBF_tikh; psnr_CBF_TIPS_bSVD; psnr_CBF_SPD; psnr_CBF_TTV; psnr_CBF_TIPS_TTV];
    PSNR = [PSNR psnr];
    lin = [lin_CBF_sSVD; lin_CBF_bSVD; lin_CBF_tikh; lin_CBF_TIPS_bSVD; lin_CBF_SPD; lin_CBF_TTV; lin_CBF_TIPS_TTV];
    Lin = [Lin lin];
    p = {sprintf('y=%.3fx+%.3f',P_sSVD(1),P_sSVD(2));sprintf('y=%.3fx+%.3f',P_bSVD(1),P_bSVD(2));sprintf('y=%.3fx+%.3f',P_tikh(1),P_tikh(2));sprintf('y=%.3fx+%.3f',P_TIPS_bSVD(1),P_TIPS_bSVD(2));sprintf('y=%.3fx+%.3f',P_SPD(1),P_SPD(2));sprintf('y=%.3fx+%.3f',P_TTV(1),P_TTV(2));sprintf('y=%.3fx+%.3f',P_TIPS_TTV(1),P_TIPS_TTV(2))};
    P = [P p];
    save('data/results.mat','PSNR','RMSE','Lin','P');
end


