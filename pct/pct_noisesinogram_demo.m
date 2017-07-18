% pct_noisesinogram_demo
%Test function pct_noisesinogram
%
% Ruogu Fang
% 4/15/2013 

close all; clear; clc;

% mu_water=0.01835;% attenuation factor of water
% mu_water=0.19;

% User-defined phantom
E = [1200 0.69 0.92 0 0 0;
    -480 0.6224 0.874 0 -0.0184 0;
    -120 0.41 0.16 -0.22 0 108;
    -120 0.31 0.11 0.22 0 72;
    60 0.25 0.21 0 0.35 90;
    60 0.046 0.046 0 0.1 0;
    60 0.023 0.046 0 -0.1 0;
    80 0.046 0.023 -0.08 -0.605 0;
    80 0.023 0.023 0 -0.605 0;
    80 0.046 0.023 0.06 -0.605 90;
    36 0.029 0.029 0 0.61 0];

C = phantom(512);

% %1/11 of low level
% I0 = 2.5e5;
% sigma = sqrt(10);

I0= 1e6; %noise free
I0=2*10^5; % low-dose
sigma = sqrt(10);

Cn = pct_noisesinogram(C, I0, sigma);

figure;imshow(C,[0 1]);
figure;imshow(Cn,[0 1]);
