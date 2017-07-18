% pct_lcurve_demo
close all; clear; clc;

% Set path
addpath(genpath('../regu'));

% Parameters
PSNR = 20;
CBV = 4; %ml/100g
CBF = 20; % 20 : 10 : 80 ml/100g/min
X = 1; Y = 1; Z = 1; T = 60; sz = [T X Y]; % volumn size
m = 1;
regs = [0.001 0.01 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.25 0.5 0.75 1 1.25 1.5 2 5 10];
N = length(regs);

% Step 1: Generate Artery Input Function (AIF)
t = (0 : T-1)'; t0 = 5; a = 1; b = 3; c = 1.5;
AIF = pct_aif(t,t0,a,b,c);

% Step 2: Generate Impulse Residue Function RIF using Gamma Family
MTT = CBV / CBF * 60; %sec (3 ~ 12 s)
alpha = 10; beta = MTT / alpha;
RIF = pct_irf(t, alpha, beta); % Gamma RIF

% Step 3: Generate Tissue Time Enhancement Curve
TEC = pct_tec(AIF, RIF, CBF);

% Step 4: Correlated Gaussian Noise Generation
load acf;
sigma = max(TEC(:)) / (10^(PSNR/20));
% sigma = 0;
Cv = repmat(TEC,[1 X Y Z]);
Cn = pct_noise(Cv, acf, sigma,'s');
A = pct_circ(AIF,[],[],m); % block-circulant version
b = reshape(Cn,T,[]);

[px,py,regs] = pct_lcurve(A,b,sz,regs);

%% Make L-curve plot
plot(px,py,'-*','LineWidth',5), ax = axis;
xlabel('residual norm || CaK - C ||^2');
ylabel('solution norm || K ||_{TV}');
title('L-curve')

% Find the corner of L-curve
[reg_c,px_c,py_c] = l_corner(px,py,regs);