%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Cn = pct_noisesinogram(C,I0,sigma)
%
% Simulate low-dose recontructed CTP data by adding Poisson and Gaussian
% noise to the sinogram generated from fanbeam
%
% INPUT:
%   C       - CTP time series data [X x Y x T]
%   I0      - Incident X-ray intensity, or incident photon number (10^6 for normal-dose, 2.5*10^5 for 1/7 of normal dose)
%   sigma   - Standard deviation of Gaussian noise in low-dose (0 for normal-dose, 10 for 1/7 of normal dose)
%
% OUTPUT:
%   Cn      - The noisy CTP volume [X x Y x T]
%
% Ruogu Fang
% 4/15/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cn = pct_noisesinogram(C, I0, sigma)

C = double(C);
Mask = mean(C(:,:,1:min(size(C,3),10)),3) > 1e-3;
Cn = zeros(size(C));
theta = 0:179;
mu_water=0.19; % attenuation coefficient of water
scale = 1e12; % scaling factor for poisson noise

for t = 1 : size(C,3)
    F = C(:,:,t);
    R = radon(F,theta);
    mu = R/1000*mu_water+mu_water; % transform from CT number to attenuation coefficient
    I=I0*exp(-mu);
%     In=scale*imnoise(I/scale,'poisson')+randn(size(mu))*sigma;
    In = poissrnd(I) + randn(size(mu))*sigma;
%     In = I;
    mu_n=-log(In/I0);
    Rn=1000*(mu_n-mu_water)/mu_water; % tranform back from attenuation coefficient to CT number
    Fn=iradon(Rn,theta,'linear','Shepp-Logan',1,size(F,1));
    Fn(~Mask)=0;
    
    Cn(:,:,t) = Fn;
end
end