%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Sn = pct_noisesino(S,I0,sigma)
%
% Simulate low-dose recontructed CTP data by adding Poisson and Gaussian
% noise to the sinogram generated from fanbeam
%
% INPUT:
%   S       - Input high-dose Sinogram [T x X x Y]
%   I0      - Incident X-ray intensity, or incident photon number (10^6 for normal-dose, 2.5*10^5 for 1/7 of normal dose)
%   sigma   - Standard deviation of Gaussian noise in low-dose (0 for normal-dose, 10 for 1/7 of normal dose)
%
% OUTPUT:
%   Sn      - Low-dose sinogram [T x X x Y]
%
% Ruogu Fang
% 07/27/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Sn = pct_noisesino(S, I0, sigma)

S = double(S);
T = size(S,1);
Sn = zeros(size(S));
mu_water=0.19; % attenuation coefficient of water
% scale = 1e12; % scaling factor for poisson noise

for t = 1 : T
    s = squeeze(S(t,:,:));
    mu = s/1000*mu_water+mu_water; % transform from CT number to attenuation coefficient
    I=I0*exp(-mu);
%     In=scale*imnoise(I/scale,'poisson')+randn(size(mu))*sigma;
    In = poissrnd(I) + randn(size(mu))*sigma;
    mu_n=-log(In/I0);
    sn=1000*(mu_n-mu_water)/mu_water; % tranform back from attenuation coefficient to CT number
    Sn(t,:,:) = sn;
end
end