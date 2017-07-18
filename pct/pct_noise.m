%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Cn,noise] = pct_noise(C, ACF, SIGMA, TYPE, MASK)
%
% Add spectral (correlated Gaussian)/Gaussian noise to tissue
% time enhancement curve
%
% INPUT:
%       C         - CTP volume [X x Y x T]
%       ACF       - Noise auto-correlation function
%       sigma     - Standard deviation of the noise
%       type      - Type of the noise
%                       'g' - Gaussian Noise
%                       's' - Spectrum Noise (default)
%       mask      - Indicate the valid pixels to add noise
%
% OUTPUT:
%       Cn        - The noisy data [X x Y x T]
%       noise     - The added noise [X x Y x T]
%
% Ruogu Fang
% 7/2/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Cn, noise] = pct_noise(C,acf,sigma,type,mask)

[X,Y,T] = size(C);

%Read Inputs
if nargin < 2 || isempty(acf)
    acf = [];
end
if nargin < 3
    sigma = 1;
end
if nargin < 4
    if ~isempty(acf)
        type = 's';
    else
        type = 'g';
    end
end
if nargin < 5
    mask = ones(X,Y);
end

Cn = zeros(size(C));

if strcmp(type,'s') && ~isempty(acf)
    for t = 1 : size(C,3)
        noise = randn(X,Y);
        %         noise = conv2(noise,acf,'same');
        noise = imfilter(noise,acf,'symmetric','conv');
        noise = sigma * (noise - mean(noise(:)))/std(noise(:));
        
        noise(~mask)=0;
        Cn(:,:,t) = C(:,:,t) + noise;
    end
elseif strcmp(type,'g')
    Cn = C + randn(X,Y,T) * sigma;
else
    error('Unknown noise type.');
end

end