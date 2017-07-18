%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [Vn] = pct_noise(V, SIGMA, MASK, type)
%FUNCTION Cn = PCT_NOISE(C,SIGMA,MASK) adds Gaussian noise to the Concentration Tissue Curves
%
% INPUT:
%       V         - MRP volume [T x X x Y x Z]
%       SIGMA     - Standard deviation of the noise
%       MASK      - Indicate the valid pixels to add noise
%
% OUTPUT:
%       Vn        - The noisy data [T x X x Y x Z]
%
% The function adds Gaussian noise to the CTC with a random variable
% Z(t)~N(0,sigma_v(t)^2)  with sigma_v^2(t) \propto
% \sigma^2(1/A_v(t)^2+1/A_v(0)^2)
% Same as add Gaussian noise with standard deviation of sigma to the real and
% imaginary part of the MR signal
%
%   Ruogu Fang 11/06/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vn = mrp_noise(V,sigma,mask,type)

thresh = 0.2;

[T,X,Y,Z] = size(V);

%Read Inputs
if nargin < 3
    mask = ones(X,Y,Z);
end

if nargin < 4
    type = 'rician';
end

Vn = zeros(size(V));

switch type
    case 'rician'
        for i = 1 : X
            for j = 1 : Y
                for k = 1 : Z
                    if mask(i,j,k)
                        C = squeeze(V(:,i,j,k));
                        if any(C) && length(unique(C))~=1 % if C is not same for all time points
                            % find the end time point of baseline
                            t0 = find(C(1:round(length(C)/2))<thresh*max(C(1:length(C))),1,'last');
                            S0 = mean(C(1:t0));
                        else
                            S0 = 0;
                        end
                        for t = 1 : T
                            sigma_v = sigma*sqrt(1/(C(t)+1)^2+1/(S0+1)^2);
                            Vn(t,i,j,k) = C(t) + randn(1)*sigma_v;
                        end
                    end
                end
            end
        end
    case 'gaussian'
        Vn(mask) = V(mask) + sigma * randn(size(V(mask)));
    otherwise
        display('Unknown noise type.');
end

end