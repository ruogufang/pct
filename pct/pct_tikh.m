function [ k ] = pct_tikh(inmap,aif,dt,lambda,m,mask)
%PCT_TIKH Deconvolution of the Indicator-Dilution equation using Tikhonov
%regularization
%
%   USAGE:  K = PCT_TIKH(INMAP,AIF,DT,LAMBDA,M,MASK)
%
%   INPUT:
%       INMAP   - A [T x Y x X] matrix
%       AIF     - The Arterial Input Function [T x 1]
%       DT      - The sampling interval in seconds [Scalar]
%       LAMBDA  - Truncation parameter for the SVD. It denotes the fraction of
%                 the lowest singular values that are set to zero. [Scalar]
%       M       - How often to extend the input. If input is of length T, then
%                 it will be zeropadded to length MT. [Scalar]
%       MASK    - A logical [Y x X] mask. The computation will only be performed
%                 for the voxels that are logical 1.
%
%   OUTPUT:
%       K       - BF * R. The impulse residue function scaled by the cerebral
%                 blood flow. [T x Y x X].
%
%   This function solves the Indicator-Dilution equation
%
%       C = F * conv( C_a,R )
%
%   using Tikhnov regularization. 
%   k_lambda = \sum_1^r(f_{lambda,i} * u'_i*c*v_i/sigma_i)
%   f_{lambda,i}=sigma_^2/(sigma_i^2+lambda^2)
%
%   See "Assessment of Perfusion by Dynamic Contrast-Enhanced Imaging Using
%   a Deconvolution Approach Based on Regression and Singular Value
%   Decomposition T. S. Koh*, X. Y. Wu, L. H. Cheong, and C. C. T. Lim
%
%   Ruogu Fang, 09/29/13
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

% %Temporal interpolation
% if dt ~= 1
%     inmap = pct_tinterp(inmap,dt);
%     aif = pct_tinterp(aif,dt);
% end

if nargin < 6
    mask = ones(size(inmap));
end

%Get the input dimensions
[L, height, width] = size(inmap);

%Extend the AIF by zero-padding
N = m*L;
Ca = zeros(N,1);
Ca(1:L) = aif;

%Create the block-circulant matrix
D = gallery('circul',Ca)'; 

%Calculate the inverse of D using SVD
[U, S, V] = svd(D);
maxS = max(diag(S));
diag_new = diag(S)./(diag(S).^2./(diag(S).^2+(lambda*maxS).^2)); % Tikhonov regularization
S = diag(diag_new);
invD = V * pinv(S) * U';

k = zeros(N,height,width);

%Deconvolve
for i = 1:height
    for j = 1:width
        if mask(i,j)
            k(:,i,j) = invD*[inmap(:,i,j); zeros(N-L,1)];
        end
    end
end

k = k(1:L,:,:);

%Correct sampling rate
k = k/dt;

end

