function [ K ] = mrp_tikh(inmap,aif,dt,lambda,m,mask)
%MRP_TIKH Deconvolution of the Indicator-Dilution equation using Tikhonov
%regularization
%
%   USAGE:  K = MRP_TIKH(INMAP,AIF,DT,LAMBDA,M,MASK)
%
%   INPUT:
%       INMAP   - A [T x X x Y x Z] tensor
%       AIF     - The Arterial Input Function [T x 1]
%       DT      - The sampling interval in seconds [Scalar]
%       LAMBDA  - Truncation parameter for the SVD. It denotes the fraction of
%                 the lowest singular values that are set to zero. [Scalar]
%       M       - How often to extend the input. If input is of length T, then
%                 it will be zeropadded to length MT. [Scalar]
%       MASK    - A logical [X x Y x Z] mask. The computation will only be performed
%                 for the voxels that are logical 1.
%
%   OUTPUT:
%       K       - BF * R. The impulse residue function scaled by the cerebral
%                 blood flow. [T x X x Y x Z].
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
%   Ruogu Fang 11/06/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University


% %Temporal interpolation
% if dt ~= 1
%     inmap = pct_tinterp(inmap,dt);
%     aif = pct_tinterp(aif,dt);
% end

%Get the input dimensions
[T,X,Y,Z] = size(inmap);

if nargin<6
    mask = ones(X,Y,Z);
end

%Extend the AIF by zero-padding
N = m*T;
Ca = zeros(N,1);
Ca(1:T) = aif;

%Create the block-circulant matrix
D = gallery('circul',Ca)'; 

%Calculate the inverse of D using SVD
[U, S, V] = svd(D);
maxS = max(diag(S));
diag_new = diag(S)./(diag(S).^2./(diag(S).^2+(lambda*maxS).^2)); % Tikhonov regularization
S = diag(diag_new);
invD = V * pinv(S) * U';

K = zeros(size(inmap));

%Deconvolve
for i = 1:X
    for j = 1:Y
        for k = 1 : Z
        if mask(i,j,k)
            K(:,i,j,k) = invD*[inmap(:,i,j,k); zeros(N-T,1)];
        end
        end
    end
end

K = K(1:T,:,:,:);

%Correct sampling rate
K = K/dt;

end

