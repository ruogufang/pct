function [ K ] = mrp_bsvd(inmap,aif,dt,lambda,m,mask)
%MRP_BSVD Deconvolution of the Indicator-Dilution equation using bSVD by
%circular deconvolution and circulating back the last values
%
%   USAGE:  K = MRP_BSVD(INMAP,AIF,DT,LAMBDA,M,MASK)
%
%   INPUT:
%       INMAP   - A [T x X x Y x Z] matrix
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
%                 blood flow. [T x X x Y x Z].
%
%   This function solves the Indicator-Dilution equation
%
%       C = F * conv( C_a,R )
%
%   using block-circulant SVD. See "Tracer Arrival Timing-Insensitive Technique
%   for Estimating Flow in MR Perfusion-Weighted Imaging Using Singular Value
%   Decomposition With a Block-Circulant Deconvolution Matrix" by Ona Wu,
%   Leif ?stergaard, Robert M. Weisskoff, Thomas Benner, Bruce R. Rosen, and
%   A. Gregory Sorensen for more detail.
%
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%
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
S(S < lambda*maxS) = 0; %Truncate the lowest singular values
invD = V * pinv(S) * U';

K = zeros(N,X,Y,Z);

%Deconvolve
for i = 1:X
    for j = 1:Y
        for k = 1:Z
            if mask(i,j,k)
                R = invD*[inmap(:,i,j,k); zeros(N-T,1)];
                if R(end) > R(1) % there is bolus delay
                    [~,idx] = max(R);
                    delay = N - idx + 1;
                    R = circshift(R,delay);
                end
                K(:,i,j,k) = R;
            end
        end
    end
end

K = K(1:T,:,:,:);

%Correct sampling rate
K = K/dt;

end

