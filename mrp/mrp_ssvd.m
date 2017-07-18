function [ K ] = mrp_ssvd(inmap,aif,dt,lambda,mask)
%MRP_SSVD Deconvolution of the Indicator-Dilution equation using standard SVD
%
%   USAGE:  K = MRP_SSVD(INMAP,AIF,DT,LAMBDA,MASK)
%
%   INPUT:
%       INMAP   - A [T x X x Y x Z] tensor
%       AIF     - The Arterial Input Function [T x 1]
%       DT      - The sampling interval in seconds [Scalar]
%       LAMBDA  - Truncation parameter for the SVD. It denotes the fraction of
%                 the lowest singular values that are set to zero.
%       MASK    - A logical [X x Y x Z] mask. The computation will only be performed
%                 for the voxels that are logical 1.
%
%   OUTPUT:
%       K       - BF * R. The impulse residue function scaled by the cerebral
%                 blood flow. [T x X x Y x Z].
%
%   This function solves the Indicator-Dilution equation
%
%       C = F * conv( C_art,R )
%
%   using standard SVD. See "High Resolution Measurement of Cerebral Blood Flow
%   using Intravascular Tracer Bolus Passages. Part I: Mathematical Approach
%   and Statistical Analysis" by Leif  Ostergaard, Robert M. Weisskoff,
%   David A. Chesler, Carsten Gyldensted, Bruce R. Rosen for more details.
%
%   Ruogu Fang 11/06/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University

if nargin < 5
    [T,X,Y,Z] = size(inmap);
    mask = ones(X,Y,Z);
end

%We will solve the equation Ab=c, where A=dt*conv_matrix(aif), b=BF*R, and c
%is the TAC for that pixel.

% %Temporal interpolation
% if dt ~= 1
%     inmap = pct_tinterp(inmap,dt);
%     aif = pct_tinterp(aif,dt);
% end

%Create A, the convolution matrix of the AIF
A = tril(toeplitz(aif));
%Perform the SVD
[U,S,V] = svd(A);
%Truncate the singular values
maxS = max(diag(S));
S(S < lambda*maxS) = 0;
%Calculate the pseudoinverse of A
invA = V*pinv(S)*U';

%Pre-allocate the output
[T,X,Y,Z] = size(inmap);
K = zeros(size(inmap));

%Go through each pixel
for i = 1:X
    for j = 1:Y
        for k = 1 : Z
            if mask(i,j,k)
                K(:,i,j,k) = invA*inmap(:,i,j,k);
            end
        end
    end
end

K = K/dt;

end

