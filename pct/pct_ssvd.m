function [ k ] = pct_ssvd(inmap,aif,dt,lambda,mask)
%PCT_SSVD Deconvolution of the Indicator-Dilution equation using standard SVD
%
%   USAGE:  K = PCT_SSVD(INMAP,AIF,DT,LAMBDA,MASK)
%
%   INPUT:
%       INMAP   - A [T x Y x X] matrix
%       AIF     - The Arterial Input Function [T x 1]
%       DT      - The sampling interval in seconds [Scalar]
%       LAMBDA  - Truncation parameter for the SVD. It denotes the fraction of
%                 the lowest singular values that are set to zero.
%       MASK    - A logical [Y x X] mask. The computation will only be performed
%                 for the voxels that are logical 1.
%
%   OUTPUT:
%       K       - BF * R. The impulse residue function scaled by the cerebral
%                 blood flow. [T x Y x X].
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
%   Kolbeinn Karlsson, 08/02/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%We will solve the equation Ab=c, where A=dt*conv_matrix(aif), b=BF*R, and c
%is the TAC for that pixel.

% %Temporal interpolation
% if dt ~= 1
%     inmap = pct_tinterp(inmap,dt);
%     aif = pct_tinterp(aif,dt);
% end

if nargin < 5
    mask = ones(size(inmap));
end

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
[l,h,w] = size(inmap);
k = zeros([l h w]);

%Go through each pixel
for i = 1:h
    for j = 1:w
        if mask(i,j)
            k(:,i,j) = invA*inmap(:,i,j);
        end
    end
end

k = k/dt;

end

