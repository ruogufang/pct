function [ k ] = pct_bsvd(inmap,aif,dt,lambda,m,mask)
%PCT_BSVD Deconvolution of the Indicator-Dilution equation using bSVD
%
%   USAGE:  K = PCT_BSVD(INMAP,AIF,DT,LAMBDA,M,MASK)
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
%   using block-circulant SVD. See "Tracer Arrival Timing-Insensitive Technique
%   for Estimating Flow in MR Perfusion-Weighted Imaging Using Singular Value
%   Decomposition With a Block-Circulant Deconvolution Matrix" by Ona Wu,
%   Leif ?stergaard, Robert M. Weisskoff, Thomas Benner, Bruce R. Rosen, and 
%   A. Gregory Sorensen for more detail.
%
%   Kolbeinn Karlsson, 08/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

% %Temporal interpolation
% if dt ~= 1
%     inmap = pct_tinterp(inmap,dt);
%     aif = pct_tinterp(aif,dt);
% end

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
S(S < lambda*maxS) = 0; %Truncate the lowest singular values
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

