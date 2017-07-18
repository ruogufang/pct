function [ k ] = pct_dsvd(inmap,aif,dt,lambda,mask)
%PCT_DSVD Solves the Indicator-Dilution equation using delay-corrected SVD
%
%   USAGE:  K = PCT_DSVD(INMAP,AIF,DT,LAMBDA,MASK)
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
%   using delay-corrected SVD. The delay correction consints of 
%       t_d = t_aif - t_c
%   where t_aif is the bolus arrival time of the AIF and t_c is the bolus
%   arrival time of the voxel of interest. The time-attenuation curve (TAC) for
%   each voxel is shifted by t_d before deconvolving. For additional details,
%   see "Difference in Tracer Delay?induced Effect among Deconvolution 
%   Algorithms in CT Perfusion Analysis: Quantitative Evaluation with Digital 
%   Phantoms" by K. Kudo et al.
%
%   Kolbeinn Karlsson, 08/02/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

function [out] = tdelay(in,t)
    if t == 0
        out = in;
    else
        out = [in(1+t:end); zeros(t,1)];
        % out = [zeros(t,1); in(1:end-t)]; % Shift right (Ruogu)
    end
end

%We will solve the equation Ab=c, where A=dt*conv_matrix(aif), b=BF*R, and c 
%is the TAC for that pixel.

%Bolus arrival time threshold is 5%
threshold = 0.05;

%Create A, the convolution matrix of the AIF
A = tril(toeplitz(aif));
%Perform the SVD
[U S V] = svd(A);
%Truncate the singular values
maxS = max(diag(S));
S(S < lambda*maxS) = 0;
%Calculate the pseudoinverse of A
invA = V*pinv(S)*U';

%Pre-allocate the output
[l h w] = size(inmap);
k = zeros([l h w]);

%Compute the delay-correction
delay = zeros(h,w);
for i = 1:h
    for j = 1:w
        if mask(i,j)
            B = inmap(:,i,j);
            bat = find(B > threshold*max(B),1,'first');
            if ~isempty(bat)
                delay(i,j) = bat;
            end
            
        end
    end
end
th_aif = threshold*max(aif);
aif_delay = find(aif > th_aif,1,'first');
delay(delay < aif_delay) = aif_delay;
        
%Go through each pixel and deconvolve
for i = 1:h
    for j = 1:w
        if mask(i,j)
            k(:,i,j) = invA*tdelay(inmap(:,i,j),delay(i,j));
        end
    end
end


end

