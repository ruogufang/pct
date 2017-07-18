function [ k, aif ] = pct_bd(inmap,aif_init,dt,maxiter,inneriter,mask)
%PCT_BD Blind Deconvolution of the Indicator-Dilution equation using
%Richardson-Lucy algorithm
%
%   USAGE:  [K,AIF] = PCT_BD(INMAP,AIF_INIT,DT,MASK)
%
%   INPUT:
%       INMAP       - A [T x Y x X] matrix
%       AIF_INIT    - The Initial Arterial Input Function [T x 1]
%       DT          - The sampling interval in seconds [Scalar]
%       MAXITER     - Maximum numbere of outer iteration
%       INNERITER   - Inner iteration number of Lucy-Richardson algorithm
%       MASK        - A logical [Y x X] mask. The computation will only be performed
%                   for the voxels that are logical 1.
%
%   OUTPUT:
%       K       - BF * R. The impulse residue function scaled by the cerebral
%                 blood flow. [T x Y x X].
%       AIF     - Arterial Input Function after blind deconvolution [T x 1]
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
%   Ruogu Fang, 10/20/2014
%   SCIS, FIU

%We will solve the equation Ab=c, where A=dt*conv_matrix(aif), b=BF*R, and c
%is the TAC for that pixel.

% %Temporal interpolation
% if dt ~= 1
%     inmap = pct_tinterp(inmap,dt);
%     aif = pct_tinterp(aif,dt);
% end

if nargin<3
    dt = 1;
end
if nargin<4
    maxiter = 10;
end
if nargin<5
    inneriter = 10;
end
if nargin<6
    mask = ones(size(inmap,2),size(inmap,3));
end

%Pre-allocate the output
[l,h,w] = size(inmap);
k = zeros([l h w]);
aif = zeros([l h w]);

% Richardson-Lucy algorithm
for i = 1 : h
    for j = 1 : w
        if mask(i,j)
            a0 = aif_init;
            r0 = ones(l,1)*sum(inmap(:,i,j))/l;
            af = a0;
            rf = r0;
            for t = 1 : maxiter
                as = af;
                rs = rf;
                c = squeeze(inmap(:,i,j));
                
                B = max(as,0);% PSF positivity constraint
                sumAIF = sum(B(:));
                B = B/(sumAIF + (sumAIF==0)*eps);% normalization is a necessary constraint,
                as = B;
                
                %                 a. Update a
                for tt = 1 : inneriter
                    a = conv((c./conv(as,rf,'same')),rf(end:-1:1),'same').*as;
                    as = a;
                    plot(a,'r');
                    hold on;
                end
                af = a;
                
                %                 b. Update r
                for tt = 1 : inneriter
                    r = conv((c./conv(af,rs,'same')),af(end:-1:1),'same').*rs;
                    rs = r;
                    plot(r);
                end
                
                rf = r;
                
                %                 % Convolution in frequency using FFT
                %                 % a. Update r
                %                 H = psf2otf(as,size(a0));
                %                 scale = real(ifftn(conj(H))) + sqrt(eps);
                %                 rf = max(c.*real(ifftn(conj(H)))./scale,0);
                %                 clear scale;
                %
                %                 % 2. Update a
                %                 H = fftn(rs);
                %                 scale = otf2psf(conj(H),size(a0)) + sqrt(eps);
                %                 af = max(B.*otf2psf(conj(H),size(a0))./scale,0);
                %                 clear CC H;
                
                % c. apply positivity constraint
                af = max(af,0);
                rf = max(rf,0);
                
                sumAIF = sum(af(:));
                af = af/(sumAIF + (sumAIF==0)*eps);
                
            end
            k(:,i,j) = rf;
            aif(:,i,j) = af;
        end
    end
end

k = k/dt;
aif = aif/dt;

end

