function Ca = pct_bcaif(aif,m)
%PCT_BCAIF generates block-circulant AIF matrix
%
%   Ruogu Fang 07/16/2012
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  [CA] = PCT_BCAIF(AIF,M);
%
%   PRE:
%       AIF     - A time-attenuation curve of a AIF pixel [T x 1]
%       M       - Extend AIF by m times [Scalar]
%
%   POST:
%       Ca      - Block circulant AIF matrix [mT x mT]
%
%   This function generates the block circulant AIF matrix

if nargin < 2;
    m = 1;
end

aif = aif(:);
Naif = length(aif);
N = m * Naif;
Nc = N + Naif - 1;

if N < Naif,
    error('The circular convolution matrix size must be greater or equal the channel length');
end

H = convmtx(aif,N);

Ca = H(1:N, 1:N);

if m > 1
    Ca(1:(Naif-1), 1:N) = Ca(1:(Naif-1), 1:N) + H((N+1):Nc, 1:N);
end

end