function [ psnr ] = pct_psnr(in, ref, mask)
%PCT_PSNR computes the peak signal to noise ratio (PSNR) between input and reference
%
%   Ruogu Fang 05/03/2014
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  SNR = PCT_SNR(IN, REF, MASK);
%
%   PRE:
%       IN      - An input signal of N-dimension [N-D]
%       REF     - A reference signal of same size as IN [N-D]
%       MASK    - A mask for valid pixels [N-D] (Logical, optional)
%
%   POST:
%       PSNR    - Peak signal-to-noise ratio
%
%   SNR is defined as 10*log10(power(in)/var(noise))
%

if nargin < 3
    mask = true(size(ref));
end

psnr = 20*log10(max(ref(mask(:)))/sqrt(mean((in(mask(:))-ref(mask(:))).^2)));

end

