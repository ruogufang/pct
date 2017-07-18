function [ snr ] = pct_snr(in, ref, mask)
%PCT_SNR computes the signal to noise ratio between input and reference
%
%   Ruogu Fang 03/20/2014
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
%       SNR     - Signal-to-noise ratio
%
%   SNR is defined as 10*log10(power(in)/var(noise))
%
if nargin < 3
    mask = true(size(ref));
end

snr = 10*log10(mean(in(mask(:)).^2)/var(in(mask(:))-ref(mask(:))));

end

