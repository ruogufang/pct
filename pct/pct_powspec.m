function ps = pct_powspec(input,norm)
%PCT_POWERSPECTRUM computes the power spectrum of the input
%
%   USAGE:  PS = PCT_POWSPEC(INPUT,NORM)
%
%   INPUT:
%       INPUT   - A input image or volume
%       NORM    - Flag to indicate whether to normalize the input. 
%                   0-no normalization, 1 - normalize (default)
%
%   OUTPUT:
%       PS      - Power spectrum of the input
%
%   This function computes the power spectrum of the input by FFT.
%
%
%   Ruogu Fang, 10/01/13
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if nargin < 2
    norm = 1;
end

if norm == 1 % normalize image
    input = (input-mean(input(:)))./std(input(:));
end

input_fft = fftshift(fft2(input));
ps = abs(input_fft).^2;

end