function [ lsnr ] = pct_lsnr(in)
%PCT_RMSE computes the root mean-square-error between in and ref
%
%   Ruogu Fang 04/16/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  LSNR = PCT_LSNR(IN);
%
%   PRE:
%       IN      - An input signal of N-dimension [N-D]
%
%   POST:
%       LSNR    - Local signal to noise ratio
%
%

lsnr = mean(in(:))/std(in(:));

end

