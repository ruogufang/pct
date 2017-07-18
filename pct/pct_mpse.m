function [ mpse ] = pct_mpse(in, ref)
%PCT_RMSE computes the mean percent square error between in and ref
%
%   Ruogu Fang 04/16/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  MPSE = PCT_MPSE(IN, REF);
%
%   PRE:
%       IN      - An input signal of N-dimension [N-D]
%       REF     - A reference signal of same size as IN [N-D]
%
%   POST:
%       MPSE    - Mean percent square error between input signal IN and reference
%
%   MPSE = 100/mean(ref)*sqrt(1/(Q-1)*sum_1^Q((in-ref).^2))
%
%

mpse = 100/mean(ref(:))*sqrt((sum(in(:)-ref(:)).^2)/(numel(in)-1));

end

