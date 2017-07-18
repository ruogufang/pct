function [ rmse ] = pct_rmse(in, ref, mask)
%PCT_RMSE computes the root mean-square-error between in and ref
%
%   Ruogu Fang 03/14/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  RMSE = PCT_RMSE(IN, REF, MASK);
%
%   PRE:
%       IN      - An input signal of N-dimension [N-D]
%       REF     - A reference signal of same size as IN [N-D]
%       MASK    - A mask for valid pixels [N-D] (optional)
%
%   POST:
%       RMSE    - Root mean-square-error between input signal IN and reference
%
%

if nargin < 3
    mask = true(size(ref));
end

rmse = sqrt(mean((in(mask(:))-ref(mask(:))).^2));

end

