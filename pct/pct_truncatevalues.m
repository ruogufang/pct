function [ out ] = pct_truncatevalues(in, lo, hi)
%PCT_TRUNCATEVALUES Truncates the values of a matrix
%
%   USAGE:  OUT = PCT_TRUNCATEVALUES(IN, LO, HI);
%
%   PRE:
%       IN      - A matrix of any dimension
%       LO      - The lowest allowed value [Scalar]
%       HI      - The highest allowed value [Scalar]
%
%   POST:
%       OUT     - A matrix of the same dimensions as IN where all the 
%                 values <LO are set to LO and the values >HI are set to HI.
%
%   This function changes all values in IN that are higher than HI to HI and all
%   values that are lower than LO to LO. To only truncate from below, pass Inf
%   as the "hi" argument and to only truncate from above, pass -Inf as the "lo"
%   argument.
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

out = in;

if lo > -Inf
    out(out < lo) = lo;
end

if hi < Inf
    out(out > hi) = hi;
end

end

