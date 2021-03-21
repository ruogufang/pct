function [ out ] = pct_truncatetime(in, pre, post)
%PCT_TRUNCATETIME Truncates a TAC [T x 1] or image sequence [T x Y x X] in time
%
%   USAGE:  OUT = PCT_TRUNCATETIME(IN, PRE, POST);
%
%   PRE:
%       IN      - A time-attenuation curve (TAC) [T x 1] or 
%                 image sequence [T x Y x X]
%       PRE     - The first frame to include. 1<= PRE < T
%       POST    - The last frame to include. PRE < POST <= T.
%
%   POST:
%       OUT     - A TAC [T_new x 1] or image sequence containing the PRE:POST frames of 
%                 INMAP [T_new x Y x X].
%
%   To only truncate the beginning and not the end, set POST=-1.
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if ndims(in) == 3
    if post == -1
        out = in(pre:end,:,:);
    else
        out = in(pre:post,:,:);
    end
elseif ndims(in) == 2
    if post == -1
        out = in(pre:end,:);
    else
        out = in(pre:post,:);
    end
end

