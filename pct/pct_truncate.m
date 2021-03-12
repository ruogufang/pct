function [ outmap ] = pct_truncate(inmap, pre, post)
%PCT_TRUNCATE Truncates a  PCT map (in the time-dimension)
%
%   Kolbeinn Karlsson 06/06/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUTMAP = PCT_TRUNCATE(INMAP, PRE, POST);
%
%   PRE:
%       INMAP   - A PCT map [T x Y x X]
%       PRE     - The first frame to include. 1<= PRE < T
%       POST    - The last frame to include. PRE < POST <= T.
%
%   POST:
%       OUTMAP  - A PCT series containing the PRE:POST frames of INMAP.
%                 [TT x Y x X].
%

outmap = inmap(pre:post,:,:);


end

