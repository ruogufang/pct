function [ outmap ] = mrp_truncate(inmap, pre, post)
%MRP_TRUNCATE Truncates a MRP map (in the time-dimension)
%
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University

%
%   USAGE:  OUTMAP = MRP_TRUNCATE(INMAP, PRE, POST);
%
%   PRE:
%       INMAP   - A PCT map [T x X x Y x Z]
%       PRE     - The first frame to include. 1<= PRE < T
%       POST    - The last frame to include. PRE < POST <= T.
%
%   POST:
%       OUTMAP  - A PCT series containing the PRE:POST frames of INMAP.
%                 [(POST-PRE+1) x X x Y x Z].
%

outmap = inmap(pre:post,:,:,:);


end

