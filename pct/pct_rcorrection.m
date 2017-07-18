function [ rout ] = pct_rcorrection(r)
%PCT_RCORRECTION Corrects the residue function so that it is physically sound
%
%   USAGE:  ROUT = PCT_RCORRECTION(R);
%
%   INPUT:
%       R   - The residue function [T x Y x X]
%
%   OUTPUT:
%       ROUT- The corrected residue function [T x Y x X]
%
%   Once the residue function hits zero, it stays zero. Strictly speaking, the
%   residue function is a non-negative, decreasing function of time.
%
%   Kolbeinn Karlsson, 08/13/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

[T H W] = size(r);

for i = 1:H
    for j = 1:W
        i_0 = find(r(:,i,j)<=0,1,'first');
        r(i_0:end,i,j) = 0;
    end
end

rout = r;

end

