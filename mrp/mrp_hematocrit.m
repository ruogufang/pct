function [ out ] = mrp_hematocrit(in, k_H)
%MRP_HEMATOCRIT Corrects a CT map for hematocrit differences.
%
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%
%   USAGE:  OUTMAP = MRP_HEMATOCRIT(CTMAP, k_H);
%
%   PRE:
%       IN      - A contrast concentration map / MR map [X x Y x Z x T]
%       k_H     - Hematocrit correction factor
%
%   POST:
%       OUT     - A contrast concentration map corrected for hematocrit
%                 differences
%
%   This function simply multiplies each pixel point-wise with k_H.

out = k_H .* in;

end

