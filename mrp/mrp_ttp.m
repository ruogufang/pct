function [ ttpmap ] = mrp_ttp(inmap,dt,mask)
%MRP_TTP Calculates a Time-To-Peak map from a MR perfusion sequence
%
%   USAGE:  TTPMAP = MRP_TTP(INMAP);
%
%   PRE:
%       INMAP   - A MR / Tracer concentration map sequence [T x X x Y x Z]
%       DT      - Time interval between samples in seconds [Scalar]
%
%   POST:
%       TTPMAP  - A Time-To-Peak map [X x Y x Z]
%
%   This function calculates a TTP map from a MR / tracer concentration
%   map. It basically shows you the time index of
%   the peak of each pixel. 
%
%   Ruogu Fang 11/06/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University

[T,X,Y,Z] = size(inmap);

if nargin < 2
    dt = 1;
end
if nargin < 3
    mask = ones(X,Y,Z);
end

[~,ttpmap] = max(inmap,[],1);
    
% ttpmap = zeros(X,Y,Z);
% for i = 1:X
%     for j=1:Y
%         for k=1:Z
%             if mask(i,j,k)
%                 ttpmap(i,j,k) = find(rmap(:,i,j,k)==rmax(i,j,k),1,'last');
%             end
%         end
%     end
% end

ttpmap = squeeze(ttpmap);

%Convert to seconds
ttpmap = dt * ttpmap;
ttpmap(~mask) = 0;

end

