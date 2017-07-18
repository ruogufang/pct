function [out] = mrp_segment(in, low, high, frames)
%MRP_SEGMENT Segments out soft tissue and bone from a PCT map
%   
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%
%   USAGE:  OUT = MRP_SEGMENT(IN, LOW, HIGH, FRAMES)
%   
%   PRE:
%       IN     - a MR series [T x X x Y x Z]
%       LOW    - The lower cutoff value [Scalar]
%       HIGH   - The higher cutoff value [Scalar]
%       FRAMES - No. of frames to calculate the average from [Scalar]
%
%   POST:
%       OUT    - A MRP map like IN but with voxels below LOW and above HIGH
%                segmented out (changed to zero)
%
%   This function calculates the average of the first FRAME frames of the
%   CT series. All voxels with value below LOW are set to zero and all
%   voxels with value above HIGH are set to zero.

if nargin < 4
    frames = 10;
end

%Calculate the average of the first frames
avg = squeeze(mean(in(1:frames,:,:,:)));

%Find the indexes of the tissue to be segmented out
[ind_lo] = find(avg < low);
[ind_hi] = find(avg > high);

%Set chosen voxels to zero
in(ind_lo) = 0;
in(ind_hi) = 0;

%All done
out = in;

end

