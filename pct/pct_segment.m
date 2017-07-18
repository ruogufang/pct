function [out] = pct_segment(in, low, high, frames)
%PCT_SEGMENT Segments out soft tissue and bone from a PCT map
%   
%   Kolbeinn Karlsson 06/04/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUT = PCT_SEGMENT(IN, LOW, HIGH, FRAMES)
%   
%   PRE:
%       IN     - a CT series in Hounsfield Units [T x X x Y]
%       LOW    - The lower cutoff value [Scalar]
%       HIGH   - The higher cutoff value [Scalar]
%       FRAMES - No. of frames to calculate the average from [Scalar]
%
%   POST:
%       OUT    - A CT map like IN but with voxels below LOW and above HIGH
%                segmented out (changed to zero)
%
%   This function calculates the average of the first FRAME frames of the
%   CT series. All voxels with value below LOW are set to zero and all
%   voxels with value above HIGH are set to zero.

%Calculate the average of the first frames
avg = squeeze(mean(in(1:frames,:,:),1));

%Find the indexes of the tissue to be segmented out
[row_lo,col_lo] = find(avg < low);
[row_hi,col_hi] = find(avg > high);
row = [row_lo; row_hi];
col = [col_lo; col_hi];

%Set chosen voxels to zero
for i = 1:length(row)
    in(:,row(i),col(i)) = 0;
end

%All done
out = in;

end

