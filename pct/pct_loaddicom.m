function [ data ] = pct_loaddicom( dirname, nSlices )
%PCT_LOADDICOM Loads dicom images from a directory into a variable
%
%   USAGE: DATA = PCT_LOADDICOM(DIRNAME);
%
%   INPUT:
%       DIRNAME     - A string containing the path to the directory containing
%                     the desired DICOM images. The directory should contain
%                     only DICOM images.
%       nSlices     - A number of slices in DICOM images [Scalar]
%
%   OUTPUT:
%       DATA        - A variable with the data contained in the DICOM files.
%
%   This function is still very basic, but I will hopefully expand it at some
%   point. It assumes four slices and that the files are arranged in temporal
%   order and slice order. Later on I will expand this to automatically detect
%   number of slices, the correct order, size of the images, and temporal
%   information.
%
%   Kolbeinn Karlsson, 08/13/12
%   Ruogu Fang 01/29/13 Revise
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if ~strcmp(dirname(end),'/')
    dirname = [dirname '/'];
end

if nargin < 2
    nSlices = 4;
end

dirlist = dir(dirname);
dirlist = dirlist(3:end);
% nSlices = 4;
nFiles = length(dirlist);
nLoops = nFiles/nSlices;
imageHeight = 512;
imageWidth  = 512;
rescaleIntercept = -1024;

%Pre-allocate the output variable
data = zeros(imageHeight,imageWidth,nSlices,nLoops);

%Load the data from the DICOM files
for i = 1:nLoops
    for j = 1 : nSlices
        k = nSlices*(i-1);
        data(:,:,j,i) = dicomread([dirname dirlist(k+j).name]) + rescaleIntercept;
    end
end

end

