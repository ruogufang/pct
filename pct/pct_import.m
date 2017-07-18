%pct_import.m
%Kolbeinn Karlsson 05/29/12
%Advanced Multimedia Processing (AMP) Lab, Cornell University
%This script imports the images from the folder DICOM and 
%divides them into four series, one for each slice. The images
%are rescaled to Hounsfield Units (HU).
clear all; close all; clc;

nslice = 4;
raw_list = dir(cd);
list = raw_list(3:end);
file_name = list(1).name;
series_length = floor(length(list)/nslice);

info = dicominfo(file_name);
V = zeros(info.Height, info.Width, nslice, series_length,'int16');

for i = 0:series_length-1
    for j = 1 : nslice
        V(:,:,j,i+1) = info.RescaleSlope * dicomread(list(j+nslice*i).name) + info.RescaleIntercept;
    end
end

clear raw_list list file_name i j
