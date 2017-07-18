%mrp_import.m
%Import MR Perfusion Image
% Image source: https://public.cancerimagingarchives.net
% Shared List: TCGA-GBM DSC T2* MR Perfusion
% Subject ID: TCGA-06-0133
%
% File order: S1(T1,...,Tn) S2(T1,...,Tn) Sk(T1,...,Tn)
%
% Ruogu Fang
% 2014-11-05

clear all; close all; clc;

nslice = 32;
raw_list = dir(cd);
list = raw_list(3:end);
file_name = list(1).name;
series_length = floor((length(list))/nslice);

info = dicominfo(file_name);
V = zeros(info.Height, info.Width, nslice, series_length,'int16');

for j = 0 : nslice-1
    for i = 1:series_length
        V(:,:,j+1,i) = info.RescaleSlope * dicomread(list(i+series_length*j).name) + info.RescaleIntercept;
    end
end

clear raw_list list file_name i j

