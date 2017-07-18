% Demo for function pct_irb
%
%   Ruogu Fang 07/06/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

clear; clc; close all;

ctp_path = '~/Dropbox/Research/DICOM/CT_MAT/';
irb_path = '~/Dropbox/Research/DICOM/IRB/';
meta = importdata('~/Dropbox/Research/DICOM/CTP_meta.csv');

files = dir(ctp_path);
files = files(4:end);
for i = 1 : length(files)
    filename = files(i).name;
    load([ctp_path filename]);
    irb_file = [irb_path 'IRB_' filename(5:7) '.mat'];
    pct_irb(V,meta.data(i,:),irb_file);
end
