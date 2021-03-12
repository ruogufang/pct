function pct_vol2gif_demo
% Convert dicom files to GIF files
files = dir(cd);
files = files(4:end);

for i = 1 : length(files)
    load(files(i).name);
    V = permute(double(squeeze(V(:,:,1,:))),[3 1 2]);
    name = ['../GIF/' files(i).name(1:end-4)];
    pct_vol2gif(V,name);
end

end