function [x, y] = pct_tacplot(data)
%Plot time attenuation curve of data
% Ruogu Fang 06/14/2012

T = size(data, 1);
figure(1);
imshow(squeeze(data(floor(T/2),:,:)),[0 100]);
[x y] = ginput(1);
x = round(x)
y = round(y)
while x > 0 && y > 0
    tac = data(:,y,x);
    figure(2);
    plot(tac);
    figure(1);
    [x y] = ginput(1);
    x = round(x)
    y = round(y)
end


