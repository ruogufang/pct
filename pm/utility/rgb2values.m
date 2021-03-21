function image_values = rgb2values(image, colormap, modality)
% Function Header
% Description:
%   Converts an RGB-image to the values which correspond to the colors
%       based on a colormap that is passed in.
%
% Inputs:
%   Image : An RGB image of any row and col size, the color pixel must be 
%       RGB. That is, [x][y][1] must corresond to red component of the pixel, [x][y][2] must
%       correspond to green component of the pixel and [x][y][3] must correspond to the blue
%       component of the pixel
%
%   Colormap : Must be a color look up table, any colormap in
%       select_colormap() will work, as well as the colormaps from running
%       load('Rainbow_CBP')
%
%   Modality : Modality of the image which is being passed in. Values can
%       be editted in the switch - case area. Can also be used to generate
%       indexed image.
%
% Output:
%   image_values : returns a 2-D image of the same size of the original
%       image columns and rows.
%
% Written by : Simon Kato
%              Smile-LAB @UF

% Function
switch modality
    case 'CBF'
        maxV = 60;
        minV = 0;
    case 'CBV'
        maxV = 30 %4;
        minV = 0;
    case 'MTT'
        maxV = 12;
        minV = 0;
    case 'TTP'
        maxV = 25;
        minV = 0;
    case 'DLY'
        maxV = 10;
        minV = 0;
    case 'scale'
        maxV = 1;
        minV = 0;
    case 'gray' %Colormaps that have look-up-table 0-255
        maxV = 255;
        minV = 0;
    case 101 %Colormaps that have look-up table 0-100
        maxV = 100;
        minV = 0;
    case 201 %Colormaps that have look-up table 0-200
        maxV = 200;
        minV = 0;
    case 401 %Colormaps that have look-up table 0-400 ie. Rainbow_4
        maxV = 400;
        minV = 0;
    otherwise
        ME = MException('MyComponent:noSuchModality', ...
            'modality %s not found',modality);
        throw(ME)
end
        
Red = image(:,:,1);
Green = image(:, :, 2);
Blue = image(:, :, 3);

colormap_Red = (colormap(:,1)*256);
colormap_Green = (colormap(:,2)*256);
colormap_Blue = (colormap(:,3)*256);

dims = size(image); x = dims(2) ;y = dims(1);

image_values = zeros(y,x);

scaler = size(colormap,1)-1;

for i = 1:y
    for j = 1:x
        if Red(j,i) == 0 && Green(j,i) == 0 && Blue(j,i) == 0 %Meant to increase performance
            image_values(j,i) = 0;
        else
            if isempty(find(colormap_Red == Red(j,i) & colormap_Green == Green(j,i) & colormap_Blue == Blue(j,i))) %#ok<EFIND>
                values = [colormap_Red - double(Red(j,i))*ones(scaler + 1,1), colormap_Green - double(Green(j,i))*ones(scaler + 1,1), colormap_Blue - double(Blue(j,i))*ones(scaler + 1,1)]';
                [~, closestValue] = min(vecnorm(values));
                image_values(j,i) = (closestValue-1)/scaler*(maxV-minV) + minV;
            else
                image_values(j,i) = ((find(colormap_Red == Red(j,i) & colormap_Green == Green(j,i) & colormap_Blue == Blue(j,i)) - 1)/scaler)*(maxV-minV) + minV;
            end
        end
    end
end
end
