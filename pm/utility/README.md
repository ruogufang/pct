# `rgb2values`

## Description:
    rgb2values - Takes a RGB, [x][y][3], image tied to a colormap and produces an 2-D array with the values at each pixel.
   ### Specifics:
    If an image with values has a colormap applied to it, it generates a RGB image which indicate those colors. The reason for 
    this code is if the image was somehow turned into the RGB version of the image. That is, the image that is had is the RGB 
    image corresponding to the values indicated by the range of colormap. This code takes  the RGB image and the colormap which 
    was used to generate the image and uses the colormap to "reverse engineer" the values which would've given the color at each 
    pixel. Since vector operations are quick in matlab, what was done was separate the colorbar into R, G, and B and separate each 
    pixel into R, G, and B. Then, the minimum norm of the R and the value of red at a pixel was calculated. This was used to then 
    normalize the image between 0-1. Afterwards, given specific bounds on modalities, was scaled so that max would be the max of 
    the modality and min was 0.
    
## Parameters:
```
    --image        An [x][y][3] image where the last array of size 3 corresponds to [R, G, B]
    --colormap     A loaded colormap, for using Rainbow_CBP_4, load('Rainbow_CBP_1.mat');
    --modality     The modality or indexing options which have the following min-max:
        'rCBF' : [0, 60]
            "Blood Flow"
        'rCBV' : [0, 4]
            "Blood Volume"
        'MTT' : [0, 12]
            "Mean Transit Time"
        'DLY' : [0, 10] 
            "Delay"
        'TTP' : [0, 25]
            "Time to Peak"
        'scale' : [0, 1]
            "Scales the values between 0 to 1
        'gray' : [0,255]
            "Colormaps that have look-up-table 0-255
 ```
        
## Example:
    value_image = rgb2values(rgb_image, my_colormap, 'rCBF');
## Note:
    How to reapply a color map for the sake of verifying the colors are the same:
    (1) You must create an indexed image according to your colormap. To make an indexed image, this is easily done with the code: 
    see how many rows there are in the colormap. Subtract one from the colormap, within the switch-case add the num_rows - 1 as a 
    modality, make maxV = num_rows - 1 and minV = 0, calling rgb2values(image, colormap, numRows-1) will produce a correctly indexed image for the colormap.
    (2) Call ind2rgb(image_created, colormap)
    (*) You may need to cast your rgb2values to int8,int16, double, or some other number type before calling ind2rgb, I recommend using matlab im2double()
    
## Contribution and Contanct:
    SMILE LAB @UF
    By: Simon Kato 
