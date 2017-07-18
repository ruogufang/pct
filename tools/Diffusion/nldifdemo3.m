% NLDIFDEMO3   Non linear diffusion demonstration 3
%
%  Type NLDIFDEMO3.
%
clc

disp('_____________________________Non Linear diffusion demonstration  3________________________')
disp(' ')
disp('  Now let`s see the nonlinear diffusion in color images. There are two ways to diffuse a ')
disp('color image: 1) each channel can be diffused separately, as if there were three different images')
disp('or 2) the gradient can be computed by a norm of the gradients in each channel, so that the same')
disp('gradient and diffusivity are used for all channels.')
disp('  Let`s start with the image of a house. The objective here is to simplify this image so that we')
disp('get something like a picewise constant image (this is ideal for a latter segmentation).')
disp(' ')
disp('Press any key to continue...')
pause

%$$$$$ y = nldifc(im,8,1,20,20,50,2,1,'aos', 'grad', 'dfstep', 5, 'alt1','imscale');

im1_path = which('dif_house.jpg');
im1 = imread(im1_path);
figure(1)
colormap(gray(256))
clf
image(im1)
title('Original Image')

disp(' ')
disp(' ')
disp('  Let`s use the 1st method of color diffusion: each channel individually.')
disp(' ')
disp('Press any key to continue...')
pause

y = nldif(rgb2gray(im1),8,0,10,0,1,1,1,'aos', 'grad');
figure(1)
subplot(1,2,1)
image(im1)
title('Original Image')
subplot(1,2,2)
image(im1)
title('Nonlinear diffusion')
figure(2)
colormap(hot)

disp(' ')
disp(' ')
disp('  Adjust the windows to better see the diffusion.')
disp(' ')
disp('Press any key to continue...')
pause

y = nldifc(im1,5,1,20,20,6,1,2,'aos', 'grad', 'dfstep', 5);

disp(' ')
disp(' ')
disp('  In the 2nd method, the difussivity is calculated based on a norm of the gradients')
disp('in each channel. This norm can be any p-norm, including "inf"-norm.')
disp('Let`s use it to see the results')
disp(' ')
disp('Press any key to continue...')
pause

y = nldifc(im1,7,1,20,linspace(20,100,10),10,1,2,'aos', 'grad', 'dfstep', 5, 'alt1','imscale', 'norm', 'inf');

disp(' ')
disp(' ')
disp('  The 2nd method provided better results then the first. So let`s keep with it a little more.')
disp('  Let`s take as another example the image of a plane and use the 2nd method of diffusion.')
disp(' ')
disp('Press any key to continue...')
pause

im1_path = which('dif_plane.jpg');
im1 = imread(im1_path);

y = nldifc(im1,linspace(1,2,50),linspace(.5,.5,50),5,linspace(200,2000,50),50,1,2,'aos', 'grad', 'dfstep', 5, 'alt1','imscale', 'norm', 'inf');

disp(' ')
disp(' ')
disp('  We can increase even more the simplificatoin, but there is a price to pay. See it')
disp(' ')
disp('Press any key to continue...')
pause

y = nldifc(im1,6,1,8,linspace(200,1000,30),30,1,2,'aos', 'grad', 'dfstep', 5, 'alt1','imscale', 'norm', 'inf');

disp(' ')
disp(' ')
disp('  OK, this simplification did not result in a very usefull image. However it serves to show how')
disp('high contrast areas (the plane shadow) can be preserved while everything else it blurred.')
disp('  Let`s see how we can perform on a human hand.')
disp(' ')
disp('Press any key to continue...')
pause

im1_path = which('dif_hand.jpg');
im1 = imread(im1_path);

y = nldifc(im1,linspace(3,9,30),linspace(3,.5,30),14,1000,30,1,2,'aos', 'grad', 'dfstep', 2, 'alt1','imscale', 'norm', 2);
disp(' ')
disp(' ')
disp('  As a final example, let`s see the diffusion applied to some medical images. Medical image processing')
disp('is a very interesting field. Many automations and improvements can be done in this area, some of them')
disp('can take the diffusion filters as a tool.')
disp(' ')
disp('  The first medical image is a microscopic view of a tissue. Let`s see how we can simplify this image for')
disp('further segmentation and feature extraction.')
disp(' ')
disp('Press any key to continue...')
pause

im1_path = which('dif_tissue.jpg');
im1 = imread(im1_path);

y = nldifc(im1,linspace(3.5,6,30),linspace(2.5,.5,30),10,1000,30,1,2,'aos', 'grad', 'dfstep', 2, 'alt1','imscale', 'norm', 1);

disp(' ')
disp(' ')
disp('  Now let`s see the results on a tomography slice')
disp(' ')
disp('Press any key to continue...')
pause

im1_path = which('dif_tomography.jpg');
im1 = imread(im1_path);

y = nldif(im1,linspace(3,3,30),linspace(.1,.1,30),10,1000,30,1,2,'aos', 'grad', 'dfstep', 2,'imscale');


