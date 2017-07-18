
% NLDIFDEMO2   Nonlinear diffusion demonstration 2
%
%  Type NLDIFDEMO2.
%
clc

disp('_____________________________Non Linear diffusion demonstration  2________________________')
disp(' ')
disp('  The previous demo (NLDIFDEMO1) showed the basics about the nonlinear diffusion. Now it`s time to')
disp('go a step further and see the real power of this technique.')
disp('  Observe the following synthetic image')
disp(' ')
disp('Press any key to continue...')
pause

im1_path = which('dif_im1.jpg');
im1 = imread(im1_path);
figure(1)
clf
image(im1)
colormap(gray(256))
colorbar
title('Original Image')
xlabel('No noise')

disp(' ')
disp(' ')
disp('  Lets put some noise on it. Starting whith an additive gaussian noise stand. dev. = 25% of image amplitude.')
disp('Watch the result.')
disp('  Note that the image scale has changed as now there are points outside the [0,255] interval; that`s')
disp('why the white area has become gray.')
disp(' ')
disp('Press any key to continue...')
pause

figure(2)
clf
colormap(gray(256))
im1n = noise(im1,'ag','25%');
imagesc(im1n);
title('Noisy image')
xlabel('25% gaussian noise')
colorbar

disp(' ')
disp(' ')
disp('  Let`s run the nonlinear diffusion filter to see what we can do to improve this image')
disp(' ')
disp('Press any key to continue...')
pause

y = nldif(im1n,[linspace(3,15,40) linspace(15,15,10)],[linspace(4,1,40) linspace(1,1,10)],12,10,50,2,4,'dfstep',4,'aos','imscale');

disp(' ')
disp(' ')
disp('  That was nice, wasn`t it? Well, the only problem is that the circle in the lower part of the cross')
disp('has disapeared. That happend because the contrast between that circle and the cross itself was too low.')
disp('This could be avoided by putting lambda (contrast parameter) dow. But then we would have problems to')
disp('eliminate the noise. So, this is one limitation of this technique. However sometimes this is a good')
disp('feature, when you want to perform not a denoising but a image simplification, we`ll see it later.')
disp(' ')
disp('  OK. Now let do it harder. How about gaussian noise std. dev = 50% of amplitude?')
disp(' ')
disp('Press any key to continue...')
pause

y = nldif(noise(im1,'ag','50%'),[linspace(3,15,40) linspace(15,15,10)],[linspace(4,1,40) linspace(1,1,10)],12,10,50,2,4,'dfstep',4,'aos','imscale');

disp(' ')
disp(' ')
disp('  The final test. Gaussian noise std. dev = 100% image amplitude.')
disp(' ')
disp('Press any key to continue...')
pause

y = nldif(noise(im1,'ag','100%'),[linspace(3,15,40) linspace(15,15,10)],[linspace(4,1,40) linspace(1,1,10)],12,10,50,2,4,'dfstep',4,'aos','imscale');

disp(' ')
disp(' ')
disp('  Now let`s exam the image simplification capacity of the nonlinear diffusion. As a first example')
disp('let`s take a sinthetic shark image. Say we want to simplify this image so that we get rid of any small detail.')
disp('That`s easy.')
disp(' ')
disp('Press any key to continue...')
pause

clear im1
im2_path = which('dif_im2.jpg');
im2 = imread(im2_path);
y = nldif(im2,4,1,12,linspace(10,100,14),14,2,1,'aos','grad', 'dfstep', 3);

disp(' ')
disp(' ')
disp('  Despite the eye and teeth are small (size) they have a too big contrast to the rest of the image, so')
disp('they are not eliminated. However, any small contrast structure is removed while the main shape in kept intact.')
disp('We can achieve different simplifications (segmentations) by properly adjusting the parameters.')
disp('Watch it.')
disp(' ')
disp('Press any key to continue...')
pause

y = nldif(im2,[linspace(2,3,10) linspace(3.5,3.75,40)],[linspace(1,1,10) linspace(1,1,40)],20,linspace(100,1000,50),50,2,3,'aos','grad', 'dfstep',2);

disp(' ')
disp(' ')

