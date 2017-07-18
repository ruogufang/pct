
% CEDIFDEMO   Coherence enhancing diffusion demonstration
%
%  Type CEDIFDEMO
%
clc

disp('_____________________________Coherence enhancing diffusion demonstration________________________')
disp(' ')
disp('  Coherence enhancing diffusion is an advanced type of non linear diffusion where the problem')
disp('of not removing noise at image borders is solved. This solution consists in making the diffusion')
disp('at the borders perpendicular to the image gradient (and not parallel to it as usual) or following')
disp('the borders. This kind of diffusion is able to eliminate noise at borders without loosing them.')
disp('  In order to make this kind of diffusion possible, it necessary to modify a little the diffusion')
disp('equation, because we are going to work with two different diffusivities: one parallel to the gradient')
disp('and a new one perpendicular to it.')
disp(' ')
disp('          dH/dt = div( D * grad(H) )  ,  where : H is the concentration, ')
disp('                                                 D is the diffusion tensor ')
disp(' ')
disp('         D = transpose(R) x /d1   0\ x R    ,  where : R is the rotation matrix with')
disp('                            \ 0  d2/                      the gradient orientation,')
disp('                                                       d1 is the diffusivity parallel to the grad,')
disp('                                                       d2 is the diffusivity perp. to the grad.')
disp('         R = ______ 1 ______ * /Gx  -Gy\')
disp('             sqrt(Gx^2+Gy^2)   \Gy   Gx/')
disp(' ')
disp('  The existence of a diffusivity perpendicular to the gradient may sound weird at first as the multiplying')
disp('of perpendicular vectors results zero. However we must remember that, as the orientation of the image is ')
disp('calculated based on a blurred (rho) version of the gradient of a blurred (sigma) version of the image, this')
disp('orientation and the real image gradient ( grad(y) ) will almost never have the same direction.')
disp('  The blurring (sigma) of the image has already been explained in the nhdidemo files. It is necessary to ')
disp('supress small image details that would cause large gradient magnitudes in noisy regions. The ideia of the')
disp('blurring (rho) of the gradient to calculate its orientation follows the same principle. We do not want to')
disp('have the orientation affected be small variations in the gradient. Instead, we want to extract the main flow')
disp('directions of the image.')
disp(' ')
disp('  Before using the Coherence Enhancing diffusion, let`s take a look at the orientations and the effect of the rho')
disp('parameter. Look at the following image')
disp(' ')
disp('Press any key to continue...')
pause

clf
fig = 1;
figure(fig)
colormap(gray(256))
b = zeros(100,100);
b(:,1:75) = 255;
b(30:70,50:100) = 0;
image(b);

disp(' ')
disp(' ')
disp('  Let`s calculate the structure orientation for some different rho values starting with rho = 0. We will use')
disp('sigma = 1 for the following examples. We will add a very small amount of noise to the image to better see the')
disp('effects of rho, otherwise the orientation would be the same for almost the entire image.')
disp(' ')
disp('Press any key to continue...')
pause

y = cedif(noise(b,'ag','1%'),4,1,0,10,0,1,1,1,'struc', 'grad');

disp(' ')
disp(' ')
disp('  Using very small values (zero) for rho cause the image structure orientation to be too sensitive too noise.')
disp('We added gaussian noise with std.dev. = 0.1% of the image amplitude (invisible) and the structure get extremelly')
disp('noisy. Let`s gradually increase the value of rho to see its effect on the orientation.')
disp(' ')
disp('  Using rho = 1.')
disp(' ')
disp('Press any key to continue...')
pause

y = cedif(noise(b,'ag','1%'),4,1,1,10,0,1,1,1,'struc', 'grad');

disp(' ')
disp(' ')
disp('  We can see that the orientation is not changing so fast now. Let`s increase rho even more.')
disp(' ')
disp('  Using rho = 3.')
disp(' ')
disp('Press any key to continue...')
pause

y = cedif(noise(b,'ag','1%'),4,1,3,10,0,1,1,1,'struc', 'grad');

disp(' ')
disp(' ')
disp('  These examples have shown that increasing rho makes the structure orientation to be less sensitive to small')
disp('changes in the gradient direction. However, increasing rho also makes the correct orientations to appear as ')
disp('a larger line.')
disp(' ')
disp('  Let`s now use the coherence enhancing diffusion and see its features. Let`s put some more noise in this image')
disp('to make it a better example. Using rho = 3 and gaussian noise  std.dev. = 25% amplitude')
disp(' ')
disp('Press any key to continue...')
pause

y = cedif(noise(b,'ag','25%'),5,1,3,10,.24,150,1,10,'struc', 'grad', 'dfstep',10);

disp(' ')
disp(' ')
disp('  We can observe that, differently from what happend with the nonlinear diffusion, the noise on the')
disp('borders is quickly eliminated. However, as diffusion is not inhibited on borders, a rounding effect occours.')
disp(' ')
disp('  Let`s use another image to see the results. This image has already been used on the nldifdemo2. We added')
disp('gaussian noise std.dev. 50% amplitude to it.')
disp(' ')
disp('Press any key to continue...')
pause

im1_path = which('dif_im1.jpg');
im1 = imread(im1_path);
y = cedif(noise(im1,'ag','50%'),1.5,4,1,10,.24,50,1,1,'dfstep',25,'imscale', 'grad', 'struc');

disp(' ')
disp(' ')
disp('  Once we can see the quick noise removing on borders and the rounding effect. This rounding causes small')
disp('structures - like the small inner cross - to be too blurred.')
disp('  The observation of the properties of the coherence Enhancing diffusion makes us think about using it togheter')
disp('with the nonlinear diffusion to achieve better results. The first would be used just a little to remove')
disp('noise on borders and the latter would take the job from this point to the end.')
disp(' ')
