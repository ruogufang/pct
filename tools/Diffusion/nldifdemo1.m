
% NLDIFDEMO1   Non linear diffusion demonstration 1
%
%  Type NLDIFDEMO1.
%
clc

disp('_____________________________Non Linear diffusion demonstration  1________________________')
disp(' ')
disp('  Non linear diffusion, also known as nonlinear diffusion, is a very powerfull image processing')
disp('technique which can be used to image denoising and image simplification.')
disp('  The non linear diffusion is based on an analogy of physical diffusion processes, like the')
disp('temperature diffusion on a metal bar, or the diffusion between two fluids put together.')
disp('These physical diffusion processes are modeled by the following differtial equation:')
disp(' ')
disp('          dH/dt = div( d * grad(H) )  ,  where : H is the concentration (temperature)')
disp('                                                 d is the diffusivity (termal conductance)')
disp(' ')
disp('That means that assuming the diffusivity is constant, which is true for the hot metal bar example and')
disp('for many others, the concentration variation will be faster where the concentration gradient is higher.')
disp('This phenomena is the linear diffusion - see ldifdemo for more details.')
disp('  This kind of behavior is not very interesting for some image processing tasks, as sometimes there is')
disp('the need to preserve image borders (which are, by definition, high gradient areas). Using an linear')
disp('diffusion would quickly destroy the borders.')
disp('  To solve this problem, the nonlinear diffusion makes the diffusivity parameter (d) no longer a')
disp('constant value, but insted the diffusivity becomes a function of the concentration gradient which')
disp('decreases for high gradients')
disp(' ')
disp('          dH/dt = div( d(grad(H)) * grad(H) )')
disp(' ')
disp('  This new formulation allows to perform a image denoising while preserving the borders.')
disp(' ')
disp('  To better understand take a look at the following linear and the non linear diffusion examples.')
disp(' ')
disp('Press any key to continue...')
pause

clf
%fig = gcf;
fig = 1;
figure(fig)
clf
colormap(gray(256))
b = zeros(1,100);
b(1:75) = 255;
image(b);
title('Temperature profile on metal bar')
xlabel('hot = light     cold = dark')

disp(' ')
disp(' ')
disp('  This image represents a metal bar. The light side at the left is colder then the dark part at the right.')
disp('As the time passes, the temperature on the bar tends to equalize: this is an example of linear diffusion.')
disp('Watch the linear diffusion process applied to this image.')
disp(' ')
disp('Press any key to continue...')
pause

y = ldif(b,linspace(1,10,150),150,fig,20);


disp(' ')
disp(' ')
disp('  In oposition to what happened in the linear diffusion, if we apply a nonlinear diffusion to this')
disp('same image the result is that nothing will happen! Why? Because as the diffusivity is very low (even zero) on')
disp('high contrast areas, the diffusion (heat propagation) between these areas is inhibited. Watch the nonlinear')
disp('diffusion in progress.')
disp(' ')
disp('Press any key to continue...')
pause

b = zeros(2,100);
b(:,1:75) = 255;
y = nldif(b,4,0.1,16,.2,50,fig,20,'dfstep',500);

disp(' ')
disp(' ')
disp('  So, what is the good thing about this non linear diffusion if it does not change the image at all? The')
disp('point is, the non linear diffusion does not affect image borders - high gradient (contrast) areas -, but ')
disp('it does affect other parts of the image where the gradient is not so high. Take this example : imagine a noisy')
disp('version of the metal bar image as seen in the picture.')
disp(' ')
disp('Press any key to continue...')
pause

b = zeros(50,100);
b(:,1:75) = 255;
bn = noise(b,'ag','50%');
clf
colormap(gray(256))
image(bn)
title('Temperature profile on metal bar - 50% gaussian noise')
xlabel('hot = light     cold = dark')

disp(' ')
disp(' ')
disp('  If we want to know which part of the bar is hot and which is cold it would be much better if we could remove')
disp('the noise and recover the previous image where the two distinct areas appear. As previoully seen, the linear')
disp('diffusion can remove this noise but it will also blur the entire image (or let the heat go from one side to the')
disp('other). What would happen if we use the non linear diffusion? We expect that the low gradient areas (noise) -')
disp('note that the gradient is taken on a blured version of the original image so that the random noise does not affect')
disp('very much this gradient measure - have a strong diffusion (because their diffusivity is high) while the high gradient')
disp('areas (the temperature interface) have a much weaker diffusion. Take a look at the results')
disp(' ')
disp('Press any key to continue...')
pause

y = nldif(bn,15,2.5,16,2,20,fig,3, 'aos');

disp(' ')
disp(' ')
disp('  Now its clear the power of the nonlinear diffusion filters - they can be used to denoise an image while')
disp('preserving its borders !!!')
disp('  If you take a close look at the image right at the temperature interface you will see that the noise on it')
disp('has not been removed. This happens because, as previouly said, the interface is a high gradient area, and the')
disp('diffusivity is near zero on high gradient areas. So there is almost no diffusion (or no diffusion at all) on ')
disp('the interface, therefore the noise cannot be removed from it.')
disp('  To solve this problem you can use one of both approaches: 1) Make the sigma parameter decrease and the lambda increase')
disp('as the diffusion approaches its end or 2) use the edge enhance diffusion.')
disp('  Decreasing the sigma parameter makes the gradient to be calculated at a less blured version of the image. As a')
disp('consequence, the temperature interface effect gets more and more concentrated to the exact points of the interface. The')
disp('increase on the lambda parameter is needded as while using a less blurred image the noise effect becomes more sensible,')
disp('however, increasing lambda may make the border to be lost (lambda controls which gradient intensities will be diffused')
disp('and which will not), so this strategy is a little complicated and a fine tune on the parameter modification must be done.)')
disp('To help to see this effect look at the image gradient and diffusivity as the diffusion goes on. After the 20th step the')
disp('sigma parameter begins to decrease and lambda begins to increase. Note the effect on the calculated image gradient and ')
disp('on the diffusivity. ')
disp(' ')
disp('Press any key to continue...')
pause

y = nldif(bn,15,2.5,12,2,1,fig,1,'grad', 'flux', 'aos');
figure(fig+1)
colormap(hot(256))
figure(fig+2)
colormap(hot(256))

disp(' ')
disp(' ')
disp(' Adjust the figures to get a nice view and press any key...')
pause


y = nldif(bn,[linspace(10,10,20) linspace(20,50,20)],[linspace(2.5,2.5,20) linspace(2,.01,20)],12,2,40,fig,2,'grad', 'flux', 'aos');