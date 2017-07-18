
% NLDIF1DEMO   1D Non linear diffusion demonstration
%
%  Type NLDIF1DEMO.
%
clc
warning off

disp('_____________________________ 1D Non Linear diffusion demonstration ________________________')
disp(' ')
disp('  This demonstration shows the power of the Nonlinear diffusion filters applied to 1D signals.')
disp('We start with a binary digital signal show on figure 1.')
disp(' ')
disp(' ')

figure(1)
clf
x = square(0:0.01:6*pi);
plot(x)
ini_axis = axis;
ini_axis = ini_axis.*[1 1 1.2 1.2];
axis(ini_axis)
title('Original Signal')
disp('Press any key to continue...')
pause

disp(' ')
disp(' ')
disp('  We add to our original signal an additive gaussian noise with variance 10% of the signal.')
disp('amplitude. The noisy signal is show superimposed to the original one.')
disp(' ')
disp(' ')

xn = noise(x,'ag','10%');
plot(xn,'r');
hold
plot(x,'b')
title('Original ans Noisy Signals')
disp(' ')
disp(' ')
disp('Press any key to continue...')
pause

disp(' ')
disp(' ')
disp('  We add now a multiplicative noise, replacing 10% of our measured samples by a random value')
disp('taken from a gaussian distribuition with variance set to the signal amplitude.')
disp(' ')
disp(' ')

xn = noise(xn,'mg','10%');
clf
plot(xn,'r');
hold
plot(x,'b')
title('Original ans Noisy Signals')


disp(' ')
disp(' ')
disp('Press any key to continue...')
pause

disp(' ')
disp(' ')
disp('  This noisy signal will be used as the input to the 1D Nonlinear Diffusion filter. We use')
disp('lambda 0.5, sigma 10 -> 1, stepsize .5 -> 50, 300 steps. (x -> y means the parameter')
disp('changes linearly from x to y as the steps are performed.)')
disp(' ')
disp(' ')

n = 100;
y = nldif1(xn, .019, linspace(50,1,n), 10, linspace(.5,1000,n), n,1,3, 'scale','aos', 'grad');


disp(' ')
disp(' ')
disp('Press any key to continue...')
pause

disp(' ')
disp(' ')
disp('  The next figure shows a comparison between the original, noisy and the restored signals.')
disp('Note that not only the signal shape is almost completelly restored, but also the signal')
disp('transisions are placed exactly at the right places.')
disp(' ')
disp(' ')


clf
subplot(2,1,1)
plot(xn,'y')
hold
plot(x,'b');
plot(y,'r')
%axis(ini_axis)
title('Original (blue), Noisy (yellow) and Restored (red) Signals')
subplot(2,1,2)
plot(xn-x,'y')
hold
plot(y-x,'r')
title('Noisy (yellow) and Restored (red) Error')


disp(' ')
disp(' ')
disp('Press any key to continue...')
pause

disp(' ')
disp(' ')
disp('  Lets try to recover the signal after adding a gaussian noise with variance 100% of the')
disp('signal amplitude.')
disp(' ')
disp(' ')

clf
xn = noise(x,'ag','100%');
plot(xn)
title('Noisy Signal')

disp(' ')
disp(' ')
disp('Press any key to continue...')
pause
n = 200;
y = nldif1(xn, .019, linspace(50,1,n), 10, linspace(.5,1000,n), n,1,3, 'scale','aos', 'grad');

disp(' ')
disp(' ')
disp('Press any key to continue...')
pause

clf
subplot(2,1,1)
plot(xn,'y')
hold
plot(x,'b');
plot(y,'r')
%axis(ini_axis)
title('Original (blue), Noisy (yellow) and Restored (red) Signals')
subplot(2,1,2)
plot(xn-x,'y')
hold
plot(y-x,'r')
title('Noisy (yellow) and Restored (red) Error')

disp(' ')
disp(' ')
disp('Press any key to continue...')
pause

disp(' ')
disp(' ')
disp('  The results again are impressive. The signal is recovered, and the transitions are')
disp('again correctely placed.')
disp(' ')
disp(' ')