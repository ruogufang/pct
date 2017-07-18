
% DIFFDEMO     Diffusivity demonstration
%
%  Type DIFFDEMO
%
clc

disp('____________________________________Diffusivity demonstration________________________________')
disp(' ')
disp('  The diffusivity is the parameter that relates the concentration gradient and the flux. This is')
disp('a general concept and can be applied to heat - where the diffusivity would be the thermal conductance,')
disp('the concentration gradient in the thermic gradient and the flux is the heat flux -, to electric charges')
disp('- where the diffusivity is the electric conductivity, the concentration gradient is the electric potential')
disp('difference and the flux is the electric current - and to many other physical examples.')
disp('  The main difference in the nonlinear diffusion filtering techniques is that this diffusivity is')
disp('made to be a function of the concentration gradient. In physical problems generally this diffusivity is')
disp('either constant or constant for objects made of the same material and dimensions. The effect of variating')
disp('the diffusivity with the concentration gradient - the diffusivity gets smaller as the gradient rises - can')
disp('be visually verified in the NLDIFDEMOs demonstration files. This demonstration is intended to show the')
disp('behavior of the diffusivity for different gradients, lambda and m.')
disp('  The expression of the diffusivity is given by:')
disp(' ')
disp('                         d(grad) = 1 - exp( -Cm / (grad/lambda)^m  )')
disp(' ')
disp('  The constant Cm is automatically calculated for each value m to make the flux ascendig for x<lambda')
disp('and descending for x>=lambda, so the only free variables are grad, lambda and m. As grad is the image')
disp('gradient, only to parameters are to be set: lambda and m.')
disp(' ')
disp('  Look at the following graphic')
disp(' ')
disp('Press any key to continue...')
pause

figure(1)
clf
diffusivity(10,2)

disp('  This graphic represents the diffusivity calculated for different gradient values (x-label), lambda=10 and m=2.')
disp('Note that the diffusivity is monotonically decreasing as the gradient increases. The flux (diffusivity * gradient)')
disp('however increases for grad<lambda and decreases for grad>=lambda. This property is achieved by calculating the Cm')
disp('constant. This constant is found by solving the following expression:')
disp(' ')
disp('                       Cm = root(  1 - exp(-x) - m*x*exp(-x)  )')
disp(' ')
disp('  Note that the expression in parenthesis is the derivate of the flux with respect to the gradient when the ')
disp('gradient is equal to lambda.')
disp(' ')
disp('                  F = grad*d(grad) = grad - grad*exp(-Cm/(grad/lambda)^m)')
disp(' ')
disp('       dF/dgrad = 1 - exp(-Cm/(grad/lambda)^m) - m*Cm(lambda/grad)^m*exp(-Cm/(grad/lambda)^m)')
disp(' ')
disp('  if grad = lambda')
disp(' ')
disp('                             dF/dgrad = 1 - exp(-Cm) - m*Cm*exp(-Cm)')
disp(' ')
disp(' ')
disp(' Let`s see what would happen if we used lambda = 20 and again m = 2.')
disp(' ')
disp('Press any key to continue...')
pause

figure(1)
clf
diffusivity([10,20],2)

disp('  It`s clear to see the effect of lambda in the diffusivity and flux calculations.')
disp('  Now let`s exam the influence of the parameter m. This parameter determines how fast the diffusivity')
disp('changes (allways in the neighborhood of lambda) and also how fast the flux changes. Look at this new plot')
disp('of diffusivities and fluxes for lambda = 20 and m = 2 and 4.')
disp(' ')
disp('Press any key to continue...')
pause

figure(1)
clf
diffusivity(20,[2,4])

disp('  Just to make it clearer the role of m let`s do it again. lambda is fixed 20 and m = 2, 4 and 20.')
disp(' ')
disp('Press any key to continue...')
pause

figure(1)
clf
diffusivity(20,[2,4,20])