% Diffusion filtering toolbox.
% Version 1.1    Mar-2004
% Frederico D'Almeida - DEE/Federal University of Bahia - Brazil
% Comments, improvements, bugs, etc. -> frederico.dalmeida@terra.com.br
%
% Diffusion filters.
%   cedif       - Coherence enhancing diffusion.
%   ldif        - Linear isotropic diffusion.
%   ldifc       - Color Linear isotropic diffusion.
%   nldif1      - 1D Nonlinear isotropic diffusion.
%   nldif       - 2D Nonlinear isotropic diffusion.
%   nldif3      - 3D Nonlinear isotropic diffusion.
%   nldifc      - Color nonlinear isotropic diffusion.
%   pmdif       - Perona & Malick difusion.
%
% Support functions.
%   anidifstep  - Anisotropic diffusion step calculation.
%   aosiso1     - 1D Aditive operator splitting isotropic interation.
%   aosiso      - Aditive operator splitting isotropic interation.
%   conv2br     - 2D convolution with border repetition.
%   diffusivity - Plots the difusivity and flux functions.
%   difplot     - Plots the diffused image.
%   difplot1    - Plots the diffused 1D signal.
%   fluxplot    - Plots the diffusivity and flux in a diffusing image.
%   fluxplot1   - Plots the 1D diffusivity and flux in a diffusing image.
%   isodifstep1 - 1D Isotropic diffusion step calculation.
%   isodifstep  - Isotropic diffusion step calculation.
%   isodifstep3 - 3D Isotropic diffusion step calculation.  
%   orideriv    - Optimized rotation invariance image gradient.
%   gsderiv     - Gaussian smoothed derivate.
%   gsdplot     - Plots gaussian smoothed derivates.
%   stplot      - Plots the structure tensor.
%   thomas      - Fast tridiagonal linear system solver algorithm.
%
% Demonstrations.
%   cedifdemo   - Coherence enhancing diffusion demo.
%   diffdemo    - Diffusivity demonstration.
%   nldif1demo  - 1D Nonlinear diffusion demo.
%   nldifdemo1  - Nonlinear diffusion demo : Basic concepts.
%   nldifdemo2  - Nonlinear diffusion demo : Noise reduction and simplification.
%   nldifdemo3  - Nonlinear diffusion demo : Color images simplification.
%
% Images.
%   dif_hand    - Human hand and big blue ring.
%   dif_house   - House, grass and trees.
%   dif_im1     - Synthetic image (cross, triangle, etc.).
%   dif_im2     - Synthetic shark image.
%   dif_plane   - Airplane image.
%   dif_tissue  - Microscopic tissue image.
%   dif_tomography - Tomography slice image.
%
% Extra functions 
% (These functions are part of the General Extral Toolbox and are supplied
%  because they are used in the demonstrations)
%   colormapc   - Creates circular colormaps.
%   conv2br     - 2D convolution with border repetition.
%   grow        - Expands a matrix by border repetition.
%   noise       - Adds noise to an image/matrix.
%   roll        - Rolls or rotates matrix elements.
%   scale       - Scales matrix elements to a new range.
%