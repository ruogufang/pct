README for Tensor Total Variation (TTV) Package

Ruogu Fang
Advanced Multimedia Laboratory
Department of Electrical and Computer Engineering
Cornell University
June 21, 2014


Tensor total variation deconvolution package uses a tensor total variation to robustly estimate the perfusion parameters in low-dose, low sampling rate CT perfusion with delay in tracer arrival time. In this package, we estimate the perfusion parameter cerebral blood flow (CBF), cerebral blood volume (CBV) and mean transit time (MTT), as described in the following papers. 


CITATION:
--------------
Please cite the following papers if you use code in this SPD package.

Ruogu Fang, Pina Sanelli, Shaoting Zhang, Tsuhan Chen.
Tensor Total-Variation Regularized Deconvolution for Efficient Low-Dose CT Perfusion.
MICCAI'14, The 17th Annual International Conference on Medical Image Computing and Computer Assisted Intervention, 2014.


BIBTEX:
-----------------
@incollection{fang2014tensor
year={2014},
booktitle={Medical Image Computing and Computer-Assisted Intervention – MICCAI 2014},
series={Lecture Notes in Computer Science},
editor={Barillot, Golland, Hornegger, Howe},
title={Tensor Total-Variation Regularized Deconvolution for Efficient Low-Dose CT Perfusion},
publisher={Springer Berlin Heidelberg},
author={Fang, Ruogu and Sanelli, Pina C. and Zhang, Shaoting and Chen, Tsuhan },
}


FILES ORGANIZATION:
----------------------------------
Main_TTV.m: main file to run the demo. It tests four algorithms: sSVD, bSVD, Tikhonov and TTV. 
Exp_NPS_RIF: Noise power spectrum and residue impulse functions recovered from low-dose CTP data.
Exp_AIF_RIF_delay: Delayed arterial input function and residue impulse functions restored from low-dose CTP data.
Exp_CBFs: Estimated CBF values at different true CBF values.
Exp_PSNR: Estimated CBF values at different PSNR levels.
Exp_Uniform: Visual display of perfusion parameters (CBF, CBV and MTT) of a uniform region.
Exp_PCNR: Visual display of CBF and MTT estimation at certain PCNR level for a region with two different CBF values.
Exp_Clinical: CBF estimation of a clinical dataset at low-dose, low sampling rate.
Exp_ROI: Enlarged ROI regions of saved clinical maps in data folder.
Exp_Boxplot: Boxplots of RMSE and Lin's CCC in the paper.
Exp_para_reg: Parameter tuning for TV weight \gamma.
Exp_para_ratio: Parameter tunign for the ratio between temporal and spatial TV weight \gamma_t / \gamma_s.
TTV_FCSA: Core TTV algorithm implementation.
denoise_TTV: Solver for tensor total variation.


INSTALLATION: 
--------------------------
1.  Unzip the package. 
2.  Run Main_TTV.m. This file executes TTV algorithm and compares with sSVD, bSVD and Tikhonov algorithms on a clinical CTP dataset.
3.  For figures in the above mentioned papers, run Exp*.m files which indicated the Figure number in its headline.
4.  SPD comparison is included in all Exp*.m files. To execute SPD algorithm, follow the step below:
     (1) Go to toolbox/SPD_matlab v1.0/Utilities/, compile two packages: ompbox10 and spams-matlab. The binaries for Mac OS 10.9 is already included in the package. For different platforms and OS, 
   a. Go to ompbox10 folder, follow the readme.txt file to make files in MATLAB.
   b. Go to spams-matlab folder and compile following the instructions in HOW_TO_INSTALL.txt. 
   * Note that for newer version of Mac OS (10.7+) or gcc (XCode), issues may arise when running compile.m. Please follow the notes at the end of this README file for possible solutions.
   * If you encounter problems in compiling SPD package, you could refer to README.txt in SPD_matlab v1.0 folder. Or you can comment out all SPD related experiments in the Exp*.m files.
5. You may change the tube current-exposure time product (mAs) and sampling rate in Main_TTV.m, Exp_Clinical.m and other Exp*.m files. You may also use your own simulated or real low-dose image by loading different DICOM or MAT files as Vn (the low-dose CTP data [T x X x Y]). 


TOOLBOX INCLUDED:
-------------------------------
This package already includes the utility software packages downloaded from other website. Please properly cite the related papers if you use these packages.
a. PCT v1.0:   Perfusion computed tomography toolbox written by Ruogu Fang, Kolbeinn Karlsson.

b. SPD-matlab v1.0:  Sparse perfusion deconvolution toolbox written by Ruogu Fang. http://chenlab.ece.cornell.edu/people/ruogu/publications/SPD_matlab%20v1.0.zip



Please contact rf294@cornell.edu if you have issues or suggestions for this TTV package.
