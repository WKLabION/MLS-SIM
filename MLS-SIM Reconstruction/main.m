clear;
close all;
addpath("./subfunctions/");

%% params
Param.PSFScanDirection=[-1 1];                      % specifying the scanning direction in psf measurement with respect to the coordinate of the camera. The positive direction is defined following the convention of a scanning microscope. Dim-1 is vertical, dim-2 is horizontal.
Param.ImgScanDirection=-1;                          % specifying the y scanning direction in imaging mode with respect to the coordinate of the camera. The positive direction is defined following the convention of a scanning microscope
Param.RI=1.33;                                      % refractive of working medium
Param.LineN = 8;                                    % the vertical pixel number used on camera 
CamPixSize = 6.5;                                   % microns
Param.PixelNx = 768;                                % FOV in PSF calculation, which should be matched with experimental image size, this has to be even

Param.Y_shift = 0;
Param.PhaseShift = 0;

Param.PixelFineN=2;                                 % specify the camera pixel reduce factor for the reconstructed high res image: 2 means half pixel size
Param.PhaseFineN=100;                               % determine how fine to estimate phases in phase calibration
Param.PhaseAvgN=40;                                 % specify how to smooth the phase measure conversion curve
Param.PhaseOptimization=1;                          % whether to estimate phase based on collected images

Magnification = 62.857;                             % 62.85 70.72
Param.NA=1.15;                                      % NA of the imaging system
Param.LamdExt=0.488;                                % excitation laser wavelength
Param.LamdDet=0.520;                                % detection fluorescence central laser wavelength

Param.PeriodEst=6;                                  % estimated period in uints of pixels
Param.Bkg=100;                                      % camera bkg pixel reading
Param.PixelNy=640;                                  % FOV in PSF calculation, which should be matched with experimental image size

Param.ResolutionLimit = [80, 80] * 1e-3;            % resolution limit to reduce noise

Param.PSFXYPixelSize=CamPixSize/Magnification/2;    % x-y scanning step size in PSF measurement
Param.XPixelSize=CamPixSize/Magnification;          % equivalent physical size of camera pixel
Param.YPixelSize=CamPixSize/Magnification/10;       % physical scanning step size in imaging plane

%%
param_dir = "./params";
Param.PeriodCalibFile = fullfile(param_dir, 'period_calib');
Param.PSFCalibFile = fullfile(param_dir, 'PSF2');
Param.PSFCalibFile_Y = fullfile(param_dir, 'PSF1');
Param.PhaseCalibFile = fullfile(param_dir, 'phase_calib');
Param.PhaseN=6;                                     % number of phases for SIM illumination


%% reconstruction
Param.Bkg=120;
Param.PSFZ = [1, 2, 3];                             % the z depth used in imaging reconstruction
Param.PSFZ_Y = [1, 2, 3];
Param.CFLinesX = [1, 2, 3, 4, 5];                   % the vertical pixel on camera used in reconstruction
Param.CFLinesY = [1, 2, 3, 4, 5];
Param.X_first = false;                              % the order of two excitation patterns in raw image
intensity_calib = 1.68;                             % the ratio of the peak intensity of two excitation patterns
Param.ItN = [10, 20, 30, 40];                       % the iteration number of the output reconstructed images

fname = './example/1';
reconstruction_basic_Y_memory_saved(fname, Param, intensity_calib); 


