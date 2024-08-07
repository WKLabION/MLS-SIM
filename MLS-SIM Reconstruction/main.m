clear;
close all;
addpath("./subfunctions/");

%% params
Param.PSFScanDirection=[-1 1];    % specifying the scanning direction in psf measurement with respect to the coordinate of the camera. The positive direction is defined following the convention of a scanning microscope. Dim-1 is vertical, dim-2 is horizontal.
Param.ImgScanDirection=-1;        % specifying the y scanning direction in imaging mode with respect to the coordinate of the camera. The positive direction is defined following the convention of a scanning microscope
Param.RI=1.33;              % refractive of working medium
Param.CameraLineReadTime=4.867647;          % the line read speed of the camera in unit of us
Param.LineN = 8;              % this is the numbe of lines on camera 
CamPixSize = 6.5; % microns
Param.PixelNx = 768;          % FOV in PSF calculation, which should be matched with experimental image size, this has to be even

Param.PixelFineN=2;         % specify the camera pixel reduce factor for the reconstructed high res image: 2 means half pixel size
Param.PhaseFineN=100;       % determine how fine to estimate phases in phase calibration
Param.PhaseAvgN=60;         % specify how to smooth the phase measure conversion curve
Param.PhaseOptimization=1;              % whether to estimate phase based on collected images

Magnification = 62.857;  % 62.85 70.72
Param.NA=1.15;                 % NA of the imaging system
Param.LamdExt=0.515;        % excitation laser wavelength
Param.LamdDet=0.538;         % detection fluorescence central laser wavelength

Param.PeriodEst=5;          % estimated period in uints of pixels
Param.Bkg=100;               % camera bkg pixel reading
Param.PixelNy=640;          % FOV in PSF calculation, which should be matched with experimental image size

Param.ResolutionLimit = [80, 80] * 1e-3; % x, y in microns

Param.PSFXYPixelSize=CamPixSize/Magnification/2;  % x-y scanning step size in PSF measurement
Param.XPixelSize=CamPixSize/Magnification;    % equivalent physical size of camera pixel
Param.YPixelSize=CamPixSize/Magnification/10;  % physical scanning step size in imaging plane

Param.PeriodCalibFile = './params/period_calib_G';
Param.PhaseCalibFile = './params/phase_calib_G';
Param.PhaseN=6;             % number of phases for SIM illumination

%% re-center and saturate base psf
Param.PSFPhaseN = 12;
Param.P_shift = (4-4.5) * Param.PixelFineN; % the vertical shift of illumination pattern with regard to the center of camera
Param.S_shift = (3-4.5) * Param.PixelFineN; % the vertical shift of illumination pattern with regard to the center of camera
Param.PSFCalibFile = './params/psfGX';
Param.PSFCalibFile_Y = './params/psfGY';
[PSFTot, PSFTot_Y] = psf_recenterNsaturate(Param);
Param.PSFTot = PSFTot;
Param.PSFTot_Y = PSFTot_Y;

%% reconstruction
Param.Bkg=100;
Param.PSFZ = [1, 2, 3];             % the z depth used in imaging reconstruction
Param.PSFZ_Y = [1, 2, 3];           % the z depth used in imaging reconstruction
Param.CFLinesX = [2, 3, 4, 5, 6];   % the vertical pixel on camera used in reconstruction
Param.CFLinesY = [1, 2, 3, 4, 5];   % the vertical pixel on camera used in reconstruction
Param.X_first = false;              % the order of two excitation patterns in raw image
Param.PhaseOptimization=true;       % specify whether to adjust the global phase shift according to the acquired image
Param.PhaseShift = 0;               % the global phase shift of SIM patterns
intensity_calib = 0.715;            % the ratio of the peak intensity of two excitation patterns
Param.ItN = 40;                     % the iteration number of the output reconstructed images

reconstruction_basic_Y_memory_saved("./example/36", 1, [], Param, intensity_calib); 
