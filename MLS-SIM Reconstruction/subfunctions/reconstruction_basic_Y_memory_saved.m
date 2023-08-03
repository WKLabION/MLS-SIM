function ObjRecon=reconstruction_basic_Y_memory_saved(ImgFileName, Param, illumination_calib)
close all;

cfft2 = @(x) fftshift(fftshift(fft2(ifftshift(ifftshift(x, 1), 2)), 1), 2);
icfft2 = @(x) fftshift(fftshift(ifft2(ifftshift(ifftshift(x, 1), 2)), 1), 2);
norm_mean = @(x) x./mean(x(:));
m0r = @(x) max(0, real(x));
fftshift2 = @(x) fftshift(fftshift(x, 1), 2);

%% loading params
PhaseFineN=Param.PhaseFineN;
ImgScanDirection=Param.ImgScanDirection;  
PixelFineN=Param.PixelFineN;
XPixelSize = Param.XPixelSize;
YPixelSize = Param.YPixelSize;
LineN=Param.LineN;
PhaseN=Param.PhaseN;
PixelNx=Param.PixelNx;
CFLinesX = Param.CFLinesX;
CFLinesY = Param.CFLinesY;
PSFZ=Param.PSFZ;
PSFZ_Y = Param.PSFZ_Y;
ItN=Param.ItN;
PhaseOptimization=Param.PhaseOptimization;
ResolutionLimit = Param.ResolutionLimit;
PhaseShift = Param.PhaseShift;

PeriodCalibFile=Param.PeriodCalibFile;
PhaseCalibFile=Param.PhaseCalibFile;
load(PeriodCalibFile);
load(PhaseCalibFile);
if ~isfield(Param, 'PSFTot') || isempty(Param.PSFTot)
    f = load(Param.PSFCalibFile, 'PSFTot');
    PSFTot = f.PSFTot;
else
    PSFTot = Param.PSFTot;
end

if ~isfield(Param, 'PSFTot_Y') || isempty(Param.PSFTot_Y)
    f = load(Param.PSFCalibFile_Y, 'PSFTot_Y');
    PSFTot_Y = f.PSFTot_Y;
else
    PSFTot_Y = Param.PSFTot_Y;
end

if ~isfield(Param, 'PSFXPick') || isempty(Param.PSFXPick)
    PSFXPick = 1:size(PSFTot, 4);
else
    PSFXPick = Param.PSFXPick;
end

if ~isfield(Param, 'saturate_factor') || isempty(Param.saturate_factor)
    saturate_factor = 0.05;
else
    saturate_factor = Param.saturate_factor;
end

if ~isfield(Param, 'PartDiv') || isempty(Param.PartDiv)
    PartDiv = 1;
else
    PartDiv = Param.PartDiv;
    PartN = Param.PartN;
    PartL = Param.PartL;
end
if ~isfield(Param, 'PartDiv_X') || isempty(Param.PartDiv_X)
    PartDiv_X = 1;
else
    PartDiv_X = Param.PartDiv_X;
    PartN_X = Param.PartN_X;
    PartL_X = Param.PartL_X;
end



%% load raw image
Img=imread(strcat(ImgFileName, '.tif'));
Img=reshape(double(Img),LineN,size(Img,1)/LineN,size(Img,2));
PixelNy = size(Img, 2);
if mod(PixelNy, 1920)~=0
    PixelNy = PixelNy - mod(PixelNy, 1920);
    Img = Img(:, 1:PixelNy, :);
    % to remove galvo flyback lines
end
Img=squeeze(Img(:,:,:));

% subtract bkg and pixel gain
if isfield(Param, 'fbkg') && ~isempty(Param.fbkg)
    fbkg = load(Param.fbkg);
    bkg = fbkg.bg;
    fgain = load(Param.fgain);
    gain = fgain.bg;
    gain = gain - bkg;
    gain = gain ./ mean(gain, 'all');
    gain = 1*(gain-1)+1;
    Img = single(Img);
    Img = Img - permute(bkg, [1, 3, 2]);
    Img = Img ./ permute(gain, [1, 3, 2]);
    Img = max(0, Img);
else
    Img = Img - Param.Bkg;
    Img = max(Img, 0);
end

if PartDiv~=1 && ~isempty(PartDiv)
    Img = Img(:, size(Img, 2)*PartN/PartDiv + (1:size(Img, 2)*PartL/PartDiv), :);
    PixelNy = PixelNy*PartL/PartDiv;
end
if PartDiv_X~=1 && ~isempty(PartDiv_X)
    Img = Img(:, :, size(Img, 3)*PartN_X/PartDiv_X + (1:size(Img, 3)*PartL_X/PartDiv_X));
    PixelNx = PixelNx*PartL_X/PartDiv_X;
    PhaseOffset = PhaseOffset(numel(PhaseOffset)*PartN_X/PartDiv_X + (1:numel(PhaseOffset)*PartL_X/PartDiv_X));
end
if mod(PixelNy, PhaseN*2)~=0
    Img = Img(:, 1:PixelNy-mod(PixelNy, PhaseN*2), :);
    PixelNy = PixelNy - mod(PixelNy, PhaseN*2);
end

Img=permute(Img,[2 3 1]); %[y x lines]


if Param.X_first
    xp_start = 1;
    ys_start = 2;
else
    xp_start = 2;
    ys_start = 1;
end
ImgSP = Img;
Img = Img(xp_start:2:end, :, :);


%% load phase measurement
PhaseMeasure=dlmread(strcat(ImgFileName, ".txt"));
PhaseMeasure=squeeze(PhaseMeasure);
if PartDiv~=1 && ~isempty(PartDiv)
    PhaseMeasure = PhaseMeasure(numel(PhaseMeasure)*PartN/PartDiv + (1:numel(PhaseMeasure)*PartL/PartDiv));
end
if mod(numel(PhaseMeasure), PhaseN*2)~=0
    PhaseMeasure = PhaseMeasure(1:numel(PhaseMeasure)-mod(numel(PhaseMeasure), PhaseN*2));
end

PhaseMeasure = PhaseMeasure(:, xp_start:2:end);

%% calc phase
PhaseMeasure=PhaseMeasure-mean(PhaseMeasure);
PhaseMeasure=PhaseMeasure/mean(abs(PhaseMeasure));
PhaseMeasureGroup=reshape(PhaseMeasure,PhaseN,length(PhaseMeasure)/PhaseN);
PhaseMeasureCorr=PhaseMeasureGroup'*PhaseMeasureY;
[~, TmpIndex]=max(PhaseMeasureCorr,[],2);
PhaseMeasureFit=TmpIndex*2*pi/PhaseFineN;

%% optimize phase use experimental data
if PhaseOptimization==1
    PhaseOptimizationNum = 28;
    LineIndex = round(numel(CFLinesX)/2);
    xxx=[-PixelNx/2:1:PixelNx/2-1];
    PhaseSearch=[-pi:pi/PhaseOptimizationNum:pi];
    Corr=PhaseSearch*0;
    for ii=1:length(PhaseSearch)
        for jj=1:size(Img,1)
            Ref=cos(xxx/PeriodInPixels*2*pi+PhaseOffset+PhaseSearch(ii)+PhaseMeasureFit(ceil(jj/PhaseN))+mod((jj-1),PhaseN)*2*pi/PhaseN);
            Ref = satu(Ref, saturate_factor);
            Corr(ii)=Corr(ii)+sum(Ref.*Img(jj,:,LineIndex));
        end
    end
    [~, TmpIndex]=max(Corr);
    PhaseMeasureFit=PhaseMeasureFit+PhaseSearch(TmpIndex);
end


%% calc phase for each pixel
PixelNy = PixelNy/2;
x=[-PixelNx/2:1:PixelNx/2-1];
PhaseMap = x / PeriodInPixels * 2 * pi + PhaseOffset;
PhaseMap = reshape(PhaseMap, [1, numel(PhaseMap)]);
PhaseMap = repmat(PhaseMap, [PixelNy, 1, LineN]);
LineIdx = 1:PixelNy;
temp = PhaseMeasureFit(ceil(LineIdx / PhaseN)) + mod((LineIdx'-1), PhaseN) * (2 * pi / PhaseN);
PhaseMap = PhaseMap + temp(:);
PhaseMap = PhaseMap + PhaseShift;

PhaseMap=mod(PhaseMap,2*pi);
PhaseMapL = zeros([size(PhaseMap, 1)*2, size(PhaseMap, 2), size(PhaseMap, 3)], 'double');
PhaseMapL(xp_start:2:end, :, :) = PhaseMap;
PhaseMapL(ys_start:2:end, :, :) = +Inf;
PhaseMap = PhaseMapL;


Ub = 1:size(PSFTot, 4); % [1, n]
Ub = permute(Ub, [1, 3, 4, 2]); % [1, 1, 1, n]
Ub = repmat(Ub, [size(PhaseMap), 1]); % [dim(PhaseMap), n]
Lb = Ub - 1;
Ub = Ub .* (2*pi/size(PSFTot,4));
Lb = Lb .* (2*pi/size(PSFTot,4));

PhaseMask = (PhaseMap>=Lb) & (PhaseMap<=Ub);
PhaseMask = PhaseMask(:, :, :, PSFXPick);
PhaseMask(:, :, :, end+1) = isinf(PhaseMap)&(PhaseMap>0);
PhaseMask=permute(PhaseMask,[1 2 4 3]); %[y x psfphase lines]

%% rearrange images
Img = ImgSP;
if ImgScanDirection==-1
    Img=flipud(Img);
    PhaseMask=flipud(PhaseMask);
end

Img = permute(Img, [1, 2, 4, 3]); % [y, x, 1, Nlines]
ImgFine = PhaseMask .* Img;

ImgFine=padarray(ImgFine,[0 0 0 0 PixelFineN-1],0,'post');
ImgFine=permute(ImgFine,[1 5 2 3 4]);
ImgFine=reshape(ImgFine,[size(ImgFine,1) size(ImgFine,2)*size(ImgFine,3) size(ImgFine,4) size(ImgFine,5)]);
% ImgFine [y, x, phase, lines]
ImgFineX = ImgFine(:, :, 1:end-1, CFLinesX);
ImgFineY = ImgFine(:, :, end, CFLinesY);

PhaseMask=padarray(PhaseMask,[0 0 0 0 PixelFineN-1],0,'post');
PhaseMask=permute(PhaseMask,[1 5 2 3 4]);
PhaseMask=reshape(PhaseMask,[size(PhaseMask,1) size(PhaseMask,2)*size(PhaseMask,3) size(PhaseMask,4) size(PhaseMask,5)]);
PhaseMask = logical(PhaseMask); % [y, x, phase, lines]
PhaseMask = gpuArray(PhaseMask); % [y, x, phase, lines]

PSFTot=PSFTot(:,:,PSFZ, PSFXPick, CFLinesX);
PSFTot_Y=PSFTot_Y(:,:,PSFZ_Y, :, CFLinesY);
PSFTot_Y = PSFTot_Y * illumination_calib;
PSFTot = gpuArray(single(PSFTot)); % [y, x, z, phase, lines]
PSFTot_Y = gpuArray(single(PSFTot_Y));

NormalizationX = zeros([size(ImgFineX, 1), size(ImgFineX, 2), size(PSFTot, 3)],'single', 'gpuArray');
NormalizationY = zeros([size(ImgFineY, 1), size(ImgFineY, 2), size(PSFTot, 3)],'single', 'gpuArray');
for line_idx = 1:numel(CFLinesX)
    for phase_idx = 1:numel(PSFXPick)
        for z_idx = 1:numel(PSFZ)
            padded_PSF = padarray(squeeze(PSFTot(:, :, z_idx, phase_idx, line_idx)),...
                [(size(ImgFineX,1)-size(PSFTot,1))/2, (size(ImgFineX,2)-size(PSFTot,2))/2, 0, 0], 0, 'both');

            mask_weight = m0r(ifft2(fft2(PhaseMask(:, :, phase_idx, CFLinesX(line_idx))).*conj(fft2(padded_PSF))));
            NormalizationX(:, :, z_idx) = NormalizationX(:, :, z_idx) + mask_weight;
        end
    end
end
for line_idx = 1:numel(CFLinesY)
    for z_idx = 1:numel(PSFZ_Y)
        padded_PSF = padarray(squeeze(PSFTot_Y(:, :, z_idx, :, line_idx)),...
            [(size(ImgFineY,1)-size(PSFTot_Y,1))/2, (size(ImgFineY,2)-size(PSFTot_Y,2))/2, 0, 0], 0, 'both');
        mask_weight = m0r(ifft2(fft2(PhaseMask(:, :, end, CFLinesY(line_idx))).*conj(fft2(padded_PSF))));
        NormalizationY(:, :, z_idx) = NormalizationY(:, :, z_idx) + mask_weight;
    end
end
NormalizationX = fftshift2(NormalizationX);
NormalizationY = fftshift2(NormalizationY);

clearvars fPhaseMask mask_weight;

%% Resolution Limit Filter
H = size(ImgFineX, 1); W = size(ImgFineX, 2);
X = -W/2:(W/2-1); Y = -H/2:(H/2-1);
[X, Y] = meshgrid(X, Y);
dx = XPixelSize / PixelFineN; dy = YPixelSize;
dkx = 2 * pi / dx / W; dky = 2 * pi / dy / H;
kx = dkx .* X; ky = dky .* Y;
kx_limit = 2*pi/ResolutionLimit(1); ky_limit = 2*pi/ResolutionLimit(2);
OTF_mask = ((kx./kx_limit).^2+(ky./ky_limit).^2) <= 1 ;

%% Reconstruction
PhaseMask=gpuArray(PhaseMask);
ImgFineX=gpuArray(single(ImgFineX));
ImgFineY = gpuArray(single(ImgFineY));
ImgFineEst=gpuArray(zeros(size(ImgFineX,1),size(ImgFineX,2),'single'));
ObjRecon=gpuArray(ones(size(ImgFineX,1),size(ImgFineX,2),length(PSFZ),'single'));
ObjReconAvgX=gpuArray(zeros(size(ObjRecon),'single'));
ObjReconAvgY=gpuArray(zeros(size(ObjRecon),'single'));
OTF_mask = gpuArray(single(OTF_mask));

if PartDiv~=1 && ~isempty(PartDiv)
    ImgFileName = [ImgFileName, '_p', num2str(PartN), 'd', num2str(PartDiv)];
end

for ii=1:max(ItN)
    disp(['iteration: ' num2str(ii) ' of ' num2str(max(ItN))]);
    ObjReconAvgX=ObjReconAvgX*0;
    ObjReconAvgY=ObjReconAvgY*0;
    
    for psf_idx = 1:numel(PSFXPick)
        for Lineidx = 1:numel(CFLinesX)
            PSFTemp = padarray(PSFTot(:, :, :, psf_idx, Lineidx),...
                                [(size(ImgFineX,1)-size(PSFTot,1))/2, (size(ImgFineX,2)-size(PSFTot,2))/2, 0], 0, 'both');
            OTFTemp = fft2(PSFTemp);
    
            ImgFineEst = sum(m0r(ifft2(fft2(ObjRecon) .* OTFTemp)), 3);
            ImgFineEst = fftshift(ImgFineEst);
    
            Ratio=ImgFineX(:,:,psf_idx,Lineidx)./ImgFineEst;
            Ratio = Ratio .* PhaseMask(:,:,psf_idx,CFLinesX(Lineidx));
            ObjReconAvgX = ObjReconAvgX + ObjRecon.*fftshift2(m0r(ifft2(fft2(Ratio).*conj(OTFTemp))));              
        end
    end

    for Lineidx = 1:numel(CFLinesY)
        PSFTemp = padarray(PSFTot_Y(:, :, :, :, Lineidx),...
            [(size(ImgFineY,1)-size(PSFTot_Y,1))/2, (size(ImgFineY,2)-size(PSFTot_Y,2))/2, 0], 0, 'both');
        OTFTemp = fft2(PSFTemp);

        ImgFineEst = sum(m0r(ifft2(fft2(ObjRecon) .* OTFTemp)), 3);
        ImgFineEst = fftshift(ImgFineEst);

        Ratio=ImgFineY(:,:,end,Lineidx)./ImgFineEst;

       
        Ratio = Ratio .* PhaseMask(:,:,end,CFLinesY(Lineidx));
        ObjReconAvgY = ObjReconAvgY + ObjRecon.*fftshift2(m0r(ifft2(fft2(Ratio).*conj(OTFTemp))));
    end
    
    ObjRecon = norm_mean(ObjReconAvgX ./ NormalizationX + ObjReconAvgY ./ NormalizationY);
    
    for z_idx = 1:size(ObjRecon, 3)
        ObjRecon(:, :, z_idx) = max(0, real(icfft2(cfft2(ObjRecon(:, :, z_idx)) .* OTF_mask)));
    end

    subplot(121);imagesc(ObjRecon(:,:,1));
    subplot(122);imagesc(ImgFineEst);
    if ismember(ii, ItN)
        ObjReconFinal = gather(ObjRecon)*200;
        imstackwrite(strcat(ImgFileName, '_itr', num2str(ii), '.tif'), ...
            uint16(imbin(ObjReconFinal(:, :, :), 5, 1)));
    end
end
clearvars ImgFineEst ImgFineX ImgFineY ...
    ObjRecon ObjReconAvgX ObjReconAvgY ...
    padded_PSF PSFTemp Ratio PhaseMask;
end

