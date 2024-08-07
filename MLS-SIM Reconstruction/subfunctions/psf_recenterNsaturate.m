function [PSFTot, PSFTot_Y] = psf_recenterNsaturate(Param, ~, dilation_factor, zero_bias)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in psf measurement, the camera FOV should be square, the scan FOV should also be square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XPixelSize=Param.XPixelSize;                % equivalent physical size of camera pixel
PixelFineN=Param.PixelFineN;
YPixelSize=Param.YPixelSize;
LineN=Param.LineN;
PeriodCalibFile=Param.PeriodCalibFile;
PSFPhaseN = Param.PSFPhaseN;
PSFCalibFile_Y = Param.PSFCalibFile_Y;
PSFCalibFile=Param.PSFCalibFile;
P_shift = Param.P_shift;
S_shift = Param.S_shift;
if ~isfield(Param, 'saturate_factor') || isempty(Param.saturate_factor)
    saturate_factor = 0.05;
else
    saturate_factor = Param.saturate_factor;
end
if nargin<=2 || isempty(dilation_factor)
    dilation_factor = [1, 1];
end
if nargin<=3
    zero_bias = 0;
end
load(PeriodCalibFile, 'PeriodInPixels');
Nxy=400;                % resulting PSF size in calculation

%% prepare x psf
load(PSFCalibFile, 'PSFExt', 'PSFDet');
PSFExt = satu(PSFExt, saturate_factor);
Nxy = 400;
PSFDet = PSFDet(size(PSFDet, 1)/2+1+(-Nxy/2:1:Nxy/2-1), size(PSFDet, 1)/2+1+(-Nxy/2:1:Nxy/2-1), :);
PSFExt = PSFExt(size(PSFExt, 1)/2+1+(-Nxy/2:1:Nxy/2-1), size(PSFExt, 1)/2+1+(-Nxy/2:1:Nxy/2-1), :);
x1=[-size(PSFDet,1)/2:1:size(PSFDet,1)/2-1]*XPixelSize/PixelFineN;
[x1, y1]=meshgrid(x1,x1);
x2=[-size(PSFExt,1)/2:1:size(PSFExt,1)/2-1]*XPixelSize/PixelFineN;
y2 = ((-size(PSFExt,1)/2:1:size(PSFExt,1)/2-1)*dilation_factor(1)+P_shift)*XPixelSize/PixelFineN;
[x2, y2]=meshgrid(x2,y2);
x3=[-Nxy/2:1:Nxy/2-1]*XPixelSize/PixelFineN;
y3=[-Nxy/2:1:Nxy/2-1]*YPixelSize;
[x3, y3]=meshgrid(x3,y3);

PSFTot = zeros(Nxy, Nxy, size(PSFDet, 3), PSFPhaseN, LineN, 'double');
for ii=1:size(PSFDet,3) % iteration of z plane
    for jj=1:LineN % iteration of camera line
        Tmp1=interp2(x1,y1,PSFDet(:,:,ii),x3,y3+(jj-(1+LineN)/2)*XPixelSize,'linear',0); % magnified PSFDetection
        for kk=1:PSFPhaseN
            Tmp2=interp2(x2,y2,PSFExt(:,:,ii),x3+(kk-1/2)/PSFPhaseN*PeriodInPixels*XPixelSize,y3,'linear',0); % magnified PSFExcitation
            PSFTot(:,:,ii,kk,jj)=Tmp1.*circshift(rot90(Tmp2,2),[1 1]);  %[y x z psfphase lines]
        end
    end
end

%% prepare y psf
load(PSFCalibFile_Y, 'PSFExt', 'PSFDet');
PSFExt = satu(PSFExt, saturate_factor);
Nxy = 400;
PSFDet = PSFDet(size(PSFDet, 1)/2+1+(-Nxy/2:1:Nxy/2-1), size(PSFDet, 1)/2+1+(-Nxy/2:1:Nxy/2-1), :);
PSFExt = PSFExt(size(PSFExt, 1)/2+1+(-Nxy/2:1:Nxy/2-1), size(PSFExt, 1)/2+1+(-Nxy/2:1:Nxy/2-1), :);
PSFExt = PSFExt.*(1-zero_bias)+zero_bias;
y1=(-size(PSFDet,1)/2:1:size(PSFDet,1)/2-1)*XPixelSize/PixelFineN;
x1=(-size(PSFDet,1)/2:1:size(PSFDet,1)/2-1)*XPixelSize/PixelFineN;
[x1, y1]=meshgrid(x1,y1);
x2=[-size(PSFExt,1)/2:1:size(PSFExt,1)/2-1]*XPixelSize/PixelFineN;
y2 = ((-size(PSFExt,1)/2:1:size(PSFExt,1)/2-1)*dilation_factor(2)+S_shift)*XPixelSize/PixelFineN;
[x2, y2]=meshgrid(x2,y2);
x3=[-Nxy/2:1:Nxy/2-1]*XPixelSize/PixelFineN;
y3=[-Nxy/2:1:Nxy/2-1]*YPixelSize;
[x3, y3]=meshgrid(x3,y3);

PSFTot_Y = zeros(Nxy, Nxy, size(PSFDet, 3), 1, LineN, 'double');
for ii=1:size(PSFDet,3) % iteration of z plane
    for jj=1:LineN % iteration of camera line
        Tmp1=interp2(x1,y1,PSFDet(:,:,ii),x3,y3+(jj-(1+LineN)/2)*XPixelSize,'linear',0); % magnified PSFDetection
        Tmp2=interp2(x2,y2,PSFExt(:,:,ii),x3,y3,'linear',0); % magnified PSFExcitation
        PSFTot_Y(:,:,ii,1,jj)=Tmp1.*circshift(rot90(Tmp2,2),[1 1]);  %[y x z psfphase lines]
    end
end

end