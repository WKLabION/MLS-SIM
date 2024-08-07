function [PSFDet, PSFExt] = psf_calib(Param, PSFDet, PSFExt, NLfactor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in psf measurement, the camera FOV should be square, the scan FOV should also be square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XPixelSize=Param.XPixelSize;                % equivalent physical size of camera pixel
PixelFineN=Param.PixelFineN;
YPixelSize=Param.YPixelSize;
PSFXYPixelSize=Param.PSFXYPixelSize;        % x-y scanning step size in PSF measurement
PSFZPixelSize=Param.PSFZPixelSize;          % z scanning step size in PSF measurement
NA=Param.NA;                                % NA of the imaging system
LamdExt=Param.LamdExt;                      % excitation laser wavelength
LamdDet=Param.LamdDet;                      % excitation laser wavelength
RI=Param.RI;                                % refractive of working medium
LineN=Param.LineN;
PeriodCalibFile=Param.PeriodCalibFile;
PSFPhaseN = Param.PSFPhaseN;
PSFCalibFile=Param.PSFCalibFile;
PSFZ = Param.PSFZ;
ExtPhaseRtvZ = Param.ExtPhaseRtvZ;
DetPhaseRtvZ = Param.DetPhaseRtvZ;
PSFScanDirection = Param.PSFScanDirection;
Y_shift = Param.P_shift;
CameraSize = size(PSFExt, 1);

load(PeriodCalibFile, 'PeriodInPixels');

ItN=200;                % Iteration number for phase retrieval
Nxy=400;                % resulting PSF size in calculation

z = ((1:size(PSFDet, 3)) - PSFZ(1)) * PSFZPixelSize;

%% psf center calibration with maximum intensity
% flatten = @(x) x(:);
% centers = zeros(2, numel(z));
% r = 10;
% for indZ = 1:numel(z)
%     [~, I] = max(flatten(PSFDet(:, :, indZ)));
%     [Y, X] =  ind2sub([size(PSFDet, 1), size(PSFDet, 2)], I);
%     [deltaY, deltaX] = meshgrid(-r:r, -r:r);
%     roi = PSFDet(Y+(-r:r), X+(-r:r), indZ);
%     deltaY = sum(flatten(deltaY .* roi)) ./ sum(roi(:));
%     deltaX = sum(flatten(deltaX .* roi)) ./ sum(roi(:));
%     centers(:, indZ) = [Y+deltaY, X+deltaX];
% end
% centers = centers - [floor(size(PSFDet, 1)/2)+1; floor(size(PSFDet, 2)/2)+1];

%% move psfdet to center
% H = size(PSFDet, 1); W = size(PSFDet, 2);
% x = (-W/2:1:W/2-1)*PSFXYPixelSize;
% y = (-H/2:1:H/2-1)*PSFXYPixelSize;
% [x, y] = meshgrid(x,y);
% dkx = 2 * pi / W / PSFXYPixelSize;
% dky = 2 * pi / H / PSFXYPixelSize;
% kx=x*dkx/PSFXYPixelSize;
% ky=y*dky/PSFXYPixelSize;
% 
% OTFDet=fftshift(fftshift(fft2(ifftshift(ifftshift(PSFDet,1),2)),1),2);
% for indZ = 1:numel(z)
%         PhaseTerm = exp(1i*(ky*centers(1, indZ)+kx*centers(2, indZ))*PSFXYPixelSize);
%         OTFDet(:, :, indZ) = OTFDet(:, :, indZ) .* PhaseTerm;
% end
% PSFDet=max(0,real(fftshift(fftshift(ifft2(ifftshift(ifftshift(OTFDet,1),2)),1),2)));
% 
%% PSFDet phase retrieval
xn=-Nxy/2:1:Nxy/2-1;
[xn, yn]=meshgrid(xn,xn);
dkxy=2*pi/Nxy/(XPixelSize/PixelFineN);
kx=xn*dkxy;
ky=yn*dkxy;
kz=sqrt((2*pi/LamdDet*RI).^2-kx.^2-ky.^2);
PSFDet=padarray(PSFDet,[(Nxy-size(PSFDet,1))/2 (Nxy-size(PSFDet,2))/2 0],0,'both');
PupilDetMask=(kx.^2+ky.^2)<(2*pi/LamdDet*NA)^2;

AmpExp=sqrt(PSFDet);
PupilDetAvg=PupilDetMask;
FieldEst=zeros([size(PupilDetAvg), size(PSFDet, 3)]);
FieldExp=FieldEst;
PupilDetEst = FieldExp;

figure(1);subplot(2,3,1);imagesc(PSFDet(:,:,PSFZ(1)));axis image;
figure(1);subplot(2,3,2);imagesc(PSFDet(:,:,PSFZ(2)));axis image;
figure(1);subplot(2,3,3);imagesc(PSFDet(:,:,PSFZ(3)));axis image;

for ii=1:ItN
    for jj=1:size(PSFDet, 3)
        FieldEst(:,:,jj)=fftshift(ifft2(ifftshift(PupilDetAvg.*exp(1i*kz*z(jj)))));
    end
    for jj=DetPhaseRtvZ
        FieldExp(:,:,jj)=AmpExp(:,:,jj).*exp(1i*angle(FieldEst(:,:,jj)));
        PupilDetEst(:,:,jj)=fftshift(fft2(ifftshift(FieldExp(:,:,jj)))).*exp(-1i*kz*z(jj)).*PupilDetMask;
    end
    PupilDetAvg=sum(PupilDetEst,3)/size(PSFDet,3);
    PupilDetAvg=PupilDetMask.*exp(1i*angle(PupilDetAvg));
    figure(1);subplot(2,3,4);imagesc((abs(FieldEst(:,:,PSFZ(1))).^2));axis image;
    figure(1);subplot(2,3,5);imagesc((abs(FieldEst(:,:,PSFZ(2))).^2));axis image;
    figure(1);subplot(2,3,6);imagesc((abs(FieldEst(:,:,PSFZ(3))).^2));axis image;
    figure(2);subplot(1,2,1);imagesc(abs(PupilDetAvg));axis image;
    figure(2);subplot(1,2,2);imagesc(angle(PupilDetAvg));axis image;
    drawnow;
    disp(num2str(ii));
end
PSFDet=(abs(FieldEst)).^2;

%% center the PSFDet
ROI=1;
[tmp tmpindex]=max(PSFDet(:));
[yn xn ZCenter]=ind2sub(size(PSFDet),tmpindex);
tmp=PSFDet(yn-ROI:yn+ROI,xn-ROI:xn+ROI,ZCenter);
tmp=tmp-min(tmp(:));
[tmpx tmpy]=meshgrid([-ROI:ROI],[-ROI:ROI]);
YCenter=sum(tmp(:).*tmpy(:))/sum(tmp(:))+yn;
XCenter=sum(tmp(:).*tmpx(:))/sum(tmp(:))+xn;
for ii=1:size(PSFDet,3)
    Phase=exp(j*(ky*(YCenter-Nxy/2-1)*PSFXYPixelSize+kx*(XCenter-Nxy/2-1)*PSFXYPixelSize+kz*z(ii)));    
    PSFDet(:,:,ii)=(abs(fftshift(ifft2(ifftshift(PupilDetAvg.*Phase))))).^2;
end
figure(3);imagesc(PSFDet(:,:,ZCenter));axis image;
PupilDet = PupilDetAvg;

%% phase retrieval of PSFExt
PSFExt=padarray(PSFExt,[(Nxy-size(PSFExt,1))/2 (Nxy-size(PSFExt,2))/2 0],0,'both');

xn=[-Nxy/2:1:Nxy/2-1];
[xn yn]=meshgrid(xn,xn);
dkxy=2*pi/Nxy/PSFXYPixelSize;
kx=xn*dkxy;
ky=yn*dkxy;
kz=sqrt((2*pi/LamdExt*RI).^2-kx.^2-ky.^2);
PupilExtMask=(kx.^2+ky.^2)<(2*pi/LamdExt*NA)^2;
PupilExtMask=PupilExtMask.*((abs(xn)==round(Nxy/(PeriodInPixels*XPixelSize/PSFXYPixelSize)))+(xn==0));
        
figure(4);subplot(2,3,1);imagesc(PSFExt(:,:,PSFZ(1)));axis image;
figure(4);subplot(2,3,2);imagesc(PSFExt(:,:,PSFZ(2)));axis image;
figure(4);subplot(2,3,3);imagesc(PSFExt(:,:,PSFZ(3)));axis image;            

AmpExp=sqrt(PSFExt);
PupilExtAvg=PupilExtMask;
FieldEst=[];
FieldExp=[];
for ii=1:ItN
    for jj=1:size(PSFExt,3)
        FieldEst(:,:,jj)=fftshift(ifft2(ifftshift(PupilExtAvg.*exp(j*kz*z(jj)))));
    end
    for jj=ExtPhaseRtvZ
        FieldExp(:,:,jj)=AmpExp(:,:,jj).*exp(j*angle(FieldEst(:,:,jj)));
        PupilExtEst(:,:,jj)=fftshift(fft2(ifftshift(FieldExp(:,:,jj)))).*exp(-j*kz*z(jj)).*PupilExtMask;
%          PupilExtEst(:,:,jj)=exp(j*angle(fftshift(fft2(ifftshift(FieldExp(:,:,jj)))).*exp(-j*kz*z(jj)))).*PupilExtMask;
    end
    PupilExtAvg=sum(PupilExtEst,3)/size(PSFExt,3);
    figure(4);subplot(2,3,4);imagesc((abs(FieldEst(:,:,PSFZ(1))).^2));axis image;
    figure(4);subplot(2,3,5);imagesc((abs(FieldEst(:,:,PSFZ(2))).^2));axis image;
    figure(4);subplot(2,3,6);imagesc((abs(FieldEst(:,:,PSFZ(3))).^2));axis image;
    figure(5);subplot(1,2,1);imagesc(abs(PupilExtAvg));axis image;
    figure(5);subplot(1,2,2);imagesc(angle(PupilExtAvg));axis image;
    drawnow;
%     pause;
    pause(0.01);
    disp(num2str(ii));
end
PSFExt=(abs(FieldEst)).^2;
PupilExt = PupilExtAvg;

%% center the PSFExt
ROI=1;
PSFExtROI=PSFExt(Nxy/2+1-CameraSize:Nxy/2+1+CameraSize,Nxy/2+1-CameraSize:Nxy/2+1+CameraSize,ZCenter);
[tmp tmpindex]=max(PSFExtROI(:));
[yn xn]=ind2sub(size(PSFExtROI),tmpindex);
yn=yn+(Nxy/2-CameraSize);
xn=xn+(Nxy/2-CameraSize);
tmp=PSFExt(yn-ROI:yn+ROI,xn-ROI:xn+ROI,ZCenter);
tmp=tmp-min(tmp(:));
[tmpx tmpy]=meshgrid([-ROI:ROI],[-ROI:ROI]);
YCenter=sum(tmp(:).*tmpy(:))/sum(tmp(:))+yn-Y_shift;
XCenter=sum(tmp(:).*tmpx(:))/sum(tmp(:))+xn;

for ii=1:size(PSFExt,3)
    Phase=exp(j*(ky*(YCenter-Nxy/2-1)*PSFXYPixelSize+kx*(XCenter-Nxy/2-1)*PSFXYPixelSize+kz*z(ii)));    
    PSFExt(:,:,ii)=(abs(fftshift(ifft2(ifftshift(PupilExtAvg.*Phase))))).^2;
end
figure(6);imagesc(PSFExt(:,:,ZCenter));axis image;

%% saturate PSFExt
if nargin == 4
    PSFExt = satu(PSFExt, NLfactor);
    figure(6);imagesc(PSFExt(:,:,ZCenter));axis image;
end

%% normalize PSF
PSFExt = PSFExt ./ max(PSFExt(:));
PSFDet = PSFDet ./ max(PSFDet(:));

%% prepare different psf
x1=[-size(PSFDet,1)/2:1:size(PSFDet,1)/2-1]*XPixelSize/PixelFineN;
[x1, y1]=meshgrid(x1,x1);
x2=[-size(PSFExt,1)/2:1:size(PSFExt,1)/2-1]*PSFXYPixelSize;
[x2, y2]=meshgrid(x2,x2);
x3=[-size(PSFDet,1)/2:1:size(PSFDet,1)/2-1]*XPixelSize/PixelFineN;
y3=[-size(PSFDet,1)/2:1:size(PSFDet,1)/2-1]*YPixelSize;
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
disp(['Z center is #', num2str(ZCenter)]);
save(PSFCalibFile,'PSFDet','PSFExt','PSFTot', '-v7.3');
end