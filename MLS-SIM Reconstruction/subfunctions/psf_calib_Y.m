function [PSFDet, PSFExt] = psf_calib_S(Param, PSFDet, PSFExt, NLfactor)
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
PSFCalibFile_Y = Param.PSFCalibFile_Y;
PSFZ = Param.PSFZ_Y;
PSFZCenter = Param.PSFZCenter;
ExtPhaseRtvZ = Param.ExtPhaseRtvZ;
DetPhaseRtvZ = Param.DetPhaseRtvZ;
Y_shift = Param.S_shift;
dip_zero = Param.dip_zero;
CameraSize = size(PSFExt, 1);

ItN=200;                % Iteration number for phase retrieval
Nxy=400;                % resulting PSF size in calculation

z = ((1:size(PSFDet, 3)) - PSFZ(1)) * PSFZPixelSize;

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
[~, tmpindex]=max(PSFDet(:));
[yn, xn, ZCenter]=ind2sub(size(PSFDet),tmpindex);
tmp=PSFDet(yn-ROI:yn+ROI,xn-ROI:xn+ROI,ZCenter);
tmp=tmp-min(tmp(:));
[tmpx, tmpy]=meshgrid([-ROI:ROI],[-ROI:ROI]);
YCenter=sum(tmp(:).*tmpy(:))/sum(tmp(:))+yn;
XCenter=sum(tmp(:).*tmpx(:))/sum(tmp(:))+xn;
for ii=1:size(PSFDet,3)
    Phase=exp(1i*(ky*(YCenter-Nxy/2-1)*(XPixelSize/PixelFineN)+kx*(XCenter-Nxy/2-1)*(XPixelSize/PixelFineN)+kz*z(ii)));    
    PSFDet(:,:,ii)=(abs(fftshift(ifft2(ifftshift(PupilDetAvg.*Phase))))).^2;
end

figure(3);imagesc(PSFDet(:,:,ZCenter));axis image;


%% phase retrieval of PSFExt
Nxy = 400;
PSFExt=padarray(PSFExt,[(Nxy-size(PSFExt,1))/2 (Nxy-size(PSFExt,2))/2 0],0,'both');
PSFExt = sum(PSFExt, 2);
PSFExt = permute(PSFExt, [1, 3, 2]);
PSFExt = PSFExt ./ sum(PSFExt, 1);

yn=[-Nxy/2:1:Nxy/2-1].';
dky=2*pi/Nxy/PSFXYPixelSize;
ky=yn*dky;
kz=real(sqrt((2*pi/LamdExt*RI).^2-ky.^2));
PupilExtMask=(ky.^2)<((2*pi/LamdExt*NA)^2);

figure(4);subplot(2,3,1);plot(PSFExt(:,PSFZ(1)));
figure(4);subplot(2,3,2);plot(PSFExt(:,PSFZ(2)));
figure(4);subplot(2,3,3);plot(PSFExt(:,PSFZ(3)));           

AmpExp=sqrt(PSFExt);
PupilExtAvg=PupilExtMask;
FieldEst=complex(zeros(Nxy, numel(z)));
FieldExp=complex(zeros(Nxy, numel(z)));
for ii=1:ItN
    for jj=1:size(PSFExt,2)
        FieldEst(:,jj)=fftshift(ifft(ifftshift(PupilExtAvg.*exp(j*kz*z(jj)), 1)), 1);
    end
    for jj=ExtPhaseRtvZ
        FieldExp(:,jj)=AmpExp(:,jj).*exp(j*angle(FieldEst(:,jj)));
        PupilExtEst(:,jj)=fftshift(fft(ifftshift(FieldExp(:,jj), 1)), 1).*exp(-j*kz*z(jj)).*PupilExtMask;
%          PupilExtEst(:,:,jj)=exp(j*angle(fftshift(fft2(ifftshift(FieldExp(:,:,jj)))).*exp(-j*kz*z(jj)))).*PupilExtMask;
    end
    PupilExtAvg=sum(PupilExtEst,2)/size(PSFExt,2);
    figure(4);subplot(2,3,4);plot((abs(FieldEst(:,PSFZ(1))).^2));
    figure(4);subplot(2,3,5);plot((abs(FieldEst(:,PSFZ(2))).^2));
    figure(4);subplot(2,3,6);plot((abs(FieldEst(:,PSFZ(3))).^2));
    figure(5);subplot(1,2,1);plot(abs(PupilExtAvg));
    figure(5);subplot(1,2,2);plot(angle(PupilExtAvg));
    drawnow;
%     pause;
    pause(0.01);
    disp(num2str(ii));
end
PupilExt = PupilExtAvg;
PupilZero = complex(zeros(Nxy));
PupilZero(:, Nxy/2+1) = PupilExt;
PupilExt = PupilZero;

%% expansion and z recentering of PSFExt
PupilExt=padarray(PupilExt,(PixelFineN/(XPixelSize/PSFXYPixelSize)-1)*[Nxy/2, Nxy/2], 0, 'both');
FieldExt = complex(zeros(size(PupilExt)));
Nxy = Nxy * PixelFineN / (XPixelSize/PSFXYPixelSize);
xn=[-Nxy/2:1:Nxy/2-1];
[xn, yn]=meshgrid(xn,xn);
dkxy=2*pi/Nxy/(XPixelSize/PixelFineN);
kx=xn*dkxy;
ky=yn*dkxy;
kz=sqrt((2*pi/LamdExt*RI).^2-kx.^2-ky.^2);
zc = (PSFZCenter - PSFZ(1)) * PSFZPixelSize;
z = z + zc;
for jj = 1:numel(z)
    FieldExt(:,:,jj)=fftshift(ifft2(ifftshift(PupilExt.*exp(1j*kz*z(jj)))));
end
PSFExt=(abs(FieldExt)).^2;

%% center the PSFExt
PSFExtROI=PSFExt(Nxy/2+1-CameraSize:Nxy/2+1+CameraSize,Nxy/2+1-CameraSize:Nxy/2+1+CameraSize,ZCenter);
PSFExtLOI = sum(PSFExtROI, 2);
peak = PSFExtLOI .* islocalmax(PSFExtLOI);
[~, peakidx]=max(peak);
YCenter=peakidx+(Nxy/2-CameraSize)-Y_shift;
XCenter=(Nxy/2+1);

for ii=1:size(PSFExt,3)
    Phase=exp(1j*(ky*(YCenter-Nxy/2-1)*(XPixelSize/PixelFineN)+kx*(XCenter-Nxy/2-1)*(XPixelSize/PixelFineN)+kz*z(ii)));    
    PSFExt(:,:,ii)=(abs(fftshift(ifft2(ifftshift(PupilExt.*Phase))))).^2;
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
Nxy = 400;
PSFDet = PSFDet(size(PSFDet, 1)/2+1+(-Nxy/2:1:Nxy/2-1), size(PSFDet, 1)/2+1+(-Nxy/2:1:Nxy/2-1), :);
PSFExt = PSFExt(size(PSFExt, 1)/2+1+(-Nxy/2:1:Nxy/2-1), size(PSFExt, 1)/2+1+(-Nxy/2:1:Nxy/2-1), :);
x1=[-size(PSFDet,1)/2:1:size(PSFDet,1)/2-1]*XPixelSize/PixelFineN;
[x1, y1]=meshgrid(x1,x1);
x2=[-size(PSFExt,1)/2:1:size(PSFExt,1)/2-1]*XPixelSize/PixelFineN;
[x2, y2]=meshgrid(x2,x2);
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
disp(['Z center is #', num2str(ZCenter)]);
save(char(PSFCalibFile_Y),'PSFDet','PSFExt','PSFTot_Y', '-v7.3');
end