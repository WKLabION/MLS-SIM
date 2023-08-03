function [PSFDet, PSFExt] = psf_preview(Param,PSFMeasureFile)

PSFScanDirection=Param.PSFScanDirection;    % specifying the scanning direction in psf measurement with respect to the coordinate of the camera. The positive direction is defined following the convention of a scanning microscope
Bkg=Param.Bkg;                              % camera bkg pixel reading
XPixelSize=Param.XPixelSize;                % equivalent physical size of camera pixel
PixelFineN=Param.PixelFineN;
PSFXYPixelSize=Param.PSFXYPixelSize;        % x-y scanning step size in PSF measurement
LamdExt=Param.LamdExt;                      % excitation laser wavelength
RI=Param.RI;                                % refractive of working medium

Img=imstackread(PSFMeasureFile);
RealCamH = 64;
RealCamW = 128;
CutCamH = 64;
CutCamW = 64;
ScanNx = 64;
ScanNy = 64;

ScanNxCut=14;            % due to the limited speed of scanning mirror, first several pixels in each row are not sampled at right positions and need to be discarded 

Img=max(0,double(Img)-Bkg);
Img=reshape(Img,[RealCamH,size(Img,1)/RealCamH,size(Img,2),size(Img,3)]);
Img=permute(Img,[1 3 2 4]);
Img=reshape(Img,size(Img,1),size(Img,2),[],ScanNx,ScanNy,size(Img,4));
Img = permute(Img, [1, 2, 4, 5, 3, 6]);
Img = reshape(Img, size(Img, 1), size(Img,2), ScanNx, ScanNy, []);
PSFDetRaw=permute(Img,[1 2 5 4 3]); %[camera_y camera_x scan_z scan_y scan_x]
PSFDetRaw=PSFDetRaw(RealCamH/2+(1-CutCamH/2:CutCamH/2), RealCamW/2+(1-CutCamW/2:CutCamW/2), :, :, ScanNxCut+1:ScanNx);

%% average PSFDet
x=[-CutCamW/2:1:CutCamW/2-1]*XPixelSize;
y=[-CutCamH/2:1:CutCamH/2-1]*XPixelSize;
[x, y]=meshgrid(x,y);
dkx=2*pi/CutCamW/XPixelSize;
dky=2*pi/CutCamH/XPixelSize;
kx=x*dkx/XPixelSize;
ky=y*dky/XPixelSize;
kz=sqrt((2*pi/LamdExt*RI).^2-kx.^2-ky.^2);
OTFDetRaw=fftshift(fftshift(fft2(ifftshift(ifftshift(PSFDetRaw,1),2)),1),2);
OTFDetAvg=0;
for ii=1:size(PSFDetRaw,4)
    for jj=1:size(PSFDetRaw,5)
        PhaseTerm=repmat(exp(-1i*(...
                ky*(ii-size(PSFDetRaw,4)/2)*PSFScanDirection(1)+...
                kx*(jj-size(PSFDetRaw,5)/2)*PSFScanDirection(2))...
            *PSFXYPixelSize),[1 1 size(PSFDetRaw,3)]);
        OTFDetAvg=OTFDetAvg+OTFDetRaw(:,:,:,ii,jj).*PhaseTerm;
    end
end
OTFDetAvg=padarray(OTFDetAvg,(PixelFineN-1)*[size(OTFDetAvg,1)/2 size(OTFDetAvg,2)/2 0],0,'both');
PSFDet=max(0,real(fftshift(fftshift(ifft2(ifftshift(ifftshift(OTFDetAvg,1),2)),1),2)));

f = figure(201); title('raw psf detection');
for ind_z = 1:size(PSFDet,3)
    subplot(3, ceil(size(PSFDet,3)/3), ind_z);
    imagesc(PSFDet(:,:,ind_z)); axis image;
end
linkaxes(f.Children, 'xy');

%% Show PSFExt
PSFExt = PSFDetRaw;
PSFExt = squeeze(sum(PSFExt, [1, 2]));
PSFExt = permute(PSFExt, [2, 3, 1]);

% transform measured PSFExt into Camera coordinates
if PSFScanDirection(1)==1
    PSFExt=flipud(PSFExt);
end
if PSFScanDirection(2)==1
    PSFExt=fliplr(PSFExt);
end

f = figure(202); title('raw psf excitation');
for ind_z = 1:size(PSFExt,3)
    subplot(3, ceil(size(PSFExt,3)/3), ind_z);
    imagesc(PSFExt(:,:,ind_z)); 
    axis image;
end
linkaxes(f.Children, 'xy');

end
