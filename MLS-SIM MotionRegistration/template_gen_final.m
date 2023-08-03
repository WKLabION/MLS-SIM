clear all;
close all;


%% set parameters
ZN=20;                  % number of z planes in a stack
EdgeSize=10;                    
MaxShift=[15 15 5];
XCorrThd=0.7;
FltCoef=[0.1 0.1 1];

%% set file loading parameters
FilePrefix=['E:\figure 1\'];
FileName=['raw_1'];
FileSuffix=['.tif'];
StackIdxList=[1:3:50];


%% load selected images
ImStack=[];
for ii=1:length(StackIdxList)
    StackIdx=StackIdxList(ii);
    ImStack(:,:,:,ii)=imstackread([FilePrefix FileName FileSuffix],[(StackIdx-1)*ZN+1:StackIdx*ZN]);
end

YN=size(ImStack,1);
XN=size(ImStack,2);
ZN=size(ImStack,3);

%% eliminate edges of the images
Mask=padarray(ones(YN-2*EdgeSize,XN-2*EdgeSize),[EdgeSize EdgeSize],0,'both');
ImStack=ImStack.*Mask;

%% generate first template simply by averaging
Template1=mean(single(ImStack),4);

%% pad zeros to increase the size of the template so that larger shift can be corrected
Mask1=(Template1*0+1).*Mask;
Template1=padarray(Template1,MaxShift,0,'both');
Mask1=padarray(Mask1,MaxShift,0,'both');
ImStack1=padarray(ImStack,[MaxShift 0],0,'both');


%% using 3D xcorr to find out the shift of the image stack and generate second template by correcting these shifts.
FMask1=fftn(Mask1);

XCorrCoef=[];
Shifts=[];
Center=floor(size(Template1)/2+1);
Template2=0;
Template2Mask=0;
for ii=1:size(ImStack1,4)
    XCorr=fftshift(abs(ifftn(fftn(Template1).*conj(fftn(single(ImStack1(:,:,:,ii)))))));
    Ref=fftshift(sqrt(abs(ifftn(FMask1.*conj(fftn(single(ImStack1(:,:,:,ii)).^2)))).*abs(ifftn(fftn(Template1.^2).*conj(FMask1)))));
    XCorr=XCorr./(Ref+eps);
    [XCorrMax,Idx]=max(XCorr(:));
    XCorrCoef(ii)=XCorrMax;
    [yn xn zn]=ind2sub(size(XCorr),Idx);
    Shifts(ii,:)=Center-[yn xn zn];
    figure(1);imagesc(mip_view(XCorr));
    figure(2);plot(Shifts);
    figure(3);plot(XCorrCoef);

    if XCorrMax>XCorrThd
        Template2=Template2+circshift(ImStack1(:,:,:,ii),-Shifts(ii,:));
        Template2Mask=Template2Mask+circshift(Mask1,-Shifts(ii,:));
    end
    drawnow;
end

Template2=Template2./(Template2Mask+0.01);
imstackwrite([FilePrefix 'Template.tif'],uint16(Template2));