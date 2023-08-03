clear all;
close all;

%% set parameters
ZN=20;
Param.MaxShift=[20 20 10];
Param.StepSize=2;
Param.StripeWF=3;

%% set file loading parameters
FilePrefix=['E:\figure 1\'];
FileName=['raw'];
FileSuffix=['.tif'];

Template=imstackread([FilePrefix 'Template.tif']);


%% loop over the time point to do row-by-row correction
gpuDevice();

for StackIdx=1:50
    ImStack=imstackread([FilePrefix FileName FileSuffix],[(StackIdx-1)*ZN+1:StackIdx*ZN]);
    ImStackReg=volume_reg(Template,ImStack);  % rough registration by 3D xcorr 
    for ii=1:size(ImStackReg,3)
        [StackIdx ii]
        RegParamTmp=frame_reg(ImStackReg(:,:,ii),Template,Param); % row by row registration in a frame
        RegParam(ii)=RegParamTmp;
    end

save([FilePrefix num2str(StackIdx+StackIdxOffset) '_reg_13_lines.mat'],'RegParam'); %save the image registration parameters found by xcorr
imstackwrite([FilePrefix num2str(StackIdx+StackIdxOffset) '.tif'],uint16(ImStackReg));
end


%% rough registration by 3D xcorr
function ImStackReg=volume_reg(Template,ImStack)
PadSize=floor((size(Template)-size(ImStack))/2);
ImStackMask=ImStack*0+1;
TemplateMask=(Template>0);

ImStack=padarray(ImStack,PadSize,0,'pre');
ImStack=padarray(ImStack,size(Template)-size(ImStack),0,'post');
ImStackMask=padarray(ImStackMask,PadSize,0,'pre');
ImStackMask=padarray(ImStackMask,size(Template)-size(ImStackMask),0,'post');

Center=floor(size(Template)/2+1);

XCorr=fftshift(abs(ifftn(fftn(Template).*conj(fftn(single(ImStack))))));
Ref=fftshift(sqrt(abs(ifftn(fftn(TemplateMask).*conj(fftn(single(ImStack).^2)))).*abs(ifftn(fftn(Template.^2).*conj(fftn(ImStackMask))))));
XCorr=XCorr./(Ref+eps);
[XCorrMax,Idx]=max(XCorr(:));
[yn xn zn]=ind2sub(size(XCorr),Idx);
Shifts=Center-[yn xn zn]
ImStackReg=circshift(ImStack,-Shifts);
end

%% row-by-row registration in a frame
function RegParam=frame_reg(Frame,Template,Param)
StepSize=Param.StepSize;
StripeWF=Param.StripeWF;
MaxShift=Param.MaxShift;

if sum(Frame(:)>0)
    StripeW=StripeWF*StepSize;

    Target=padarray(double(Frame),[0 MaxShift(2)],0,'both');
    TemplateExt=padarray(double(Template),[MaxShift(1) MaxShift(2) 0],0,'both');
 
    Shifts=[];
    XCorrMax=[];

    YIdx=[1:StepSize:size(Target,1)];

    Target=gpuArray(single(Target));
    TemplateExt=gpuArray(single(TemplateExt));

    for ii=1:length(YIdx)
        ymin=max(1,YIdx(ii)-StripeW);
        ymax=min(size(Target,1),YIdx(ii)+StripeW);
        TargetTmp=padarray(Target(ymin:ymax,:),[MaxShift(1) 0],0,'both');
        TemplateTmp=TemplateExt(ymin:ymax+MaxShift(1)*2,:,:);

        XCorr=abs(ifft2(fft2(TargetTmp).*conj(fft2(TemplateTmp))));
        Ref=sqrt(abs(ifft2(fft2(TargetTmp>0).*conj(fft2(TemplateTmp.^2)))).*abs(ifft2(fft2(TargetTmp.^2).*conj(fft2(TemplateTmp>0)))));
        XCorr=XCorr.*(Ref>(max(Ref(:))/1000));

        XCorr=XCorr./(Ref+eps);
        XCorr=fftshift(fftshift(XCorr,1),2);

        Center=round((size(XCorr,[1 2])+1)/2);
        XCorr=XCorr(Center(1)-MaxShift(1):Center(1)+MaxShift(1),Center(2)-MaxShift(2):Center(2)+MaxShift(2),:);

        [Tmp,Idx]=max(XCorr(:));
        [yn xn zn]=ind2sub(size(XCorr),Idx);
        Shifts(ii,:)=[yn xn zn]-[MaxShift(1)+1 MaxShift(2)+1 0];
        XCorrMax(ii)=Tmp;
    end
    figure(2);subplot(2,1,1);plot(Shifts);
    figure(2);subplot(2,1,2);plot(XCorrMax);
    pause(0.1);

    RegParam.Shifts=Shifts;
    RegParam.XCorrMax=XCorrMax;
    RegParam.Param=Param;
    RegParam.YIdx=YIdx;
else
    RegParam.Shifts=[];
    RegParam.XCorrMax=0;
    RegParam.Param=Param;
    RegParam.YIdx=[];
end


end
