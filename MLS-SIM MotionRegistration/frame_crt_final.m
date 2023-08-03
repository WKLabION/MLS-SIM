clear all;
close all;

%% set parameters and load template
Param.MaxShift=[20 20 5];
ZN=20;

FilePrefix=['E:\figure 1\'];
FileSuffix=['.tif'];

Template=imstackread([FilePrefix 'Template.tif']);

%% loop over the time series imaging stacks
for StackIdx=1:50
    load([FilePrefix num2str(StackIdx) '_reg_13_lines.mat']);  %load registration parameters found by xcorr
    ImStack=imstackread([FilePrefix num2str(StackIdx) '.tif']);
    ImStackCrt=ImStack*0;
    for ii=1:size(ImStack,3)
        [StackIdx ii]
        if ~isempty(RegParam(ii).Shifts)
            FrameCrt=frame_crt(ImStack(:,:,ii),RegParam(ii),Template);  %correct the image based on registration parameters.
        else
            FrameCrt=double(ImStack)*0;
        end
        ImStackCrt(:,:,ii)=uint16(max(FrameCrt(:,:,ZRange),[],3));
        pause(0.1);
%         pause;
    end
    imstackwrite([FilePrefix num2str(StackIdx) '_crt.tif'],ImStackCrt);
    imwrite(max(ImStackCrt,[],3),[FilePrefix ['\mip\'] num2str(StackIdx) '_mip.tif']);
end


%% correct images based on registration parameters
function FrameCrt=frame_crt(Frame,RegParam,Template)

%% set parameters
XCorrThd=0.8;
dShiftThd=5;
ShiftMax=10;
KernelSpacing=5;

ImStackSize=size(Template);

figure(100);plot(RegParam.Shifts);hold on;plot(RegParam.XCorrMax);hold off;
%% if the shift found by xcorr is too large to be correct, tease them out
XCorrMax=RegParam.XCorrMax;
YShiftRaw=RegParam.Shifts(:,1);
XShiftRaw=RegParam.Shifts(:,2);
ZShiftRaw=RegParam.Shifts(:,3);
YIdx=RegParam.YIdx;

YShift=YShiftRaw;
XShift=XShiftRaw;
ZShift=ZShiftRaw;

W=max(XCorrMax(:)-XCorrThd,0);

dYShift=YShift-circshift(YShift,-1);
dXShift=XShift-circshift(XShift,-1);
dZShift=ZShift-circshift(ZShift,-1);
dShift=sqrt(dYShift.^2+dXShift.^2+dZShift.^2);
W=W.*(abs(YShift)<=ShiftMax).*(abs(XShift)<=ShiftMax).*(dShift<=dShiftThd);

YShift=min(max(YShift,-ShiftMax),ShiftMax);
XShift=min(max(XShift,-ShiftMax),ShiftMax);
ZShift=min(max(ZShift,1),ImStackSize(3));


%% find continous regions that the shifts found by registration is reliable.
FitGroup=find_fit_group(W);
figure(1000);plot(FitGroup);


YFit=YShift*0;
XFit=YFit;
ZFit=YFit;

%% in each continous regions, curves of shifts are fitted by a collection of gaussion functions centered at different positions to yield smooth curves
FitGroupN=max(FitGroup);
for ii=1:FitGroupN
    TmpIdx=find(FitGroup==ii);
    t=[1:length(TmpIdx)];
    if round(max(t)/KernelSpacing)<=3
        t0=[1 round(max(t)/2) max(t)]';
    else
        dt=KernelSpacing;
        t0=[1:dt:floor(max(t)/dt)*dt]';
        t0=[t0;max(t)];
    end
    dt=t0(2)-t0(1);
    T=[];
    for jj=1:size(t0,1)
        T(:,jj)=exp(-(t-t0(jj)).^2/(dt)^2);
    end
    T_inv=inv(T'*diag(W(TmpIdx))*T)*T'*diag(W(TmpIdx));
    YCoef=T_inv*YShift(TmpIdx);
    XCoef=T_inv*XShift(TmpIdx);
    ZCoef=T_inv*ZShift(TmpIdx);
    YFit(TmpIdx)=T*YCoef;
    XFit(TmpIdx)=T*XCoef;
    ZFit(TmpIdx)=T*ZCoef;
end




figure(1);plot([YShift XShift ZShift],'-');hold on;plot([YFit XFit ZFit],'--');plot(W*5,'k');hold off;

YFit1=interp1(YIdx,YFit,[YIdx(1):YIdx(end)],'linear');
YFit1=[ones(YIdx(1)-1,1)*YFit(1);YFit1(:);ones(size(Frame,1)-YIdx(end),1)*YFit(end)];
XFit1=interp1(YIdx,XFit,[YIdx(1):YIdx(end)],'linear');
XFit1=[ones(YIdx(1)-1,1)*XFit(1);XFit1(:);ones(size(Frame,1)-YIdx(end),1)*XFit(end)];
ZFit1=interp1(YIdx,ZFit,[YIdx(1):YIdx(end)],'linear');
ZFit1=[ones(YIdx(1)-1,1)*ZFit(1);ZFit1(:);ones(size(Frame,1)-YIdx(end),1)*ZFit(end)];
W1=double(FitGroup>0);
W1=interp1(YIdx,W1,[YIdx(1):YIdx(end)],'nearest');
W1=[ones(YIdx(1)-1,1)*W1(1);W1(:);ones(size(Frame,1)-YIdx(end),1)*W1(end)];
ShiftsFit=[YFit1 XFit1 ZFit1];


%% interpolate the image to correct the distortion caused by shifts
[x0 y0]=meshgrid([1:size(Frame,2)],[1:size(Frame,1)]);
x1=x0(1,:);
y1=f_inv(ShiftsFit(:,1));
[x1 y1]=meshgrid(x1,y1);

Tmp=interp2(x0,y0,double(Frame),x0+ShiftsFit(:,2),y0,'linear',0);
FrameCrtTmp=uint16(interp2(x0,y0,double(Tmp),x1,y1,'linear',0));
ZFit=max(min(ShiftsFit(:,3),ImStackSize(3)),1);
FrameCrt=zeros(ImStackSize);
for ii=1:size(Frame,1)    
    z0=max(1,min(round(ShiftsFit(ii,3)),ImStackSize(3)));
    FrameCrt(ii,:,z0)=FrameCrtTmp(ii,:).*W1(ii);
end
% figure(1000);imagesc(Frame);axis image;
% figure(1001);imagesc(FrameCrtTmp);axis image;
% figure(1002);imagesc(max(FrameCrt,[],3));axis image;
% pause;

end


function y_inv=f_inv(dy)
UpSamplingFactor=5;
x0=[1:length(dy)];
y0=x0-dy(:)';
x1=linspace(x0(1),x0(end),length(x0)*UpSamplingFactor);
y1=interp1(x0,y0,x1);

for ii=1:length(dy)
    [tmp Idx]=min(abs(y1-ii));
    y_inv(ii)=x1(Idx);
end
end

function FitGroup=find_fit_group(Weight)
LongRange=15;
LongRangeThd=12;
ShortRange=3;
ShortRangeThd=3;
MinFitLength=15;

W=Weight(:)>0;

WLong=conv(W,ones(LongRange,1),'same');
WShort=conv(W,ones(ShortRange,1),'same');
WLong(1)=0;WLong(end)=0;
WShort(1)=0;WLong(end)=0;

Tmp=W*0;
Tmp(2:end)=(WLong(2:end)>=LongRangeThd).*(WLong(1:end-1)<LongRangeThd);
WLongStartIdx=find(Tmp);
Tmp=W*0;
Tmp(2:end)=(WLong(1:end-1)>=LongRangeThd).*(WLong(2:end)<LongRangeThd);
WLongEndIdx=find(Tmp);

FitGroup=W*0;
FitGroupIdx=1;
for ii=1:length(WLongStartIdx)
    Tmp=WShort(WLongStartIdx(ii):WLongEndIdx(ii))>=ShortRangeThd;
    ShortValid=find(Tmp);
    if length(ShortValid)>=MinFitLength
        FitGroup(WLongStartIdx(ii)+ShortValid(1)-1:WLongStartIdx(ii)+ShortValid(end)-1)=FitGroupIdx;
        FitGroupIdx=FitGroupIdx+1;
    end
end

end

