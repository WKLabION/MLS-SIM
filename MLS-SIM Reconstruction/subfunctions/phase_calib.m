function phase_calib(FileName ,Param)
close all;

PeriodCalibFile=Param.PeriodCalibFile;
PhaseCalibFile=Param.PhaseCalibFile;

LineN=Param.LineN;
Loi = Param.Loi;
Bkg=Param.Bkg;
PixelNx=Param.PixelNx;          
PhaseN=Param.PhaseN;
PhaseFineN=Param.PhaseFineN;
PhaseAvgN=Param.PhaseAvgN;   % a paramter to smooth the curve

load(PeriodCalibFile);

half_flag = false;

%% load and preprocess image
Img=imstackread(FileName+".tif");
Img=reshape(double(Img),LineN,size(Img,1)/LineN,size(Img,2),size(Img,3));
Img=Img-Bkg;
Img=Img.*(Img>0);
Img=squeeze(sum(Img(Loi,:,:),1));
figure(1);imagesc(Img); title('raw image');

%% load and preprocess phase measurement
PhaseMeasure=dlmread(FileName+".txt");
PhaseMeasure=PhaseMeasure-mean(PhaseMeasure);
PhaseMeasure=PhaseMeasure/mean(abs(PhaseMeasure));
figure(2);plot(PhaseMeasure);title('raw phase measure');

%% half-image process
H = size(Img, 1);
if half_flag
    Img = Img((1:H/2)+H/4, (1:PixelNx/2)+PixelNx/4);
    PhaseOffset = PhaseOffset((1:PixelNx/2)+PixelNx/4); %#ok<*NODEF>
    PixelNx = PixelNx/2;
    PhaseMeasure = PhaseMeasure(numel(PhaseMeasure)/4 + (1:numel(PhaseMeasure)/2))
end
%% phase correlation
x=[-PixelNx/2:1:PixelNx/2-1];
for ii=1:PhaseFineN
    Ref(ii,:)=cos(x/PeriodInPixels*2*pi+PhaseOffset+ii*2*pi/PhaseFineN);
end
Corr=Img*Ref';
[Tmp TmpIndex]=max(Corr,[],2);
Phase=TmpIndex/PhaseFineN*2*pi;
figure(3);plot(Phase);title('Estimated Phase from Img');
figure(4); 
plot(cos(Phase),'b');hold on;
plot(PhaseMeasure,'r');hold off;
title('Estimated Phase VS Measured Phase');

[Phase PhaseIndex]=sort(Phase);
PhaseMeasure=PhaseMeasure(PhaseIndex);
figure(7);plot(Phase,PhaseMeasure,'*'); title('PhaseFit VS PhaseMeasure');


Phase=reshape(Phase,PhaseAvgN,length(Phase)/PhaseAvgN);
PhaseMeasure=reshape(PhaseMeasure,PhaseAvgN,length(PhaseMeasure)/PhaseAvgN);
Phase=sum(Phase,1)/PhaseAvgN;
PhaseMeasure=sum(PhaseMeasure,1)/PhaseAvgN;


%%
figure(8); title('smoothed PhaseFit VS Phase Measure');
plot(Phase,PhaseMeasure,'r*');hold on;
Phase=[Phase-2*pi Phase Phase+2*pi Phase+4*pi Phase+6*pi];
PhaseMeasure=[PhaseMeasure PhaseMeasure PhaseMeasure PhaseMeasure PhaseMeasure];
Tmp=[2*pi/PhaseFineN:2*pi/PhaseFineN:2*pi];
for ii=1:PhaseN
    PhaseMeasureX(ii,:)=Tmp+(ii-1)*2*pi/PhaseN;
    PhaseMeasureY(ii,:)=interp1(Phase,PhaseMeasure,PhaseMeasureX(ii,:),'spline');
end
plot(PhaseMeasureX(1,:),PhaseMeasureY(1,:),'g');hold off;
save(PhaseCalibFile,'PhaseMeasureX','PhaseMeasureY');
end

