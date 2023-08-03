function PeriodInPixels = period_calib(FileNameRotating,FileNameStatic,Param)
close all;

PeriodCalibFile=Param.PeriodCalibFile;
PeriodEst=Param.PeriodEst;
PixelNx=Param.PixelNx;

PeriodEstRange=0.2; % a parameter specifying the range to search for period
FinePhaseN=30;      % specify how fine to estimate the phase for each pixel

Rotating=imstackread(FileNameRotating+".tif");
Static=imstackread(FileNameStatic+".tif");

Confocal=squeeze(sum(sum(Rotating,1),3));
SinglePhase=squeeze(sum(sum(Static,1),3));
Wave=SinglePhase-Confocal;
Wave=Wave-min(Wave(:));
Wave=Wave/max(Wave(:));
Wave=Wave*2-1;
figure(1);plot(Confocal);hold on;plot(SinglePhase,'r');title('Confocal & SinglePhase');

FWave=abs(fft(Wave));
PeriodIndex=round(PixelNx/PeriodEst*[1-PeriodEstRange 1+PeriodEstRange]);
[Tmp TmpIndex]=max(FWave(PeriodIndex(1):PeriodIndex(2)));
TmpIndex=TmpIndex+PeriodIndex(1)-1;
f=[-5:1:5];
PeriodIndexAvg=sum(f.*FWave(TmpIndex-5:TmpIndex+5))/sum(FWave(TmpIndex-5:TmpIndex+5))+TmpIndex;
Period=PixelNx/(PeriodIndexAvg-1);
figure(2);plot(Wave);title('Wave');
figure(3);plot(FWave);title('FWave');


%% fit phase for each spot by using its neighbouring CorrN spots
CorrN=floor(Period);
PhaseTmp=[0:1:CorrN-1]*2*pi/Period;
for ii=1:FinePhaseN
    Ref(ii,:)=cos(PhaseTmp+ii*2*pi/FinePhaseN);
end
WaveTmp=[];
for ii=1:CorrN
    WaveTmp=[WaveTmp; circshift(Wave,[0 ii-1])];
end
Corr=WaveTmp'*Ref';
[Tmp TmpIndex]=max(Corr,[],2);
Phase=TmpIndex*2*pi/FinePhaseN;
Phase=unwrap(Phase');
if Phase(end)<Phase(1)
    Phase=-Phase;
end
x=[-PixelNx/2:1:PixelNx/2-1];
FitParam=polyfit(x(CorrN:end-CorrN),Phase(CorrN:end-CorrN),3); % remove the first and last CorrN elements for fitting
PeriodInPixels=2*pi/FitParam(3);
PhaseFit=polyval(FitParam,x);
FitParam(4)=0;
FitParam(3)=0;
PhaseOffset=polyval(FitParam,x);
figure(4);plot(PhaseOffset);title('PhaseOffset');
figure(5);plot(Wave,'r');hold on;plot(cos(PhaseFit),'b');
figure(6);plot(Phase,'r');hold on;plot(PhaseFit,'b');

save(PeriodCalibFile,'PeriodInPixels','PhaseOffset');
end