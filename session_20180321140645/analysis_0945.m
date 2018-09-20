%analysis for MB14 fiber photometry, session 1

close all;clear all;
cd '/n/groups/datta/mel/d1+snc/session_20180314160945';


load('kinect_object.mat');
ext = extract_object;
phot = ext.neural_data.photometry;

animal = ext.mouse_id;


%% plotting the photometry signals

figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(1,2,1)
plot(phot.traces(1).dff, 'green')
title('Signal from dopaminergic terminals');
xlabel('Time (frame)'); ylabel('Intensity');

xlim([0 size(phot.traces(1).dff,1)]);
subplot(1,2,2)
plot(phot.traces(4).dff, 'red')
title('Signal from D1 neurons');
xlabel('Time');
xlim([0 size(phot.traces(1).dff,1)]);

% the raw traces show the bleaching, and the dff ones are already corrected

% looking at ext.projections you find the length of the different fields (ie number of frames) for the kinect data
% doing length (phot.traces(1).dff) you find the number of frames for the tdt data.
% the kinect data can have more frames and have NaNs. To overcome this, we can do

vel_frames = ext.get_original_timebase(ext.projections.velocity_mag);
% for example, to get the proper number of frames for the velocity


%% getting the scalars and their derivatives

%taking the data for the different scalars and putting them in an array
scalar(:,1) = ext.get_original_timebase(ext.projections.angle);
scalar(:,2) = ext.get_original_timebase(ext.projections.width);
scalar(:,3) = ext.get_original_timebase(ext.projections.length);
scalar(:,4) = ext.get_original_timebase(ext.projections.height_ave);
scalar(:,5) = ext.get_original_timebase(ext.projections.velocity_mag);
scalar(:,6) = ext.get_original_timebase(ext.projections.velocity_mag_3d);

%saving the names
scalarNames = {'angle','width','length','height_ave','velocity_mag','velocity_mag_3d'};

%remove the frames with Nans in the photometry from the kinect data (for example, the vel_frames)
scalar = scalar(all(~isnan(phot.traces(4).dff),2),:);
scalars = scalar(all(~isnan(scalar),2),:);

%take the derivatives
Dscalars = diff(scalars);

%zscoring them
Zscalars = zscore(scalars);
ZDscalars = zscore(Dscalars);


%% close-up to the photometry signals

%remove the rows in the photometry data with with nans:
dopamine = phot.traces(1).dff(all(~isnan(phot.traces(1).dff),2),:); 
D1 = phot.traces(4).dff(all(~isnan(phot.traces(4).dff),2),:);

%remove the same frames in the kinect data (for example, the vel_frames)
velocity = vel_frames(all(~isnan(phot.traces(4).dff),2),:);

%re-scale the time axis to seconds
time = [1/30:1/30:length(velocity)/30];


%close-up of the signals
figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(3,1,1)
plot(time,dopamine, 'green')
hold on
refline(0,mean(dopamine))
title('Close-up to signal from dopaminergic terminals');
ylabel('Intensity');
xlim([100 200]);

subplot(3,1,2)
plot(time,D1, 'red')
hold on
refline(0,mean(D1))
title('Close-up to signal from D1 neurons');
ylabel('Intensity');
xlim([100 200]);

subplot(3,1,3)
plot(time,velocity, 'black')
title('velocity');
xlabel('Time (s)'); ylabel('Velocity');
xlim([100 200]);


% 
% g=gramm('x',time,'y',dopamine);
% g.geom_line();
% figure, g.draw()

%% spectral analysis

%shortgreen = phot.traces(1).dff(5000:30000);
%shortred = phot.traces(4).dff(5000:30000);

Fs = 30; %sampling frequency (30 Hz)

FTg = fft(dopamine-mean(dopamine)); %calculate the fast fourier transform
%substract the mean to not get a huge power for a 0Hz frequency
magFTg = abs(FTg);
fx = 0:(length(dopamine)/2)-1; %build up a proper frequency axis
fx = (fx*Fs)/length(dopamine); %scaling it so that it represents the frequencies in Hz

figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(2,1,1)
plot(fx,magFTg(1:length(dopamine)/2),'g')
grid
title('Spectrum of the signal from dopaminergic terminals');
ylabel('Magnitude (dB)');

FTr = fft(D1-mean(D1)); %calculate the fast fourier transform
magFTr = abs(FTr);
subplot(2,1,2)
plot(fx,magFTg(1:length(D1)/2),'r')
grid
title('Spectrum of the signal from D1 neurons');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

%% 

%obtaining the spectrogram
[s,f,t]=spectrogram(dopamine-mean(dopamine),60,55,[],30);
[s2,f2,t2]=spectrogram(D1-mean(D1),60,55,[],30);
% this computes the spectrogram using a short fourier transform
%this is taking 60 frames (2s) as a window, and 55 frames as overlap, of a
%30 hz sampling frequency
%s is the output of the fourier transform
%f are the cyclical frequencies
%t are the time instants

figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(1,2,1), imagesc(t,f,20*log10(abs(s))); %plot the result
axis xy %to invert the y axis
caxis([-20 0]) %this sets the colormap limits
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Spectrogram for the dopamine signal');


subplot(1,2,2), imagesc(t2,f2,20*log10(abs(s2))); %plot the result
axis xy %to invert the y axis
caxis([-20 0])
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Spectrogram for the D1 signal');
c = colorbar
c.Label.String = 'Power/frequency (dB/Hz)';

%% 

figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(2,1,1),
plot(f,20*log10(mean(abs(s),2)),'g')
xlabel('Frequency (Hz)'); ylabel('Power/frequency (dB/Hz)');
hold on
plot(f2,20*log10(mean(abs(s2),2)),'r')

%take the derivative of the curves to see if they have similar
%characteristics

a=diff(20*log10(mean(abs(s),2)));
b=diff(20*log10(mean(abs(s2),2)));
dist=norm(a-b);

subplot(2,1,2),
plot(f(2:length(f)),a,'g')
xlabel('Frequency (Hz)'); ylabel('derivative of Power/frequency (dB/Hz)');
hold on
plot(f2(2:length(f2)),b,'r')
title('comparing the shape of the curves')

%% cross-correlation

Zdopamine = zscore(dopamine);
ZD1 = zscore(D1);

[XCF,lags,bounds] = crosscorr(Zdopamine,ZD1,220,3); %shows the corr values for the different signal lags.
lagtime = lags/30;

[Rgr, lags2] = xcorr(ZD1,Zdopamine,220,'coeff'); 
lagtime2 = lags2/30;
%the correlation is computed inversely that for crosscorr regarding which
%signal is considered first

%find the max crosscorr and xcorr values
[val,ind] = max(XCF);
[val2,ind2] = max(Rgr);

figure,  set(gcf, 'Position', [800, 1000, 1000, 800])
subplot(2,1,1), plot(lagtime,XCF, 'ro')
hold on
plot(lagtime,XCF, 'r')
ylim([min(XCF)-0.05,max(XCF)+0.05]); xlim([-6,6]);
line([min(lagtime),max(lagtime)],[bounds(1), bounds(1)]);
line([min(lagtime),max(lagtime)],[bounds(2), bounds(2)]);
refline(0,0), line([lagtime(ind),lagtime(ind)],[min(XCF)-0.05,max(XCF)+0.05],'color','k');
title('Cross-correlation between the signals (crosscorr)');
ylabel('Correlation');

subplot(2,1,2), plot(lagtime2, Rgr,'k')
hold on
plot(lagtime2, Rgr,'ko')
ylim([min(Rgr)-0.05,max(Rgr)+0.05]); xlim([-6,6]);
refline(0,0), line([lagtime2(ind2),lagtime2(ind2)],[min(Rgr)-0.05,max(Rgr)+0.05],'color','k');
xlabel('lag (s) [dopamine preceding]'); ylabel('Correlation');
title('Cross-correlation between the signals (xcorr)');


%% peak detection

%1) Designing a proper filter
kernel=exp(-[1:100]/30);

shift = (size(kernel,2)/2)/30; %the shift you get in the convolution equals half the size of the filter you are using


%2) Convolving the signal with it
convolvedDopamine = conv(dopamine,kernel,'same');
%the 'same flag' is to find the central part of the convolution that is the
%same size as "dopamine"

figure, set(gcf, 'Position', [800, 1000, 1000, 800])
subplot(3,2,1),
plot(time,zscore(dopamine),'g')
hold on
plot(time+shift,zscore(convolvedDopamine),'k') %crappy fix for now
xlim([500 520]);
xlabel('Time(s)');ylabel('Intensity');
legend('Original signal','Convolved signal');
title('Convolution with exponential filter');

%3) Taking a derivative of the signal
difConvolved = diff(convolvedDopamine);
subplot(3,2,3),
plot(time(2:length(time)),zscore(dopamine(2:length(time))),'g')
hold on
plot(time(2:length(time))+shift,zscore(difConvolved),'k') %again crappy fix.
%the time difference here seems to be more than in the previous case
xlim([500 520]);
xlabel('Time(s)');ylabel('Intensity');
legend('Original signal','Derivative of the convolved signal');
title('Derivative after convolution');

%4) Squaring the results
SqDifConvolved = diff(convolvedDopamine).^2;
subplot(3,2,5),
plot(time(2:length(time)),zscore(dopamine(2:length(time))),'g')
hold on
plot(time(2:length(time))+shift,zscore(SqDifConvolved),'k') %watch out with time
xlim([500 520]);
xlabel('Time(s)');ylabel('Intensity');
legend('Original signal','Square of the derivative of the convolved signal');
title('Squared derivative after convolution');


%5) Finding the peaks
zscored = zscore(SqDifConvolved);
threshDopamine = (max(SqDifConvolved)-min(SqDifConvolved))/6;
% it is actually not directly the threshold, but rather he amount above surrounding data
%for a peak to be identified


[peakLoc] = peakfinder(SqDifConvolved,threshDopamine);
subplot(3,2,[2,4,6]),
plot(time(2:length(time)),zscore(SqDifConvolved),'k')
hold on
plot(time(peakLoc+1),zscored(peakLoc),'ro') %why here does it change to +1??
xlim([500 550]);
xlabel('Time(s)');ylabel('Intensity');
title('Dopamine peaks')

%raster plot with the peaks
ones = repelem(0.03,size(peakLoc,1));

%check if it's ok to add the shift to the time axis
figure, set(gcf, 'Position', [800, 800, 800, 800]),
subplot(2,1,1),
plot(time(peakLoc+1)+shift,ones,'ko')
hold on
plot(time,D1,'r')
title('Dopamine peaks and D1 signal');
ylabel('Intensity');xlabel('Time(s)');

subplot(2,1,2)
plot(time(peakLoc+1)+shift,ones,'ko')
hold on
plot(time,D1,'r')
xlim([500 650]);
ylabel('Intensity');xlabel('Time(s)');


%6) Computing the dopamine firing rate for different "thresholds"

for i = 4:8
 threshDopamine(i-3) = (max(SqDifConvolved)-min(SqDifConvolved))/i;
 [peakLocT(i-3).peaks] = peakfinder(SqDifConvolved,threshDopamine(i-3));
 dopamineRate(i-3) = size(peakLocT(i-3).peaks)*30/size(SqDifConvolved);
end

figure, set(gcf, 'Position', [600, 600, 600, 600]),
plot(threshDopamine,dopamineRate,'ro');
xlabel('Threshold'); ylabel('Dopamine rate (Hz)');
title('Dopamine firing rate for different thresholds in the peak detection');
p = polyfit(threshDopamine,dopamineRate,2); %fitting a polynomial
x1 = linspace(min(threshDopamine)-0.25E-3,max(threshDopamine)+0.25E-3,300);
y1 = polyval(p,x1);
hold on
plot(x1,y1,'k')
xlim([min(threshDopamine)-0.25E-3,max(threshDopamine)+0.25E-3]);
l = polyfit(threshDopamine,dopamineRate,1);
yfit = l(1)*x1+l(2);
hold on;
plot(x1,yfit,'b-.');
legend('Data','Polynomial fit','Linear fit');


%% Spike-triggered average

stimulus = D1;

%we will generate a vector of events where the times with a peak are
%represented by a 1 and the rest is assigned to 0
events = zeros(size(dopamine));
events(peakLoc) = 1;


%likewise, check if it's fine to add the sift to the time axis
figure, set(gcf, 'Position', [800, 800, 800, 800]),
subplot(2,1,2), plot(time,stimulus)
ylabel('D1 intensity'), xlabel('Time (s)');
xlim([500 550]);
subplot(2,1,1), plot(time+shift,events, 'r')
ylabel('Dopamine "spikes"');
xlim([500 550]);


figure, set(gcf, 'Position', [800, 800, 800, 800]),
plot(time,stimulus)
ylabel('D1 intensity'), xlabel('Time (s)');
hold on
plot(time+shift,events*0.04, 'r')
xlim([500 550]);
legend('D1 signal','Dopamine peaks');


windowSize = 500; %number of frames to compute triggered average
[avg, avg2, nEvs] = evTrigAvg(events, stimulus, windowSize);

figure, set(gcf, 'Position', [800, 800, 800, 800]),
timeAxis = (0-windowSize/2:0+windowSize/2-1);
timeAxis = timeAxis/30;

plot(timeAxis,avg2,'r')
hold on
line([0,0],[min(min(avg2))-1E-3, max(max(avg2))+1E-3]);
ylim([min(min(avg2))-1E-3, max(max(avg2))+1E-3]);
title(['Spike triggered average of ', num2str(nEvs), ' windows'])
xlabel('Time, pts')
ylabel('Magnitude, arb.')

%the dopamine seems to be generating a deflexion in D1's activity 2s after
%a peak

%add the std

stdDopamine = std(avg,[],2);
tavg = avg';

figure,  stdshade(tavg,0.4,[],timeAxis)
hold on
line([0,0],[min(min(avg2))-1E-2, max(max(avg2))+1E-2]);
ylim([min(min(avg2))-1E-2, max(max(avg2))+1E-2]);
title(['Spike triggered average of ', num2str(nEvs), ' windows'])
xlabel('Time, pts')
ylabel('Magnitude, arb.')

%bootstrapping to get the confidence intervals:

ci = bootci(5000, @mean, tavg);
tci = ci';

figure,
plot(timeAxis,mean(avg,2),'k')
hold on
plot(timeAxis,tci(:,1),'r--');
plot(timeAxis,tci(:,2),'r--');
line([0,0],[min(min(tci))-1E-3, max(max(tci))+1E-3]);
ylim([min(min(tci))-1E-3, max(max(tci))+1E-3]);
xlim([min(timeAxis),max(timeAxis)]);
xlabel('lag time from dopamine peak (s)');
ylabel('Magnitude of the D1 signal (au)');
title(['Spike triggered average of ', num2str(nEvs), ' windows']);

%% Correlations between the velocity and the neural signals

%this part is not working.
%for some reason the cross corr are giving me complex numbers

%the velocity and 3d velocity have 1 nan that I have to take out.

Zspeed = Zscalars(:,5);
ZD1 = ZD1(all(~isnan(scalar),2),:);
Zdopamine = Zdopamine(all(~isnan(scalar),2),:);

[XCFsD1,lagssD1,boundssD1] = crosscorr(ZD1,Zspeed,220,3);
lagtimesD1 = lagssD1/30;

[RgrsD1, lags2sD1] = xcorr(Zspeed,ZD1,220,'coeff'); 
lagtime2sD1 = lags2sD1/30;
%the correlation is computed inversely that for crosscorr regarding which
%signal is considered first

%find the max crosscorr and xcorr values
[valsD1,indsD1] = max(XCFsD1);
[val2sD1,ind2sD1] = max(RgrsD1);

figure,  set(gcf, 'Position', [800, 1200, 1200, 800])
subplot(2,2,1), plot(lagtimesD1,XCFsD1, 'ro')
hold on
plot(lagtimesD1,XCFsD1, 'r')
line([min(lagtimesD1),max(lagtimesD1)],[boundssD1(1), boundssD1(1)],'Color','black','LineStyle','--');
line([min(lagtimesD1),max(lagtimesD1)],[boundssD1(2), boundssD1(2)],'Color','black','LineStyle','--');
ylim([min(XCFsD1)-0.05,max(XCFsD1)+0.05]);
refline(0,0), line([lagtimesD1(indsD1),lagtimesD1(indsD1)],[min(XCFsD1)-0.05,max(XCFsD1)+0.05],'color','k');
title('Cross-correlations for the D1 signal and speed (crosscorr)');
ylabel('Correlation');xlabel('Lag (s)');

subplot(2,2,3), plot(lagtime2sD1, RgrsD1,'k')
hold on
plot(lagtime2sD1, RgrsD1,'ko')
ylim([min(RgrsD1)-0.05,max(RgrsD1)+0.05]);
refline(0,0), line([lagtime2sD1(ind2sD1),lagtime2sD1(ind2sD1)],[min(RgrsD1)-0.05,max(RgrsD1)+0.05],'color','k');
xlabel('lag (s)'); ylabel('Correlation');
title('Cross-correlations for the D1 signal and speed (xcorr)');
xlim([-1 1]);


[XCFsDop,lagssDop,boundssDop] = crosscorr(Zdopamine,Zspeed,120,3);
lagtimesDop = lagssDop/30;

[RgrsDop, lags2sDop] = xcorr(Zspeed,Zdopamine,120,'coeff'); 
lagtime2sDop = lags2sDop/30;

[valsDop,indsDop] = max(XCFsDop);
[val2sDop,ind2sDop] = max(RgrsDop);


subplot(2,2,2), plot(lagtimesDop,XCFsDop, 'ro')
hold on
plot(lagtimesDop,XCFsDop, 'r')
line([min(lagtimesDop),max(lagtimesDop)],[boundssDop(1), boundssDop(1)],'Color','black','LineStyle','--');
line([min(lagtimesDop),max(lagtimesDop)],[boundssDop(2), boundssDop(2)],'Color','black','LineStyle','--');
ylim([min(XCFsDop)-0.05,max(XCFsDop)+0.05]);
refline(0,0), line([lagtimesDop(indsDop),lagtimesDop(indsDop)],[min(XCFsDop)-0.05,max(XCFsDop)+0.05],'color','k');
title('Cross-correlations for the dopamine signal and speed (crosscorr)');
ylabel('Correlation'); xlabel('Lag (s)');

subplot(2,2,4), plot(lagtime2sDop, RgrsDop,'k')
hold on
plot(lagtime2sDop, RgrsDop,'ko')
ylim([min(RgrsDop)-0.05,max(RgrsDop)+0.05]);
refline(0,0), line([lagtime2sDop(ind2sDop),lagtime2sDop(ind2sDop)],[min(RgrsDop)-0.05,max(RgrsDop)+0.05],'color','k');
xlabel('lag (s)'); ylabel('Correlation');
title('Cross-correlation for the dopamine signal and speed (xcorr)');
xlim([-1 1]);xlabel('Lag (s)');


%% Compute all of the pearson correlations between the signals and the scalars or their derivatives

%for normal scalars
for i = 1:size(scalars,2)
[R{i},P{i}] = corrcoef(scalars(:,i),ZD1,'rows','complete'); %the flags exclude rows with nans
[Rdop{i},Pdop{i}] = corrcoef(scalars(:,i),Zdopamine,'rows','complete');
r(i,1) = R{i}(2);
r(i,2) = Rdop{i}(2);
end

%for D1 and scalar derivatives
for i = 1:size(scalars,2)
[RD{i},PD{i}] = corrcoef(Dscalars(:,i),ZD1(2:end,:),'rows','complete');
[RDdop{i},PDdop{i}] = corrcoef(Dscalars(:,i),Zdopamine(2:end,:),'rows','complete');
rD(i,1) = RD{i}(2);
rD(i,2) = Rdop{i}(2);
end

figure,set(gcf, 'Position', [800, 1200, 1200, 800])
subplot(2,1,1);
b = bar(r)
set(b(1),'facecolor','r')
set(b(2),'facecolor','g')
ylabel('Pearson correlation coefficient');xlabel('Variables');
title('Correlation between the neural signals and different variables');
set(gca,'xtick',[1:6],'xticklabel',scalarNames)
legend('D1','Dopamine');
ylim([-0.15,0.15]);

subplot(2,1,2);
b2 = bar(rD)
set(b2(1),'facecolor','r')
set(b2(2),'facecolor','g')
ylabel('Pearson correlation coefficient');xlabel('Variables');
title('Correlation between the neural signal and different variable derivatives');
set(gca,'xtick',[1:6],'xticklabel',scalarNames)
legend('D1','Dopamine');
ylim([-0.15,0.15]);
