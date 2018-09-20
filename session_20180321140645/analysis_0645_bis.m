%analysis for MB14 fiber photometry, session 1

close all;clear all;
cd '/n/groups/datta/mel/d1+snc/session_20180321140645';


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
xlabel('Time (frame)');
xlim([0 size(phot.traces(1).dff,1)]);

% the raw traces show the bleaching, and the dff ones are already corrected

% looking at ext.projections you find the length of the different fields (ie number of frames) for the kinect data
% doing length (phot.traces(1).dff) you find the number of frames for the tdt data.
% the kinect data can have more frames and have NaNs. To overcome this, we can do

vel_frames = ext.get_original_timebase(ext.projections.velocity_mag);
% for example, to get the proper number of frames for the velocity

% save figure
saveas(gcf,'PhotSignals.jpg')

%% getting the scalars and their derivatives

%taking the data for the different scalars and putting them in an array
scalars(:,1) = ext.get_original_timebase(ext.projections.angle);
scalars(:,2) = ext.get_original_timebase(ext.projections.width);
scalars(:,3) = ext.get_original_timebase(ext.projections.length);
scalars(:,4) = ext.get_original_timebase(ext.projections.height_ave);
scalars(:,5) = ext.get_original_timebase(ext.projections.velocity_mag);
scalars(:,6) = ext.get_original_timebase(ext.projections.velocity_mag_3d);

%saving the names
scalarNames = {'Angle','Width','Length','Height','Velocity','3DVelocity'};

%remove the frames with Nans in the photometry from the kinect data (for example, the vel_frames)
scalars = scalars(all(~isnan(phot.traces(4).dff),2),:);

%convolve the scalars with a boxcar filter and then take the derivative
for i=1:size(scalars,2)
convscalars(:,i) = conv(scalars(:,i),ones(30,1)/30,'same');
Dscalars(:,i) = diff(convscalars(:,i));
end

%zscoring them
Zscalars = zscore(convscalars);
ZDscalars = zscore(Dscalars);


%% close-up to the photometry signals

%remove the rows in the photometry data with with nans:
dopamine = phot.traces(1).dff(all(~isnan(phot.traces(1).dff),2),:); 
D1 = phot.traces(4).dff(all(~isnan(phot.traces(4).dff),2),:);

%remove the same frames in the kinect data (for example, the vel_frames)
velocity = vel_frames(all(~isnan(phot.traces(4).dff),2),:);

%re-scale the time axis to seconds
time = [1/30:1/30:size(velocity,1)/30];


%close-up of the signals
figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(3,1,1)
plot(time,dopamine, 'green')
%hold on
%refline(0,mean(dopamine))
title('Close-up to signal from dopaminergic terminals');
ylabel('Intensity');
xlim([100 200]);

subplot(3,1,2)
plot(time,D1, 'red')
%hold on
%refline(0,mean(D1))
title('Close-up to signal from D1 neurons');
ylabel('Intensity');
xlim([100 200]);

subplot(3,1,3)
plot(time,velocity, 'black')
title('velocity');
xlabel('Time (s)'); ylabel('Velocity');
xlim([100 200]);

saveas(gcf,'CloseupPhotSignals.svg');

%% spectral analysis

Fs = 30; %sampling frequency (30 Hz)

FTg = fft(dopamine-mean(dopamine)); %calculate the fast fourier transform
%substract the mean to not get a huge power for a 0Hz frequency
magFTg = abs(FTg); %take the absolute value
fx = 0:(size(dopamine,1)/2)-1; %build up a proper frequency axis
fx = (fx*Fs)/size(dopamine,1); %scaling it so that it represents the frequencies in Hz

figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(2,1,1)
plot(fx,magFTg(1:size(dopamine,1)/2),'g') 
grid
title('Spectrum of the signal from dopaminergic terminals');
ylabel('Magnitude (dB)');

FTr = fft(D1-mean(D1)); %calculate the fast fourier transform
magFTr = abs(FTr);
subplot(2,1,2)
plot(fx,magFTg(1:size(D1,1)/2),'r')
grid
title('Spectrum of the signal from D1 neurons');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

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


%compute coherence between the signals
[cxy,fc] = mscohere(dopamine,D1,hamming(512),500,2048);
%Plot the coherence function and overlay the frequency responses of the filters.
[qx,fd] = freqz(dopamine);
qy = freqz(D1);

figure, %plot it
plot(fc/pi,cxy)
hold on
plot(fd/pi,abs(qx),fd/pi,abs(qy))
hold off


%another way of looking at it.
figure, set(gcf, 'Position', [800, 800, 800, 800])
subplot(2,1,1),
plot(f,20*log10(mean(abs(s),2)),'g')
ylabel('Power/frequency (dB/Hz)');
hold on
plot(f2,20*log10(mean(abs(s2),2)),'r')


%take the derivative of the curves to see if they have similar
%characteristics

a=diff(20*log10(mean(abs(s),2)));
b=diff(20*log10(mean(abs(s2),2)));
dist=norm(a-b);

subplot(2,1,2),
plot(f(2:size(f,1)),a,'g')
xlabel('Frequency (Hz)'); ylabel('Derivative of Power/frequency');
hold on
plot(f2(2:size(f2,1)),b,'r')

[ax,h3]=suplabel('Spectral content of the signals' ,'t'); 


%saveas(gcf,'Spectrum.png');
saveas(gcf,'Spectrum.svg');

%% cross-correlation

% z-score the signals
Zdopamine = zscore(dopamine);
ZD1 = zscore(D1);

% construct an Elliptic high-pass filter with a .5 Hz cutoff
% 3rd order, .2 ripple, 40 dB attenuation
max_lag=500;
[b,a]=ellip(3,.2,40,.25/(30/2),'high');

%calculate the cross-corr between the two filtered signals
[r,lags]=xcorr(filtfilt(b,a,ZD1),filtfilt(b,a,Zdopamine),max_lag,'coeff'); 

%shift the dopamine signal by a random value 100 times and compute the
%cross-corr with the normal D1 signal to determine significance boundaries
nrands=100;
rands_rs=nan(nrands,max_lag*2+1);

for i=1:nrands
    rands_rs(i,:)=xcorr(filtfilt(b,a,ZD1),...
        circshift(filtfilt(b,a,Zdopamine),randi(length(dopamine)*2,1)),max_lag,'coeff');
end

% plot the result
figure(); set(gcf, 'Position', [800, 1000, 1000, 800]),
plot(lags/30,prctile(rands_rs,[2.5 97.5])','r--','LineWidth',2) % plot the 2.5 and 97.5 percentiles of the randomized cross-corr to get confidence intervals
hold on
plot(lags/30,r,'k','LineWidth',2);
[val,ind] = max(r); % obtain the max corr value
xlabel('lags (s)'); ylabel('Correlation');
title('Cross-correlation between the dopamine and D1 signals');
xlim([min(lags/30),max(lags/30)]); ylim([-0.1,0.2]);
text((lags(ind)/30)+0.5,max(r),strcat('max corr at lag ',' ', num2str(lags(ind)/30),' s'));

%saveas(gcf,'CrossCorr.jpg');
%saveas(gcf,'CrossCorr.png');
saveas(gcf,'CrossCorr.svg');

%% Peak finder

rectZopamine = Zdopamine; 
rectZopamine = filtfilt(b,a,Zdopamine); %filter the signal
rectZopamine(rectZopamine<0) = 0; %rectify it

diffrectZopamine = diff(rectZopamine); %take the derivative
sqdiffrectZopamine = diffrectZopamine.^2; %square it

[vals,locs] = findpeaks(sqdiffrectZopamine,'minpeakheight',5*mad(sqdiffrectZopamine)/.6745,'minpeakdistance',5); %find the peaks
dopamine_vals = dopamine(locs); %find the dopamine values at the peaks

% plot the peak detection
figure, set(gcf, 'Position', [800, 1000, 1000, 800])
subplot(3,2,1),
plot(time,Zdopamine,'g')
hold on
plot(time,rectZopamine,'k')
xlim([500 520]);
xlabel('Time(s)');ylabel('Intensity');
legend('Original signal','Convolved rectified signal');
title('Convolution');

subplot(3,2,3),
plot(time(2:length(time)),Zdopamine(2:length(time)),'g')
hold on
plot(time(2:length(time)),zscore(diffrectZopamine),'k') 
xlim([500 520]);
xlabel('Time(s)');ylabel('Intensity');
legend('Original signal','Derivative of the convolved signal');
title('Derivative after convolution');

subplot(3,2,5),
plot(time(2:length(time)),Zdopamine(2:length(time)),'g')
hold on
plot(time(2:length(time)),zscore(sqdiffrectZopamine),'k') %watch out with time
xlim([500 520]);
xlabel('Time(s)');ylabel('Intensity');
legend('Original signal','Square of the derivative of the convolved signal');
title('Squared derivative after convolution');

subplot(3,2,[2,4,6]),
plot(time(2:length(time)),zscore(sqdiffrectZopamine),'k')
hold on
zscored = zscore(sqdiffrectZopamine);
plot(time(locs+1),zscored(locs),'ro') %why here does it change to +1??
xlim([500 550]);
xlabel('Time(s)');ylabel('Intensity');
title('Dopamine peaks')

saveas(gcf,'PeakFinder.svg');

% plot the peak distribution
figure,
h=histogram(vals)
xlabel('Dopamine peak amplitude');ylabel('Frequency');
title('Dopamine peak distribution');
saveas(h,'DopamineDistribution.svg');

%% STA of the D1 signal with the dopamine peaks

win_size=150; %set the window size at 150 frames
locs(locs<=win_size|locs>=length(dopamine)-win_size)=[]; %empty the peaks before and after the windowsize
thresholds=prctile(dopamine(locs),[25 50 75]); %set three different thresholds as the 25, 50 and 75 percentile of the dopamine peak distribution
dopamine_pk_vals=dopamine(locs); %get the dopamine values for all the peaks
%locs(dopamine_pk_vals<thresholds(1))=[]; %in case we want to select the
%values up to a certain percentile

wins=nan(win_size*2+1,length(locs));

filtered_d1=filtfilt(b,a,D1);

for i=1:length(locs)
    wins(:,i)=filtered_d1(locs(i)-win_size:locs(i)+win_size);  %get the filtered d1 signal 150 frames before and after the dopamine peaks  
end

time = ((1:size(wins,1))-(size(wins,1)/2))/30;

STAci = bootci(5000, @mean, wins'); %bootstrap to get confidence intervals
tSTAci = STAci'; %transpose the matrix

figure, set(gcf, 'Position', [800, 800, 800, 800]),
plot(time', mean(wins,2),'k','LineWidth',2)
xlabel('lag (s)'); ylabel('D1 signal');
title('Spike-triggered average of the D1 signal with the dopamine peaks');
hold on
plot(time',tSTAci(:,1),'r--');
plot(time',tSTAci(:,2),'r--');

saveas(gcf,'STA.svg');


%% Correlations between the velocity and the neural signals

Zspeed = zscore(velocity); %zscore speed

%get corr for the different lags between speed and D1 signal
[RgrsD1, lags2sD1] = xcorr(Zspeed,ZD1,220,'coeff'); 
lagtime2sD1 = lags2sD1/30;
[val2sD1,ind2sD1] = max(RgrsD1);

%plot it
figure,  set(gcf, 'Position', [800, 1200, 1200, 800])
subplot(1,2,1), 
plot(lagtime2sD1, RgrsD1,'k')
hold on
plot(lagtime2sD1, RgrsD1,'ko')
%ylim([min(RgrsD1)-0.05,max(RgrsD1)+0.05]);
refline(0,0);
xlabel('lag (s)'); ylabel('Correlation');
title('Cross-correlations for the D1 signal and speed (xcorr)');
xlim([-2 2]); ylim([-0.15 0.15]);

%get corr between speed and dopamine signal
[RgrsDop, lags2sDop] = xcorr(Zspeed,Zdopamine,220,'coeff'); 
lagtime2sDop = lags2sDop/30;
[val2sDop,ind2sDop] = max(RgrsDop);

%plot it
subplot(1,2,2), plot(lagtime2sDop, RgrsDop,'k')
hold on
plot(lagtime2sDop, RgrsDop,'ko')
%ylim([min(RgrsDop)-0.05,max(RgrsDop)+0.05]);
refline(0,0);
xlabel('lag (s)'); ylabel('Correlation');
title('Cross-correlation for the dopamine signal and speed (xcorr)');
xlim([-2 2]);xlabel('Lag (s)'); ylim([-0.15 0.15]);

%% Compute all of the pearson correlations between the signals and the scalars or their derivatives

%for normal scalars
[b,a]=ellip(3,.2,40,.5/15,'high'); %build an elliptic filter
rectZopamine = Zdopamine; 
rectZopamine = filtfilt(b,a,rectZopamine); %filter the signal
rectZopamine(rectZopamine<0) = 0; %rectify it

rectD1=ZD1;
rectD1=filtfilt(b,a,rectD1); %filter the D1 signal
rectD1(rectD1<0)=0; %rectify it

for i = 1:size(scalars,2) 
[R{i},P{i}] = corrcoef(convscalars(:,i),rectD1); %compute the corr coeff between the convolved scalars and the filtered, rectified D1 signal
[Rdop{i},Pdop{i}] = corrcoef(convscalars(:,i),rectZopamine); %same for the dopamine signal
co(i,1) = R{i}(2);
co(i,2) = Rdop{i}(2);
end

%for the signals and derivatives
for i = 1:size(scalars,2)
[RD{i},PD{i}] = corrcoef(Dscalars(:,i),ZD1(2:end,:));
[RDdop{i},PDdop{i}] = corrcoef(Dscalars(:,i),Zdopamine(2:end,:));
rD(i,1) = RD{i}(2);
rD(i,2) = Rdop{i}(2);
end
% why are we not using the filtered, rectified signal for this?

% plot the correlations
figure,set(gcf, 'Position', [800, 1200, 1200, 800])
subplot(2,1,1);
b = bar(co)
set(b(1),'facecolor','r')
set(b(2),'facecolor','g')
ylabel('Pearson correlation coefficient');xlabel('Variables');
title('Correlation between the neural signals and different variables');
set(gca,'xtick',[1:6],'xticklabel',scalarNames)
legend('D1','Dopamine');
ylim([-0.1,0.1]);

subplot(2,1,2);
b2 = bar(rD)
set(b2(1),'facecolor','r')
set(b2(2),'facecolor','g')
ylabel('Pearson correlation coefficient');xlabel('Variables');
title('Correlation between the neural signal and different variable derivatives');
set(gca,'xtick',[1:6],'xticklabel',scalarNames)
legend('D1','Dopamine');
ylim([-0.1,0.1]);

saveas(gcf,'CorrBarPlot.svg');

%%   correlation between scalars and d1 taking different bins

samples_idx=1:length(rectD1); %vector with frame numbers
samples_t=samples_idx/30; %vector with time values

bin_sizes=[1:200]; %vector with the different bin sizes to try out
%bin_sizes=[1:10 15:5:200];
%bin_rs=nan(1,length(bin_sizes));



for i=1:length(bin_sizes) %for every bin size
   for j=1:6 %for every scalar
    time_bins=samples_t(1):bin_sizes(i):samples_t(end);
    time_bins(end)=inf; %set the last time point as inf. why?
    [~,bin_idx]=histc(samples_t,time_bins);
    %counts the number of values in samples_t that are within each specified bin range.
    %The input, time_bins, determines the endpoints for each bin.
    %The output, bin_idx, contains 
    
    y=accumarray(bin_idx(:),rectD1(:));
    %A = accumarray(subs,val) returns array A by accumulating elements of vector val using the subscripts subs.
    y2{:,j}=accumarray(bin_idx(:),scalars(:,j));
    bin_rs{i,j}=corr(y,y2{:,j});
    end
end

bin_rs2 = cell2mat(bin_rs);

figure, ,set(gcf, 'Position', [800, 1200, 1200, 800])
for p = 1:6
    subplot(1,6,p),
    plot(bin_sizes,bin_rs2(:,p),'color',rand(1,3))
    name = scalarNames(p);
    title(name);
    ylim([-1 1]);
end
[ax,h3]=suplabel('Correlations between the D1 signal and the scalars as a function of binning' ,'t'); 
[ax,h1]=suplabel('Bin size (s)');
[ax,h3]=suplabel('Correlation coefficient','y');

saveas(gcf,'corrD1scalars.svg');

%% correlation between scalars and dopamine taking different bins


for i=1:length(bin_sizes) %for every bin size
   for j=1:6
    time_bins=samples_t(1):bin_sizes(i):samples_t(end);
    time_bins(end)=inf;
    [~,bin_idx]=histc(samples_t,time_bins);
    y_dop=accumarray(bin_idx(:),rectZopamine(:));
    y2{:,j}=accumarray(bin_idx(:),scalars(:,j));
    bin_rs_dop{i,j}=corr(y_dop,y2{:,j});
    end
end

bin_rs2_dop = cell2mat(bin_rs_dop);

figure, ,set(gcf, 'Position', [800, 1200, 1200, 800])
for p = 1:6
    subplot(1,6,p),
    plot(bin_sizes,bin_rs2_dop(:,p),'color',rand(1,3))
    name = scalarNames(p);
    title(name);
    ylim([-1 1]);
end
[ax,h3]=suplabel('Correlations between the dopamine signal and the scalars as a function of binning' ,'t'); 
[ax,h1]=suplabel('Bin size (s)');
[ax,h3]=suplabel('Correlation coefficient','y');

saveas(gcf,'corrdopscalars.svg');