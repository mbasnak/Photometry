%% FIBER PHOTOMETRY ANALYSIS
%%
close all; clear all

%open the files
load('C:\Users\melanie\Documents\Doctorado\rotations\Bob\example data/example_photometry');
%is the data collected from the fiber photometry set-up and automatically
%saved in the format of this data file?

%% traces field

%raw traces
traces = phot_obj.traces;
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
suptitle('Raw traces')
for i = 1:size(traces,2)
subplot(2,2,i)
plot(traces(i).raw);
end

%baseline
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
suptitle('Baseline')
for i = 1:size(traces,2)
subplot(2,2,i)
plot(traces(i).baseline);
end

%(raw-baseline)/baseline
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
suptitle('Baseline substracted')
for i = 1:size(traces,2)
subplot(2,2,i)
plot(traces(i).raw-traces(i).baseline,'r');
end

%it looks like dff is the difference between the raw and the baseline.
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
suptitle('dff')
for i = 1:size(traces,2)
subplot(2,2,i)
plot(traces(i).dff,'k');
end

%for the reference...
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
plot(traces(1).dff);
title('trace');
subplot(2,2,2)
plot(traces(1).reference);
title('reference');
subplot(2,2,3)
plot(traces(4).dff);
subplot(2,2,4)
plot(traces(4).reference);
%it says reference channel is channel 3. is that an extra fiber that is not
%collecting information we'll use other than as a reference?


figure,
plot(traces(1).dff_reref)
%there is some sort of temporal rearrangement? does it use the
%reference_scale somehow?


%% options field

%it looks like there are already filters that are applied to the data
%before it is saved?

%and that ICA is used to rereference?
