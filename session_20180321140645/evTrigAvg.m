function [avg, avg2, nEvs] = evTrigAvg(events, stim, windowSize)

% Discard spikes in first window
events(1:windowSize/2)=0; %the events in the first time window are set to zero
events(end-windowSize/2:end)=0; %I'm also setting the final ones to zero, 
%because I am making a window centered around the events.

%because it would be impossible to access information of stim(-750:750)
%for example if we were are analyzing at events(750).

% Find events: 
% Number of events
nEvs = sum(events);
% Indexes of events
evIdx = find(events);

% Preallocate average and std
%avg = zeros(1,windowSize);
avg = zeros(windowSize,nEvs);


% For each event
for w = 1:nEvs
    % Find the indexes of the time window preceding the event
    wIdx = evIdx(w)-(windowSize/2) : evIdx(w)+(windowSize/2)-1;
    % Add the stim from this window to the average
    %avg = avg + stim(wIdx);
    avg(:,w) = avg(:,w) + stim(wIdx);
end
% Divide by number of events to complete average
%avg=avg./sum(events);
avg2 = mean(avg,2);
