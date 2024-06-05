%% plotting optogenetic neuronexus recordings
%% extract single pulses, block starts etc from the stimulus (yStim)
% yStim = board_adc_data = stim events

%colors = 0 .66 1 for 470nm | 1 .81 0 for 595nm

clear all
close all

%SET Y LIMITS HERE IF YOU WANT TO PLOT ONLY A BLOCK OF TRIALS
stimlow = -5; %standard = -5
stimup = 180; %standard = 180, max trials

%for significance testing
testingparametersall = [];
testingparameters = cell(5,20);

%what is the distance of each channel A_008 until A_023 [um]?
channeldepthraw=[75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450];


%select file with sorted spikes
disp('Please select the SORTED SPIKES');
[filename, pathname]=uigetfile('*.mat');
spikefile=fullfile(pathname, filename);
load(spikefile);
spikefile=filename(1:end-4);

%select ADC to corresponding file
disp('Please select the ADC data (pulse trains)');
[filename, pathname]=uigetfile('*.mat');
pulsefile=fullfile(pathname, filename);
load(pulsefile);
pulsefile=filename(1:end-4);

%merge channels
sortedspikes={A_008 A_009 A_010 A_011 A_012 A_013 A_014 A_015 A_016 A_017 A_018 A_019 A_020 A_021 A_022 A_023};
sortedspikes_names={'A_008' 'A_009' 'A_010' 'A_011' 'A_012' 'A_013' 'A_014' 'A_015' 'A_016' 'A_017' 'A_018' 'A_019' 'A_020' 'A_021' 'A_022' 'A_023'}; %save channel names
sortedspikeslengthall=length(sortedspikes);
empties = find(cellfun(@isempty,sortedspikes)); % identify the empty cells
sortedspikes(empties) = []  ; % remove the empty cells
sortedspikes_names(empties)= []; % remove empty cell names
sortedspikeslengthfull=length(sortedspikes);

%how deep was each channel?
depth=input('How deep did you insert the electrode [micrometer]');
channeldepth=depth-channeldepthraw; %adjust to insertiondepth
channeldepth(empties) = []  ; % remove the empty cells

%select pre silence for better plotting
pretimesec=input('Please enter your pre pulse silence [S]: ');
pretimesamp = pretimesec * 30000; %convert to samples

%loop all channels and split units
for k=1:length(sortedspikes)
curchannel=sortedspikes{k};
[~,~,X] = unique(curchannel(:,2));
splitunits{k} = accumarray(X,1:size(curchannel,1),[],@(r){curchannel(r,:)});
end

%determine number of units and channels
numchannels=length(splitunits);
unitsperchannel=zeros(1,numchannels);
for k=1:numchannels;
    %subtract 1 to avoid classification of unsorted waveforms (unit = 0) as
    %unit
    unitsperchannel(k)=length(splitunits{1,k});
end

%convert spiketimes
samplespikes=splitunits; % create new array
for k=1:length(splitunits); %cycle through channels
    for l=1:length(splitunits{1,k}); %cycle through units within channels
        samplespikes{1,k}{l,1}(:,3) = splitunits{1,k}{l,1}(:,3)*30000; %convert seconds to samples
        round(samplespikes{1,k}{l,1}(:,3));
    end
end

%create array for averages
averagearray = splitunits; % create new array
for k=1:length(averagearray); %cycle through channels
    for l=1:length(averagearray{1,k}); %cycle through units within channels
        averagearray{1,k}{l,1}(:,1:3) = []; %empty all rows
    end
end

adc_denoised=smoothdata(board_adc_data); %smooth data to remove noise spikes for better stim onset detection

%translate spiketimes to 1 in array of 0, necessary for PSTH

%process spike events to use as PSTH
samplelength=length(board_adc_data); %length of entire recorded file
count=[1:samplelength]; %number of samples in the recording
pstharray = samplespikes; %create new array with existing structure
for k=1:length(pstharray); %cycle through channels
    for l=1:length(pstharray{1,k}); %cycle through units within channels
        samplespikes{1,k}{l,1} (:,3)= ceil(samplespikes{1,k}{l,1} (:,3)); %remove decimal 0 from sample column to allow for precise comparison
        pstharray{1,k}{l,1}(:,1:2) = []; %delete rows
        pstharray{1,k}{l,1} = zeros(length(count),1); %fill with zeros in sample length
        pstharraycurrent=(ismember(count,samplespikes{1,k}{l,1}(:,3))); %plot vector of sample length with 1 at spike event
        pstharraycurrent=+pstharraycurrent; %transform vector from logical to double
        pstharray{1, k}{l,1}=pstharraycurrent; % current vector into correct cell of new array
    end
end

thresh = 3; % threshold for detection of TTL pulses should be between 0 and 5
% find positive threshold crossings:
ups = find(adc_denoised(1:end-1)<thresh & adc_denoised(2:end)>thresh); %at which sample are stim pulses located

% the stimulus channel has some fast on-off dynamcis, likely sampling artifact - check the histogram to discard short pulses 
interThreshInterv = diff(ups);
figure
histogram(interThreshInterv,0:20:max(interThreshInterv));
tooShortThresh = 1000; % minmal sample distance of two ups - determined by looking at the histogram, should be stable now

ups([inf diff(ups)] < tooShortThresh) = []; % remove those ups that have not sufficinet samples between them and teh one before

ups = ups-pretimesamp; %adjust to prepulse interval

stimnumber=(length(ups));

figure
plot(adc_denoised), hold on, plot(ups,1,'ro') % to check correct onset detection
hold off
%%%

% chunk this in blocks with constant repetition rate
figure 
histogram(diff(ups), 50:50:max(diff(ups))) %to determine border between blocks
blockThresh = 3e4; % the border between different blocks, visually determined, shorten or lengthen if issues with block separator
blockStartIdcs = find(diff([0 ups])>blockThresh); %at which stimnumber do blocks start
nBlocks = numel(blockStartIdcs); %how many different stimulusblocks are there
temp = diff(ups);
ipi = (temp(blockStartIdcs))/30000;  % inter pulse interval in each block in s


%% plot graphs for every sorted unit
for k=1:length(samplespikes); %cycle through channels
    
    for l=2:length(1:unitsperchannel(k)); %cycle through units in channels
%align the spikes to the previous pulse start, define unit 
spikeTimes = samplespikes{1,k}{l,1}(samplespikes{1,k}{l,1}>ups(1)); % discard spikes before the first pulse
spikeTimes=transpose(spikeTimes);
upIdx = arrayfun(@(t) find((t-ups)>0,1,'last'),spikeTimes);
relSpikeTpsth =  (spikeTimes - ups(upIdx)); %spiketimes in samples per unit
relSpikeTpsthmin = min(relSpikeTpsth); %begin of current unit
relSpikeTpsthmax = max(relSpikeTpsth); %end of current unit

relSpikeT = (((spikeTimes - ups(upIdx))/30000)-pretimesec); %align to stim time & transform Sample ID to seconds

%find pulses during last third of ipi, move them to front by subtracting
%two thirds of ipi - would this work? so far: extracted stims in last third
% relSpikeTcheck = relSpikeT>(min(ipi)*(2/3));
% relSpikeTpre = relSpikeT(relSpikeTcheck);
% relSpikeTpre = relSpikeTpre-(min(ipi)*1/3);
% upIdxpre=upIdx(length(relSpikeTpre):end);
%upIdx needs to be adjusted, too


% % basic raster plot
figure
scatter(relSpikeT,upIdx','k.');
ylim([stimlow stimup])
xlim([-pretimesec (min(ipi)*2/3)]) %set xlim max to min(ipi) to make sure that plot length resembles actual stimulus course

hold on
% stimulus start marker
line([0 0],ylim,'color','b');

% stimulus offset, estimated based on known stimulus design (here based on optostimnew:
% 100ms stim, 200ms silence)
line(mean(ipi)/3*[1,1],ylim,'color','b')

set(gca,'TickDir','out')
ylabel('Pulse Number')
xlabel('Time [s]')
channeltitle = sortedspikes_names{1,k};
channeltitle=replace(channeltitle,'_',' ');
spikefile=replace(spikefile,'_',' ');
currentchanneldepth=channeldepth(k);

title({spikefile(1:end-19),channeltitle,'Unit = ' l, 'Depth = ' channeldepth(k)})

%% PSTH
blockEdges = [blockStartIdcs +stimnumber]; %where are the blocks located,
hEdges = linspace(0,0.3,max(blockEdges)); % histogram bins across all stimuli

thisSpikes = relSpikeT(upIdx >min(blockEdges) & upIdx < max(blockEdges));

exitflag = size(thisSpikes); %prevent error on empty units

if exitflag > 0
    
    %new PSTH
    psthmatrix = zeros(stimnumber,round(max(ipi)*30000));
    for m = 1:length(relSpikeTpsth)
    psthmatrix(upIdx(m),relSpikeTpsth(m)) = 1;
    end
    
    psth=30000*smooth(mean(psthmatrix(stimlow+6:stimup,:)),2000,'loess'); %add 6 in standard conditions to prevent error
   

    yyaxis right
    plot(((1/30000:1/30000:length(psth)/30000)-pretimesec),psth,'r','LineWidth',2)
    psthmax = max(psth)+1;

    ylim([0 psthmax])
    ylabel('Average Firing Rate [Hz]')
    
    ax = gca;
    ax.YColor = 'r';    
    
else
end
%% extract individual datapoints

[trialcounts,trialvalues] = groupcounts(upIdx'); %extract number of spikes within each trial, extract corresponding number of trial
maxspikes = max(trialcounts); %max number of spikes

trialswspikes = nan(stimnumber, maxspikes); %create vector stimnumber by max number of spikes

for ii = 1:stimnumber
    flags = upIdx == ii;
    zz = sum(flags);
    trialswspikes(ii, 1:zz) = relSpikeT(flags);
end

%bug: trial 90 and trial 180 (last trials of each block) have too high
%firing rate, remove rows
trialswspikes([90 180],:) = [];

%avgmatrix=[averagearray{1, k}{l, 1}([1:stimnumber],1) averagearray{1, k}{l, 1}([1:stimnumber],2)]; %merge average array into single matrix

%% plot Average boxplots

%get data for current unit
%first column of avgmatrix = during stim, second column = after/between

avgspikes = nan(stimnumber-2, 2); %correct number here for subtracted trials 90 and 180

for ii = 1:stimnumber-2 %correct number for 90 & 180  
    spikesduring = sum(trialswspikes(ii,:)> 0) + sum(trialswspikes(ii,:)<mean(ipi)/3); %spikes during stim
    spikesafter = sum(trialswspikes(ii,:) < 0) + sum(trialswspikes(ii,:)>mean(ipi)/3 & trialswspikes(ii,:)<=(2*mean(ipi)/3)); %spikes after stim
    avgspikes(ii,1)=spikesduring;
    avgspikes(ii,2)=spikesafter;
end    

avgspikeshz=avgspikes;  
    
avgspikeshz(:,1) = avgspikes(:,1)./(mean(ipi)/3); %calculate spikerate in hz during
avgspikeshz(:,2) = avgspikes(:,2)./(2*(mean(ipi))/3); %calculate spikerate in hz after, double length

%avgspikeshzfix = avgspikeshz/4;

plotmatrix(:,1)=avgspikeshz(:,1); %during left
plotmatrix(:,2)=avgspikeshz(:,2); %after right

%this scatters points in violin on y-axis by 0-1hz BUGFIX, only in plot!
for rr = (1:length(plotmatrix))
    plotmatrix(rr,1) = plotmatrix(rr,1)+(4*rand(1));
    plotmatrix(rr,2) = plotmatrix(rr,2)+(4*rand(1));
end

avgplot = [avgspikeshz(:,1)' avgspikeshz(:,2)']; %merge data into one row for scattered lines
avgplot = avgplot';

figure
% connect scatter points with lines
xline = [ones(stimnumber-2,1); ones(stimnumber-2,1)*2]; % left and right definition for connected lines, number based on number of blocks 
line([xline(1:stimnumber-2) xline(stimnumber-1:end)]',[avgplot(1:stimnumber-2) avgplot(stimnumber-1:end)]', 'Color', 'k', 'LineStyle', '--')

hold on
Violin(plotmatrix(:,1), 1, 'ViolinColor',  [1 .81 0])
Violin(plotmatrix(:,2), 2, 'ViolinColor', [0 0 0])
boxplot([plotmatrix], 'Labels', {'Light ON', 'Light OFF'}, 'Whisker', 8, 'Colors', 'k')

ylabel('Average Firing Rate [Hz]');

%% paired t-Test test to compare difference during and after light

during =avgspikeshz(:,1); 
after =avgspikeshz(:,2); 

[h,p,ci,stats] = ttest(during,after) %test for difference

stats = struct2cell(stats);
h = double(h);

annotation('textbox', [0.67, 0.83, 0.1, 0.1], 'String', "p-Value is " + p, 'LineStyle', 'none')
annotation('textbox', [0.67, 0.78, 0.1, 0.1], 'String', "T-Value is " + stats{1,1}, 'LineStyle', 'none')

hold off

testingparameters{l,k}={[p,h],stats};    
    end
end

%% save current figures

openfigs =  findobj('type','figure'); %count current figures
fignumber = length(openfigs);

folder = uigetdir();

for f=1:fignumber
    baseFileName = sprintf('improvedfig_%d.emf',f);
    figure(f); % Activate the figure again.
    saveas(gcf,fullfile(folder, baseFileName),'emf')
end

testingname = '\testingparameters.mat';
testparamsave = strcat(folder,testingname);

save(testparamsave,'testingparameters')

close all