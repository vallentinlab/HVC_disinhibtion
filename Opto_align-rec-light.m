%% plotting optogenetic neuronexus recordings
%yStim = board_adc_data = stim events
%% extract single pulses, block starts etc from the stimulus (yStim)

clear all
close all

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


adc_denoised=smoothdata(board_adc_data); %smooth data to remove noise spikes


thresh = 3; % threshold for detection of TTL pulses should be between 0 and 5
% find positiv threshold crossings:

ups = find(adc_denoised(1:end-1)<thresh & adc_denoised(2:end)>thresh);

% find negitave threshold crossings (useless, since only onsets are marked):
downs = find(adc_denoised(1:end-1)>thresh & adc_denoised(2:end)<thresh);

% the stimulus channel has some fast on-off dynamcis, likely samplin artefact - check the histogram to discard short pulses 
interThreshInterv = diff(ups);
figure
histogram(interThreshInterv,0:20:max(interThreshInterv));
tooShortThresh = 1000; % minmal sample distance of two ups - determined by looking at the histogram

ups([inf diff(ups)] < tooShortThresh) = []; % remove those ups that have not sufficinet samples between them and teh one before
downs([diff(downs) inf]< tooShortThresh) = []; % same for the downs

stimnumber=(length(ups));

figure
plot(adc_denoised), hold on, plot(ups,1,'ro'),plot(downs,1,'k*') % to check correct onset detection
%%%

% chunk this in blocks with constant repetition rate
figure 
histogram(diff(ups), 50:50:max(diff(ups))) %to determine border between blocks
blockThresh = 5e4; % the border between different blocks, visually determined
blockStartIdcs = find(diff([0 ups])>blockThresh); %at which stimnumber do blocks start
nBlocks = numel(blockStartIdcs); %how many different stimulusblocks are there
temp = diff(ups);
ipi = (temp(blockStartIdcs))/30000;  % inter pulse interval in each block in s

% %% add spikes - saved in array samplespikes, first level = channels, second
% level = channels per unit, unitsperchannel indicates sorted units

for k=1:length(samplespikes); %cycle through channels
    for l=1:length(1:unitsperchannel(k)); %cycle through units in channels
%align the spikes to the previous pulse start, define unit 
spikeTimes = samplespikes{1,k}{l,1}(samplespikes{1,k}{l,1}>ups(1)); % discard spikes before the first pulse
spikeTimes=transpose(spikeTimes);
upIdx = arrayfun(@(t) find((t-ups)>0,1,'last'),spikeTimes);
relSpikeT = (spikeTimes - ups(upIdx))/30000; %transform Sample ID to seconds

%for PSTH
numbins = round(100*0.5*max(ipi)); %set bins to 10ms
blockEdges = [blockStartIdcs +stimnumber]; %where are the blocks located, 
hEdges = linspace(0,0.5*max(ipi),numbins); % histogram bins

% % basic raster plot
figure
scatter(relSpikeT,upIdx','k.');
ylim([-5 stimnumber]) 
xlim([-0.01 0.5*max(ipi)])
hurz=histc(relSpikeT, hEdges); %define histogram data
line(hEdges, hurz, 'color', 'r') %plot histogram

% block marker lines
line(xlim,blockStartIdcs.*[1;1]-1,'color','k');
hold
% stimulus start marker
line([0 0],ylim,'color','b');

for iBl = 1:(nBlocks-1)
    % stimulus offset (estimated)
    line(ipi(iBl)/3*[1,1], [blockStartIdcs(iBl) blockStartIdcs(iBl+1)]-1,'color','b')
    
end
% last one.
line(ipi(end)/3*[1,1], [blockStartIdcs(end) stimnumber+1]-1,'color','b')
%%
%add PSTH
    
% 
% for iBl = 1:nBlocks
%     
%     h(iBl) = subplot(nBlocks,1,nBlocks-iBl+1);
%     thisSpikes = relSpikeT(upIdx >blockEdges(iBl) & upIdx < blockEdges(iBl+1));
%     histogram(thisSpikes,hEdges);
%     axis tight
%     
%     xlim([-1e3 hEdges(end)])
%     box off
%     set(gca,'xtick',0,'XTickLabel',[],'Yticklabel',[],'TickDir','out')
%     if iBl ==1
%         set(gca,'xtick', 0)
%     end
% end
% 
% linkaxes(h)
%%
set(gca,'TickDir','out')
ylabel('Pulse Number')
xlabel('Time [s]')
channeltitle = sortedspikes_names{1,k};
channeltitle=replace(channeltitle,'_',' ');
spikefile=replace(spikefile,'_',' ');
currentchanneldepth=channeldepth(k);

title({spikefile(1:end-19),channeltitle,'Unit = ' l, 'Depth = ' channeldepth(k)})
        end
end