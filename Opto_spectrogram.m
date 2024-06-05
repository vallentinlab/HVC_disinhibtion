%plot spectrograms

spectrogram(data,1024,1000,1024,fs,'yaxis')

ylim([0 8])

 colorbar off
myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
caxis([-75 -54])
%first = what is visualised (lower = less); second = intensity (lower = less)
%xlim([0 3.5])

%add waveformplot
% dt=1/fs
% t = 0:dt:(length(data)*dt)-dt;

% figure
% plot(t,data)