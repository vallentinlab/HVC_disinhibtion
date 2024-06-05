%plot change of syllable numbers ARCHT/CTRL

columnIDs={'id','treatment','sylpre','syl4wkpost'}

%create data
syllablenumbers={'m878', 'archt', 6, 9;...
    'm864','archt',6,13;...
    'm313','archt',4,9;...
    'm322','archt',5,10;...
    'm146','archt',6,10;...
    'm717','archt',5,11;...
    'm319','ctrl',5,5;...
    'm362','ctrl',7,7;...
    'm871','ctrl',5,5;...
    'm737','ctrl',7,7};

syllablenumbers=cell2table(syllablenumbers,"VariableNames",columnIDs);

%get difference between pre/post
syllablenumbers.diff=syllablenumbers.syl4wkpost-syllablenumbers.sylpre;

%create x value for scatter with jitter
x=-0.2+(0.2+0.2)*rand(length(syllablenumbers.diff),1);

%plot with different colours based on treatment
gscatter(x,syllablenumbers.diff,syllablenumbers.treatment,'yk');

xlim([-0.4 0.4])
ylim([-1 10])
ylabel('Novel elements four weeks post playback')

%plot data as histo
archtdif = table2array(syllablenumbers(1:6,5));
ctrldif = table2array(syllablenumbers(7:10,5));

figure
hnewctrl=histogram(ctrldif,'edgecolor','k', 'DisplayStyle', 'stairs', 'LineWidth',2, 'EdgeAlpha', 0.5);
hold on
hnewArchT=histogram(archtdif,'edgecolor','y', 'DisplayStyle', 'stairs', 'LineWidth',2, 'EdgeAlpha', 0.5);

binsize=1;

% hregCtrl.Normalization = 'probability';
% hregCtrl.BinWidth = binsize;
% 
% hnewArchT.Normalization = 'probability';
% hnewArchT.BinWidth = binsize;

xlabel('Difference')
ylabel('Occurence')
legend({'Syldiff Ctrl','Syldiff ArchT'},'Location','northwest')

xlim([-1 8])
ylim([0 5])
axis square