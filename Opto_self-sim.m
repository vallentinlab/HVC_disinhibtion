%Give input of excel sheet and graph labels 
disp('Please select the file containing your regular ArchT SELF SIMILARITY data')
[file,dir] = uigetfile('*.xlsx');
filename = [dir file];
Tself = readtable(filename);
Tself = rmmissing(Tself); %get rid of useless data/missing entries

Tself.Properties.VariableNames([1 2]) = {'template' 'copy'};

templatesylsself = unique(Tself.template);
copysylsself = unique(Tself.copy);

%calculate overall similarity
%Tself.similarityCORR = ((Tself.x_Similarity/100).*(Tself.Accuracy/100)).*0.76;
%Tself.similarityCORR = Tself.x_Similarity

%remove zero entries
fehler = find(Tself.x_Similarity==0);
Tself(fehler,:) = [];

datsize = size(Tself);
%%
disp('Please select the file containing your regular Ctrl self similarity data')
[file,dir] = uigetfile('*.xlsx');
filename = [dir file];
Tcross = readtable(filename);
Tcross = rmmissing(Tcross); %get rid of useless data/missing entries

Tcross.Properties.VariableNames([1 2]) = {'template' 'copy'};

templatesylscross = unique(Tcross.template);
copysylscross = unique(Tcross.copy);

%calculate overall similarity
%Tcross.similarityCORR = ((Tcross.x_Similarity/100).*(Tcross.Accuracy/100)).*0.76;
%Tself.similarityCORR = Tself.x_Similarity

%remove zero entries
fehler = find(Tcross.x_Similarity==0);
Tcross(fehler,:) = [];

datsize = size(Tcross);

%%
disp('Please select the file containing your NEW ArchT SIMILARITY data')
[file,dir] = uigetfile('*.xlsx');
filename = [dir file];
Tpb = readtable(filename);
Tpb = rmmissing(Tpb); %get rid of useless data/missing entries

Tpb.Properties.VariableNames([1 2]) = {'template' 'copy'};

templatesylspb = unique(Tpb.template);
copysylspb = unique(Tpb.copy);

%calculate overall similarity
%Tpb.similarityCORR = ((Tpb.x_Similarity/100).*(Tpb.Accuracy/100)).*0.76;
%Tself.similarityCORR = Tself.x_Similarity

%remove zero entries
fehler = find(Tpb.x_Similarity==0);
Tpb(fehler,:) = [];

datsize = size(Tpb);

%%

%plot distribution of similarity for xvsx and xvsn comparisons
%run optoArchT_plotsimilarityselfandnewsyl.m before!


figure
hregArchT=histogram(Tself.x_Similarity,'edgecolor','r', 'DisplayStyle', 'stairs', 'LineWidth',2, 'EdgeAlpha', 0.5);
hold on
hregCtrl=histogram(Tcross.x_Similarity,'edgecolor','g', 'DisplayStyle', 'stairs', 'LineWidth',2, 'EdgeAlpha', 0.5);
hnewArchT=histogram(Tpb.x_Similarity,'edgecolor','k', 'DisplayStyle', 'stairs', 'LineWidth',2, 'EdgeAlpha', 0.5);

binsize=0.5;

hregArchT.Normalization = 'probability';
hregArchT.BinWidth = binsize;

hregCtrl.Normalization = 'probability';
hregCtrl.BinWidth = binsize;

hnewArchT.Normalization = 'probability';
hnewArchT.BinWidth = binsize;

xlabel('Similarity')
ylabel('Rel. frequency')
legend({'Similarity ArchT birds reg syls','Similarity Ctrl birds', 'Similarity ArchT birds new syls'},'Location','northwest')

axis square

xlim([0 105]);
ylim([0 0.5]);
hold off

%plot cumulative density functions

vecArchTreg=Tself.x_Similarity';
vecCtrlreg=Tcross.x_Similarity';
vecArchTnew=Tpb.x_Similarity';

estimated_dist_ArchTreg = fitdist(vecArchTreg', 'Normal');
estimated_dist_Ctrlreg = fitdist(vecCtrlreg', 'Normal');
estimated_dist_ArchTnew = fitdist(vecArchTnew', 'Normal');

xRange_archtreg = linspace(min(vecArchTreg), max(vecArchTreg));
xRange_ctrlreg = linspace(min(vecCtrlreg), max(vecCtrlreg));
xRange_archtnew = linspace(min(vecArchTnew), max(vecArchTnew));

xPdf_archtreg = pdf(estimated_dist_ArchTreg, xRange_archtreg);
xPdf_ctrlreg = pdf(estimated_dist_Ctrlreg, xRange_ctrlreg);
xPdf_archtnew = pdf(estimated_dist_ArchTnew, xRange_archtnew);

xCdf_archtreg = cdf(estimated_dist_ArchTreg, xRange_archtreg);
xCdf_ctrlreg = cdf(estimated_dist_Ctrlreg, xRange_ctrlreg);
xCdf_archtnew = cdf(estimated_dist_ArchTnew, xRange_archtnew);

figure
hold on;
xlabel('Similarity [%]');
plot(xRange_archtreg, xCdf_archtreg, 'LineWidth', 2,'color','r');
plot(xRange_ctrlreg, xCdf_ctrlreg, 'LineWidth', 2,'color','g');
plot(xRange_archtnew, xCdf_archtnew, 'LineWidth', 2,'color','k');
ylabel('Cumulative Density Function');

xlim([0 105]);
ylim([0 1]);

legend({'Similarity ArchT birds reg syls','Similarity Ctrl birds', 'Similarity ArchT birds new syls'},'Location','northwest')

hold off

%plot empirical cumulative density functions

[f_archtreg,x_archtreg,flo_archtreg,fup_archtreg]=ecdf(Tself.x_Similarity');
[f_ctrlreg,x_ctrlreg,flo_ctrlreg,fup_ctrlreg]=ecdf(Tcross.x_Similarity');
[f_archtnew,x_archtnew,flo_archtnew,fup_archtnew]=ecdf(Tpb.x_Similarity');

figure
hold on;
xlabel('Similarity [%]');
plot(x_archtreg, f_archtreg, 'LineWidth', 2,'color','r');
plot(x_archtreg, flo_archtreg, 'LineWidth', 0.5,'color','r');
plot(x_archtreg, fup_archtreg, 'LineWidth', 0.5,'color','r');

plot(x_ctrlreg, f_ctrlreg, 'LineWidth', 2,'color','g');
plot(x_ctrlreg, flo_ctrlreg, 'LineWidth', 0.5,'color','g');
plot(x_ctrlreg, fup_ctrlreg, 'LineWidth', 0.5,'color','g');

plot(x_archtnew, f_archtnew, 'LineWidth', 2,'color','k');
plot(x_archtnew, flo_archtnew, 'LineWidth', 0.5,'color','k');
plot(x_archtnew, fup_archtnew, 'LineWidth', 0.5,'color','k');

ylabel('Empirical Cumulative Density Function');

xlim([0 100]);
ylim([0 1]);

legend({'Similarity ArchT birds reg syls','','','Similarity Ctrl birds','','', 'Similarity ArchT birds new syls'},'Location','northwest')

hold off

