%Give input of excel sheet and graph labels 
disp('Please select the file containing your ArchT NEW vs REGULAR CROSS SIMILARITY data')
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
disp('Please select the file containing your ArchT NEW vs PB CROSS similarity data')
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

%plot distribution of similarity for xvsx and xvsn comparisons
%run optoArchT_plotsimilarityselfandnewsyl.m before!


figure
hnewreg=histogram(Tself.x_Similarity,'edgecolor','r', 'DisplayStyle', 'stairs', 'LineWidth',2, 'EdgeAlpha', 0.5);
hold on
hnewpb=histogram(Tcross.x_Similarity,'edgecolor','g', 'DisplayStyle', 'stairs', 'LineWidth',2, 'EdgeAlpha', 0.5);
binsize=0.5;

hnewreg.Normalization = 'probability';
hnewreg.BinWidth = binsize;

hnewpb.Normalization = 'probability';
hnewpb.BinWidth = binsize;

xlabel('Similarity')
ylabel('Rel. frequency')
legend({'Similarity ArchT new vs reg','Similarity ArchT new vs PB'},'Location','northwest')

axis square

xlim([0 105]);
ylim([0 0.5]);
hold off

%plot cdf
vecArchTreg=Tself.x_Similarity';
vecArchTpb=Tcross.x_Similarity';

estimated_dist_ArchTreg = fitdist(vecArchTreg', 'Normal');
estimated_dist_ArchTPB = fitdist(vecArchTpb', 'Normal');

xRange_archtreg = linspace(min(vecArchTreg), max(vecArchTreg));
xRange_archtpb = linspace(min(vecArchTpb), max(vecArchTpb));

xPdf_archtreg = pdf(estimated_dist_ArchTreg, xRange_archtreg);
xPdf_archtpb = pdf(estimated_dist_ArchTPB, xRange_archtpb);

xCdf_archtreg = cdf(estimated_dist_ArchTreg, xRange_archtreg);
xCdf_archtpb = cdf(estimated_dist_ArchTPB, xRange_archtpb);

figure
hold on;
xlabel('Similarity [%]');
plot(xRange_archtreg, xCdf_archtreg, 'LineWidth', 2,'color','r');
plot(xRange_archtpb, xCdf_archtpb, 'LineWidth', 2,'color','g');
ylabel('Cumulative Density Function');

xlim([0 105]);
%ylim([0 0.9]);

legend({'ArchT new vs reg', 'ArchT new vs PB'},'Location','northwest')

hold off

%plot empirical cumulative density functions

[f_archtreg,x_archtreg,flo_archtreg,fup_archtreg]=ecdf(Tself.x_Similarity');
[f_archtpb,x_archtpb,flo_archtpb,fup_archtpb]=ecdf(Tcross.x_Similarity');

figure
hold on;
xlabel('Similarity [%]');
plot(x_archtreg, f_archtreg, 'LineWidth', 2,'color','r');
plot(x_archtreg, flo_archtreg, 'LineWidth', 0.5,'color','r');
plot(x_archtreg, fup_archtreg, 'LineWidth', 0.5,'color','r');

plot(x_archtpb, f_archtpb, 'LineWidth', 2,'color','g');
plot(x_archtpb, flo_archtpb, 'LineWidth', 0.5,'color','g');
plot(x_archtpb, fup_archtpb, 'LineWidth', 0.5,'color','g');

ylabel('Empirical Cumulative Density Function');

xlim([0 100]);
ylim([0 1]);

legend({'Similarity ArchTreg vs ArchTreg','','','Similarity ArchT reg vs PB'},'Location','northwest')

hold off





