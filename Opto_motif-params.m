%plot consistency, linearity, stereotypy

%select excel file
disp('Please select the file with raw data')
[file,dir] = uigetfile('*.xlsx');

%merge directory and file
filestereo = [dir file];

%import excel file separated into numbers, text and the raw file
t_raw= readtable(filestereo);
t_raw.treatment = categorical(t_raw.treatment);
t_raw.status = categorical(t_raw.status);

%split data for plotting
ctrl=t_raw(ismember(t_raw.treatment,'ctrl'),:);
archt=t_raw(ismember(t_raw.treatment,'archt'),:);
archtID=unique(archt.birdID);
ctrlID=unique(ctrl.birdID);
%%
%plot linearity
tiledlayout(3,2);
nexttile(1,[1 1])
%subplot(3,2,1,'align')
axis square
title('Song linearity')
xticks([1 2 3])
xticklabels({'Pre PB','4wks PB','4wks post PB'})
xlim([0.5 3.5])
ylabel('Linearity score')
ylim([0.1 0.8])

for i=1:length(archtID)
    hold on
    plot(archt.linearity(ismember(archt.birdID,archtID(i))),...
    '-o',...
    'Color',[1 0.81 0 0.2],...
    'LineWidth', 2,...
    'MarkerSize',5)
end

prearch=mean(archt.linearity(ismember(archt.status,'pre')));
pbarch=mean(archt.linearity(ismember(archt.status,'4wkpb')));
postarch=mean(archt.linearity(ismember(archt.status,'4wkpost')));
archmeans=[prearch pbarch postarch];
plot(archmeans,...
    '-o',...
    'Color',[1 0.81 0 1],...
    'LineWidth', 4,...
    'MarkerSize',5,...
    'MarkerEdgeColor',[1 0.81 0])

for i=1:length(ctrlID)
    hold on
    plot(ctrl.linearity(ismember(ctrl.birdID,ctrlID(i))),...
    '-o',...
    'Color',[0.5 0.5 0.5 0.2],...
    'LineWidth', 2,...
    'MarkerSize',5)
end

prectrl=mean(ctrl.linearity(ismember(ctrl.status,'pre')));
pbctrl=mean(ctrl.linearity(ismember(ctrl.status,'4wkpb')));
postctrl=mean(ctrl.linearity(ismember(ctrl.status,'4wkpost')));
ctrlmeans=[prectrl pbctrl postctrl];
plot(ctrlmeans,...
    '-o',...
    'Color',[0.5 0.5 0.5 1],...
    'LineWidth', 4,...
    'MarkerSize',5)

hold off

%%
%plot consistency
%subplot(3,2,3,'align')
nexttile(3,[1 1])
axis square
title('Song consistency')
xticks([1 2 3])
xticklabels({'Pre PB','4wks PB','4wks post PB'})
xlim([0.5 3.5])
ylabel('Consistency score')
ylim([0.1 0.8])

for i=1:length(archtID)
    hold on
    plot(archt.consistency(ismember(archt.birdID,archtID(i))),...
    '-o',...
    'Color',[1 0.81 0 0.2],...
    'LineWidth', 2,...
    'MarkerSize',5)
end

prearch=mean(archt.consistency(ismember(archt.status,'pre')));
pbarch=mean(archt.consistency(ismember(archt.status,'4wkpb')));
postarch=mean(archt.consistency(ismember(archt.status,'4wkpost')));
archmeans=[prearch pbarch postarch];
plot(archmeans,...
    '-o',...
    'Color',[1 0.81 0 1],...
    'LineWidth', 4,...
    'MarkerSize',5,...
    'MarkerEdgeColor',[1 0.81 0])

for i=1:length(ctrlID)
    hold on
    plot(ctrl.consistency(ismember(ctrl.birdID,ctrlID(i))),...
    '-o',...
    'Color',[0.5 0.5 0.5 0.2],...
    'LineWidth', 2,...
    'MarkerSize',5)
end

prectrl=mean(ctrl.consistency(ismember(ctrl.status,'pre')));
pbctrl=mean(ctrl.consistency(ismember(ctrl.status,'4wkpb')));
postctrl=mean(ctrl.consistency(ismember(ctrl.status,'4wkpost')));
ctrlmeans=[prectrl pbctrl postctrl];
plot(ctrlmeans,...
    '-o',...
    'Color',[0.5 0.5 0.5 1],...
    'LineWidth', 4,...
    'MarkerSize',5)

hold off

%%
%plot stereotypy
%subplot(3,2,5,'align')
nexttile(5,[1 1])
axis square
title('Song stereotypy')
xticks([1 2 3])
xticklabels({'Pre PB','4wks PB','4wks post PB'})
xlim([0.5 3.5])
ylabel('Stereotypy score')
ylim([0.1 0.8])

for i=1:length(archtID)
    hold on
    plot(archt.stereotypy(ismember(archt.birdID,archtID(i))),...
    '-o',...
    'Color',[1 0.81 0 0.2],...
    'LineWidth', 2,...
    'MarkerSize',5)
end

prearch=mean(archt.stereotypy(ismember(archt.status,'pre')));
pbarch=mean(archt.stereotypy(ismember(archt.status,'4wkpb')));
postarch=mean(archt.stereotypy(ismember(archt.status,'4wkpost')));
archmeans=[prearch pbarch postarch];
plot(archmeans,...
    '-o',...
    'Color',[1 0.81 0 1],...
    'LineWidth', 4,...
    'MarkerSize',5,...
    'MarkerEdgeColor',[1 0.81 0])

for i=1:length(ctrlID)
    hold on
    plot(ctrl.stereotypy(ismember(ctrl.birdID,ctrlID(i))),...
    '-o',...
    'Color',[0.5 0.5 0.5 0.2],...
    'LineWidth', 2,...
    'MarkerSize',5)
end

prectrl=mean(ctrl.stereotypy(ismember(ctrl.status,'pre')));
pbctrl=mean(ctrl.stereotypy(ismember(ctrl.status,'4wkpb')));
postctrl=mean(ctrl.stereotypy(ismember(ctrl.status,'4wkpost')));
ctrlmeans=[prectrl pbctrl postctrl];
plot(ctrlmeans,...
    '-o',...
    'Color',[0.5 0.5 0.5 1],...
    'LineWidth', 4,...
    'MarkerSize',5)

hold off
%%
%load consistency, linearity, stereotypy DELTA

%select excel file
disp('Please select the file with deltas')
[file,dir] = uigetfile('*.xlsx');

%merge directory and file
filestereo = [dir file];

%import excel file separated into numbers, text and the raw file
t_delta= readtable(filestereo) ;

t_delta.treatment = categorical(t_delta.treatment);
t_delta.status = categorical(t_delta.status);

%split data for plotting
ctrl=t_delta(ismember(t_delta.treatment,'ctrl'),:);
ctrlpb=ctrl(ismember(ctrl.status,'pretopb'),:);
ctrlpost=ctrl(ismember(ctrl.status,'pbtopost'),:);
archt=t_delta(ismember(t_delta.treatment,'archt'),:);
archtpb=archt(ismember(archt.status,'pretopb'),:);
archtpost=archt(ismember(archt.status,'pbtopost'),:);
%%
%plot delta 4wkpbdata nxt to postdata, 1 boxplot for each variable
%linearity

archtcol=[1 0.81 0];
ctrlcol=[0.5 0.5 0.5];

nexttile(2,[1 1])
%subplot(3,2,2,'align')

x=ones(length(ctrlpb.linearity),1);
x2=ones(length(archtpb.linearity),1)*2;
x3=ones(length(ctrlpost.linearity),1)*3;
x4=ones(length(archtpost.linearity),1)*4;

boxplot([ctrlpb.linearity;archtpb.linearity;ctrlpost.linearity;archtpost.linearity],[x;x2;x3;x4],...
    'Widths',0.3,...    
    'Whisker',inf,...
    'Colors',[(ctrlcol);(archtcol)])
    %'labels',{'ctrl pre to 4wk pb','archt pre to 4wk pb','ctrl 4wk pb to 4wk post','archt 4wkpb to 4wk post'})
title('Linearity delta')
ylabel('Rel. change')
ylim([-0.8 0.6])
xticks([1.5 3.5])
xticklabels({'Pre to PB','PB to post'})
hold on
scatter(x,ctrlpb.linearity,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5]);
scatter(x2,archtpb.linearity,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[1 0.81 0],...
    'MarkerEdgeColor',[1 0.81 0]);
scatter(x3,ctrlpost.linearity,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5]);
scatter(x4,archtpost.linearity,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[1 0.81 0],...
    'MarkerEdgeColor',[1 0.81 0]);
yline(0,'--')
axis square
box off
hold off

%consistency

nexttile(4,[1 1])
%subplot(3,2,4,'align')

x=ones(length(ctrlpb.consistency),1);
x2=ones(length(archtpb.consistency),1)*2;
x3=ones(length(ctrlpost.consistency),1)*3;
x4=ones(length(archtpost.consistency),1)*4;

boxplot([ctrlpb.consistency;archtpb.consistency;ctrlpost.consistency;archtpost.consistency],[x;x2;x3;x4],...
    'Whisker',inf,...
    'Widths',0.3,...
    'Colors',[(ctrlcol);(archtcol)])
    %'labels',{'ctrl pre to 4wk pb','archt pre to 4wk pb','ctrl 4wk pb to 4wk post','archt 4wkpb to 4wk post'})
title('Consistency delta')
ylabel('Rel. change')
ylim([-0.8 0.6])
xticks([1.5 3.5])
xticklabels({'Pre to PB','PB to post'})
hold on
scatter(x,ctrlpb.consistency,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5]);
scatter(x2,archtpb.consistency,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[1 0.81 0],...
    'MarkerEdgeColor',[1 0.81 0]);
scatter(x3,ctrlpost.consistency,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5]);
scatter(x4,archtpost.consistency,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[1 0.81 0],...
    'MarkerEdgeColor',[1 0.81 0]);
yline(0,'--')
axis square
box off
hold off

%stereotypy

nexttile(6,[1 1])
%subplot(3,2,6,'align')

x=ones(length(ctrlpb.stereotypy),1);
x2=ones(length(archtpb.stereotypy),1)*2;
x3=ones(length(ctrlpost.stereotypy),1)*3;
x4=ones(length(archtpost.stereotypy),1)*4;

boxplot([ctrlpb.stereotypy;archtpb.stereotypy;ctrlpost.stereotypy;archtpost.stereotypy],[x;x2;x3;x4],...
    'Widths',0.3,...    
    'Whisker',inf,...
    'Colors',[(ctrlcol);(archtcol)])
    %'labels',{'ctrl pre to 4wk pb','archt pre to 4wk pb','ctrl 4wk pb to 4wk post','archt 4wkpb to 4wk post'})
title('Stereotypy delta')
ylabel('Rel. change')
ylim([-0.8 0.6])
xticks([1.5 3.5])
xticklabels({'Pre to PB','PB to post'})
hold on
scatter(x,ctrlpb.stereotypy,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5]);
scatter(x2,archtpb.stereotypy,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[1 0.81 0],...
    'MarkerEdgeColor',[1 0.81 0]);
scatter(x3,ctrlpost.stereotypy,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5]);
scatter(x4,archtpost.stereotypy,...
    20,...
    'jitter','on', 'jitterAmount',0.1,...
    'MarkerFaceColor',[1 0.81 0],...
    'MarkerEdgeColor',[1 0.81 0]);
yline(0,'--')
axis square
box off
hold off

%%
figure
%test development
%linearity dependent on treatment and status
g1=t_raw.treatment;
g2=t_raw.birdID;
g3=t_raw.status;
y=t_raw.linearity;
[~,~,statslin]=anovan(y,{g1 g3},...
    'model','interaction','varnames',{'treatment','status'})
[resultslin,~,~,gnameslin] = multcompare(statslin,"Dimension",[1 2]);
tbllin = array2table(resultslin,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbllin.("Group A")=gnameslin(tbllin.("Group A"));
tbllin.("Group B")=gnameslin(tbllin.("Group B"));

%consistency dependent on treatment and status
y=t_raw.consistency;
[~,~,statscon]=anovan(y,{g1 g3},...
    'model','interaction','varnames',{'treatment','status'})
[resultscon,~,~,gnamescon] = multcompare(statscon,"Dimension",[1 2]);
tblcon = array2table(resultscon,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblcon.("Group A")=gnamescon(tblcon.("Group A"));
tblcon.("Group B")=gnamescon(tblcon.("Group B"));

%stereotypy dependent on treatment and status
y=t_raw.stereotypy;
[~,~,statsster]=anovan(y,{g1 g3},...
    'model','interaction','varnames',{'treatment','status'})
[resultsster,~,~,gnamesster] = multcompare(statscon,"Dimension",[1 2]);
tblster = array2table(resultsster,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblster.("Group A")=gnamesster(tblster.("Group A"));
tblster.("Group B")=gnamesster(tblster.("Group B"));

%%
%test delta
%linearity dependent on treatment and status
g1=t_delta.treatment;
g2=t_delta.birdID;
g3=t_delta.status;
y=t_delta.linearity;
[~,~,statslindel]=anovan(y,{g1 g3},...
    'model','interaction','varnames',{'treatment','status'})
[resultslindel,~,~,gnameslindel] = multcompare(statslindel,"Dimension",[1 2]);
tbllindel = array2table(resultslindel,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbllindel.("Group A")=gnameslindel(tbllindel.("Group A"));
tbllindel.("Group B")=gnameslindel(tbllindel.("Group B"));

%consistency dependent on treatment and status
y=t_delta.consistency;
[~,~,statscondel]=anovan(y,{g1 g3},...
    'model','interaction','varnames',{'treatment','status'})
[resultscondel,~,~,gnamescondel] = multcompare(statscondel,"Dimension",[1 2]);
tblcondel = array2table(resultscondel,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblcondel.("Group A")=gnamescondel(tblcondel.("Group A"));
tblcondel.("Group B")=gnamescondel(tblcondel.("Group B"));

%stereotypy dependent on treatment and status
y=t_delta.stereotypy;
[~,~,statssterdel]=anovan(y,{g1 g3},...
    'model','interaction','varnames',{'treatment','status'})
[resultssterdel,~,~,gnamessterdel] = multcompare(statscondel,"Dimension",[1 2]);
tblsterdel = array2table(resultssterdel,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tblsterdel.("Group A")=gnamessterdel(tblsterdel.("Group A"));
tblsterdel.("Group B")=gnamessterdel(tblsterdel.("Group B"));