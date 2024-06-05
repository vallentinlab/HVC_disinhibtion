%run statistics for similarity comparisons ArchT birds

%load file with data first
%self and cross separately

datacross=table2array(ALLsimilarityforTest(:,1:4));
dataself=table2array(ALLsimilarityforTest(:,5:7));

%self comparison
varnames = ALLsimilarityforTest.Properties.VariableNames;
[p,tbl,stats] = kruskalwallis(dataself,[],'off');
c=multcompare(stats);
tblmult = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%cross comparison
varnames = ALLsimilarityforTest.Properties.VariableNames;
[p,tbl,stats] = kruskalwallis(datacross,[],'off');
c=multcompare(stats);
tblmult = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%alldata
data_complete=table2array(ALLsimilarityforTest(:,1:7));
varnames = ALLsimilarityforTest.Properties.VariableNames;
[p,tbl,stats] = kruskalwallis(data_complete,[],'off');
c=multcompare(stats);
tblmult = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%comparisons of all cross to ArchT reg self
data_vsarchtreg=table2array(ALLsimilarityforTest(:,[1 2 3 4 6]));
varnames = ALLsimilarityforTest.Properties.VariableNames;
[p,tbl,stats] = kruskalwallis(data_vsarchtreg,[],'off');
c=multcompare(stats);
tblmult = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

