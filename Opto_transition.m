%function[TPnew,h,g1]=TP_calc_perbout(filename)
%this code Transition Probability per bout
clear all
close all

%Give input of excel sheet and graph labels 
disp('Please select the file containing your exported syllable data (NOTHING ELSE!)')
[file,dir] = uigetfile('*.xlsx');
filename= [dir file];

T = readtable(filename);
B=rmmissing(T); %get rid of data of unnamed notes
%Separate all files first
Notes=length(B.Label);
FileIndex=unique(B.File_);
sylls=unique(B.Label);
labels=sylls';
A={};
for m=1:length(FileIndex)
    p=1;
    for z=1:length(B.File_)
        if B.File_(z)== FileIndex(m)
            A{m}(p,1:5)=[B.File_(z),B.Syllable_(z),B.Start_s_(z),B.End_s_(z),B.Label(z)];
            p=p+1;
        end
        
    end
end
for pp=1:length(FileIndex)

test=A{1,pp}(:,5);
seqs{pp}=[test{:}];
end
%seqs contains sequence per bout

seqall=[seqs{:}];
maxnumcol=max(double(unique(seqall)))-96;

%make_regexp_matrix per bout
%then make TP matrix per bout
for q=1:length(seqs)
    seqmain=seqs{q}
    for i=1:26
        for j=1:26
            testcell{i,j}=[char(i+96), char(j+96)];
            countcell{i,j}=length(regexp(seqmain,testcell{i,j}));
        end
    end
countmat=cell2mat(countcell);
TPcount(1:maxnumcol,1:maxnumcol)=countmat(1:maxnumcol,1:maxnumcol);
tot=sum(nonzeros(TPcount));
TPnum{q}=(TPcount*100)/tot;
end

%average TP over all files
test=cat(3, TPnum{:}); %convert TPnum cell into a 3D matric
avgTP= mean(test,3); %averaging over the third dimension

%run this if the syllables are not labelled in sequence
%otherwise adjacency matrix and number of nodes do not match
%TPnew=delete_emptycols(avgTP);

%changed TPnew to avgTP
g1 = digraph(avgTP,labels);
h=plot(g1,'EdgeLabel',round(g1.Edges.Weight,3))
set(gcf,'Position',[100 100 900 700]);
%layout(h,'layered',"Direction","right","Sources",1,"Sinks",2);
layout(h,'circle','Center','c');
title(filename);
colormap jet           % select color palette 
colorbar
h.EdgeCData=g1.Edges.Weight;    % define edge colors
 h.NodeColor='k';
 h.NodeFontSize=11;
 h.EdgeFontSize=11;

%plot the transitionprobabilities in a matrixplot

labels={'a','b','c','d','e','f','g','h','i','j','k'}; %new labels according to max labels in light condition

lowestValue = min(avgTP(avgTP(:)>0));
highestValue = max(avgTP(:));

figure
imagesc(avgTP)
xticks([1:15])
yticks([1:15])
set(gca,'xticklabel',labels)
set(gca,'yticklabel',labels)
%cmap=parula(2048);
%cmap=cbrewer('seq','YlGnBu',2048);
cmap = cool(2048)
colormap(cmap);
caxis(gca,[lowestValue-0.02, highestValue]);
% Make less than lowest value white:
cmap(1,:)=[1,1,1];
colormap(cmap)
colorbar
axis square
title(file)

