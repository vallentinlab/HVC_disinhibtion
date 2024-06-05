%extract SAP comparison information from excel sheets
clear all
close all
%Give input of excel sheet and graph labels 
disp('Please select the file containing your SELF SIMILARITY data')
[file,dir] = uigetfile('*.xlsx');
filename = [dir file];
Tself = readtable(filename);
Tself = rmmissing(Tself); %get rid of useless data/missing entries

Tself.Properties.VariableNames([1 2]) = {'template' 'copy'};

templatesyls = unique(Tself.template);
copysyls = unique(Tself.copy);

%%
%calculate overall similarity
%Tself.similarityCORR = Tself.x_Similarity

%remove zero entries
fehler = find(Tself.Accuracy==0);
Tself(fehler,:) = [];

datsize = size(Tself);
%%
%Extract self similarity data

j=1;
for i = 1:datsize(1);
    
    if startsWith(Tself.template(i),'a') && startsWith(Tself.copy(i),'a')
        similarity_a_a(j)=Tself.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1)
    
    if startsWith(Tself.template(i),'b') && startsWith(Tself.copy(i),'b')
        similarity_b_b(j)=Tself.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1)
    
   if startsWith(Tself.template(i),'c') && startsWith(Tself.copy(i),'c')
        similarity_c_c(j)=Tself.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1)
    
     if startsWith(Tself.template(i),'d') && startsWith(Tself.copy(i),'d')
        similarity_d_d(j)=Tself.Accuracy(i);
        j=j+1;
        
     end
end

j=1;
for i = 1:datsize(1)
    
     if startsWith(Tself.template(i),'n') && startsWith(Tself.copy(i),'n')
        similarity_n_n(j)=Tself.Accuracy(i);
        j=j+1;
        
    end
end
%%
%extract cross similarity
%Give input of excel sheet and graph labels 
disp('Please select the file containing your CROSS SIMILARITY data')
[file,dir] = uigetfile('*.xlsx');
filename = [dir file];
Tcross = readtable(filename);
Tcross = rmmissing(Tcross); %get rid of useless data/missing entries

Tcross.Properties.VariableNames([1 2]) = {'template' 'copy'};

templatesyls = unique(Tcross.template);
copysyls = unique(Tcross.copy);

%%
%calculate overall similarity
%Tcross.similarityCORR = ((Tcross.x_Similarity/100).*(Tcross.Accuracy/100)).*0.76;
%Tcross.similarityCORR = Tcross.x_Similarity;


%remove zero entries
fehler = find(Tcross.Accuracy==0);
Tcross(fehler,:) = [];

datsize = size(Tcross);
%%
%Extract cross similarity data

j=1;
for i = 1:datsize(1);
    
    if startsWith(Tcross.template(i),'n') && startsWith(Tcross.copy(i),'a')
        similarity_n_a(j)=Tcross.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1)
    
    if startsWith(Tcross.template(i),'n') && startsWith(Tcross.copy(i),'b')
        similarity_n_b(j)=Tcross.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1)
    
   if startsWith(Tcross.template(i),'n') && startsWith(Tcross.copy(i),'c')
        similarity_n_c(j)=Tcross.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1)
    
     if startsWith(Tcross.template(i),'n') && startsWith(Tcross.copy(i),'d')
        similarity_n_d(j)=Tcross.Accuracy(i);
        j=j+1;
        
     end
end

j=1;
for i = 1:datsize(1)
    
     if startsWith(Tcross.template(i),'n') && startsWith(Tcross.copy(i),'pba')
        similarity_n_pba(j)=Tcross.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1)
    
     if startsWith(Tcross.template(i),'n') && startsWith(Tcross.copy(i),'pbc')
        similarity_n_pbc(j)=Tcross.Accuracy(i);
        j=j+1;
        
    end
end

%%
%get self similarity of ABAB PB
%Give input of excel sheet and graph labels 
disp('Please select the file containing your PLAYBACK SELF SIMILARITY data')
[file,dir] = uigetfile('*.xlsx');
filename = [dir file];
Tselfpb = readtable(filename);
Tselfpb = rmmissing(Tselfpb); %get rid of useless data/missing entries

Tselfpb.Properties.VariableNames([1 2]) = {'template' 'copy'};

templatesyls = unique(Tselfpb.template);
copysyls = unique(Tselfpb.copy);

%%
%calculate overall similarity
%Tselfpb.similarityCORR = ((Tselfpb.x_Similarity/100).*(Tselfpb.Accuracy/100)).*0.76;
%Tself.similarityCORR = Tself.x_Similarity

%remove zero entries
fehler = find(Tselfpb.Accuracy==0);
Tselfpb(fehler,:) = [];

datsize = size(Tselfpb);
%%
%Extract self similarity data

j=1;
for i = 1:datsize(1);
    
    if startsWith(Tselfpb.template(i),'a') && startsWith(Tselfpb.copy(i),'a')
        similarity_pba_pba(j)=Tselfpb.Accuracy(i);
        j=j+1;
        
    end
end

j=1;
for i = 1:datsize(1);
    
    if startsWith(Tselfpb.template(i),'c') && startsWith(Tselfpb.copy(i),'c')
        similarity_pbc_pbc(j)=Tselfpb.Accuracy(i);
        j=j+1;
        
    end
end


%%
%plotting
%%
%merge variables

%selfsimilarities
selfsimilarity = {similarity_a_a';similarity_b_b';similarity_c_c';similarity_d_d';similarity_pba_pba';similarity_pbc_pbc';similarity_n_n'};

crosssimilarity = {similarity_n_a';similarity_n_b';similarity_n_c';similarity_n_d';similarity_n_pba';similarity_n_pbc'};

x1=1;
x2=2;
figure
hold

c=jet(length(selfsimilarity));

for i = 1:length(selfsimilarity)
    e1=errorbar(x1,median(selfsimilarity{i,1}),std(selfsimilarity{i,1})/sqrt(length(selfsimilarity{i,1})),'o')
    set(e1,'Color',c(i,:))
    set(e1,'MarkerEdgeColor',c(i,:))
end

for i = 1:length(crosssimilarity)
    e1=errorbar(x2,median(crosssimilarity{i,1}),std(crosssimilarity{i,1})/sqrt(length(crosssimilarity{i,1})),'o')
    set(e1,'Color',c(i,:))
    set(e1,'MarkerEdgeColor',c(i,:))
end

for i = 1:(length(selfsimilarity)-1)
    line([x1 x2],[median(selfsimilarity{i,1}) median(crosssimilarity{i,1})])
end


legend('sim_a_a','sim_b_b','sim_c_c','sim_d_d','sim_pba_pba','sim_pbc_pbc','sim_n_n','sim_n_a','sim_n_b','sim_n_c','sim_n_d','sim_n_pba','sim_n_pbc');
title('SAP Self- and Cross-Accuracy scores');
ylabel('Accuracy [%]');
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'Self accuracy','Cross accuracy'});
hold off
