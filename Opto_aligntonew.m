%align motifs w new syls

%load file ArchT_newsyllabelsforalignment_Merge

%variables
    %ID bird ID
    %sylstate 0 = regular / 1 = new
    %on syllable onset
    %off syllable offset
    %sylid given syllable ID
    %motifid running motif number for each bird
    %phase pb = 4wks playback /post = 4wks after playbacks
    
%to do: extract data per phase, per bird, per motif, add syl length, plot,
%colourcode

%to align: for each given motif just subtract timepoint of first syllable's
%start
clc
%%%%%%%%%%%%%
%add syllable duration
ArchTnewsyllabelsforalignmentMERGE.dur=ArchTnewsyllabelsforalignmentMERGE.off-ArchTnewsyllabelsforalignmentMERGE.on;

%split data in pb and post
pbdata = ArchTnewsyllabelsforalignmentMERGE.phase=='pb';
pbdata = ArchTnewsyllabelsforalignmentMERGE(pbdata,:);

postdata = ArchTnewsyllabelsforalignmentMERGE.phase=='post';
postdata = ArchTnewsyllabelsforalignmentMERGE(postdata,:);

%get bird IDs
birdids=unique(postdata.ID);

%build matrix with each motif


%for plotIDx
plotno = 1;
sepmotifs={};
for i = 1:size((birdids),1)
    %extract individual birds
    currbird_pb=(pbdata.ID==birdids(i));
    currbird_pb=(pbdata(currbird_pb,:));
    
    currbird_post=(postdata.ID==birdids(i));
    currbird_post=(postdata(currbird_post,:));
    
    %extract individual motifs
    %pb
    
    %for motifIDx
    l = 1;
    m = 1;
    
    motifidx_pb=unique(currbird_pb.motifid);
    for j = 1:size((motifidx_pb),1)
        currmotif_pb=(currbird_pb.motifid==motifidx_pb(j));
        currmotif_pb=(currbird_pb(currmotif_pb,:));
        
        %dump in array
        sepmotifs{1,l}=currmotif_pb;
        l=l+1;
    end
    %post
    motifidx_post=unique(currbird_post.motifid);
    for k = 1:size((motifidx_post),1)
        currmotif_post=(currbird_post.motifid==motifidx_post(k));
        currmotif_post=(currbird_post(currmotif_post,:));
        
        %dump in array
        sepmotifs{2,m}=currmotif_post;
        m=m+1;
    end
    
        %get number of motifs in current dataset
        currmotifcount=sum(~cellfun(@isempty,sepmotifs),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot pb motifs

        subplot(length(birdids),2,plotno)
        titlename=string(birdids(i));
        title("PB motifs " + titlename);
        %ylim([0 85])
        %xlim([-4.5 3])
        hold on;
        %lineplots of data 
        y=ones(1,2);
        z=ones(1,2);
        for n = 1:currmotifcount(1)
            o = 1:size(sepmotifs{1,n},1);
            currline = (table2array(sepmotifs{1,n}(o,3:4)));
            
            %align to 0
                %alignfactor=table2array(sepmotifs{1,n}(1,3));
            
            %align to first new syl
            idnewsylpb=cellfun(@(x)isequal(x,1),table2cell(sepmotifs{1,n}(:,2)));
            idfirstnewsylpb=idnewsylpb(1:find(idnewsylpb == 1,1));
            alignfactorpb=table2array(sepmotifs{1,n}(idfirstnewsylpb,3));
            
            currline=currline-alignfactorpb;
            
            %plot(currline,y);
            
            %plot colour coded syls
            idxnew=find(currbird_pb.sylstate==1);
            idxreg=find(currbird_pb.sylstate==0);

            regsylids=unique(currbird_pb.sylid(idxreg));
            newsylids=unique(currbird_pb.sylid(idxnew));

            %create blueish color-array for reg syls
            lightblue=[204/256 204/256 1];
            darkblue=[0 0 153/256];
            cmap_reg = interp1([0, 1], [lightblue; darkblue], linspace(0, 1, length(regsylids)));

            %create redish  color-array for new syls
            lightred=[255/256 153/256 153/256];
            darkred=[153/256 1/256 1/256];
            cmap_new= interp1([0, 1], [lightred; darkred], linspace(0, 1, length(newsylids)));
            
            
            %actual plots
            for p = 1:length(regsylids)
                currsylid=regsylids(p);
                colflag=find(sepmotifs{1,n}.sylid==currsylid);
                for q=1:length(colflag)
                    currsyl=currline(colflag(q),:);
                    plot(currsyl,y,'Color',cmap_reg(p,:),'Linewidth',2.5);
                end
            end

            for r = 1:length(newsylids)
                currsylid=newsylids(r);
                colflag=find(sepmotifs{1,n}.sylid==currsylid);
                for s=1:length(colflag)
                    currsyl=currline(colflag(s),:);
                    plot(currsyl,y,'Color',cmap_new(r,:),'Linewidth',2.5);
                end
            end
            
            y=y+1;
        end
        
        %ylim([0 currmotifcount(1)+3])
        xlim([-4.5 3])
        %axis square
        
        hold off
        plotno=plotno+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plot post motifs
        subplot(length(birdids),2,plotno)
        title("Post PB motifs " + titlename);
        
        hold on;
        %lineplots of data 
        %y=ones(1,2);
        for n = 1:currmotifcount(2)
            %o = 1:size(sepmotifs{2,n},1);
            currline = table2array(sepmotifs{2,n}(:,3:4));
            %currline = table2array(currmotif_post(:,3:4));
            
            %align to 0
                %alignfactor=table2array(sepmotifs{2,n}(1,3));
                
            %align to first new syl
            idnewsylpost=cellfun(@(x)isequal(x,1),table2cell(sepmotifs{2,n}(:,2)));
            idfirstnewsylpost=idnewsylpost(1:find(idnewsylpost == 1,1));
            alignfactorpost=table2array(sepmotifs{2,n}(idfirstnewsylpost,3));
            
            
            currline=currline-alignfactorpost;
            
            %plot(currline,y);
            
            %plot syllable lines with colours according to identiy reg/new
                        
            %plot colour coded syls
            idxnew=find(currbird_post.sylstate==1);
            idxreg=find(currbird_post.sylstate==0);

            regsylids=unique(currbird_post.sylid(idxreg));
            newsylids=unique(currbird_post.sylid(idxnew));

            %create blueish color-array for reg syls
            lightblue=[204/256 204/256 1];
            darkblue=[0 0 153/256];
            cmap_reg = interp1([0, 1], [lightblue; darkblue], linspace(0, 1, length(regsylids)));

            %create redish  color-array for new syls
            lightred=[255/256 153/256 153/256];
            darkred=[153/256 1/256 1/256];
            cmap_new= interp1([0, 1], [lightred; darkred], linspace(0, 1, length(newsylids)));

            %actual plots
           
            for p = 1:length(regsylids)
                currsylid=regsylids(p);
                colflag=find(sepmotifs{2,n}.sylid==currsylid);
                for q=1:length(colflag)
                    currsyl=currline(colflag(q),:);
                    plot(currsyl,z,'Color',cmap_reg(p,:),'Linewidth',2.5);
                end
            end

            for r = 1:length(newsylids)
                currsylid=newsylids(r);
                colflag=find(sepmotifs{2,n}.sylid==currsylid);
                for s=1:length(colflag)
                    currsyl=currline(colflag(s),:);
                    plot(currsyl,z,'Color',cmap_new(r,:),'Linewidth',2.5);
                end
            end
            
            
            z=z+1;
        
        end
        
      
        xlim([-4.5 3])
        %ylim([0 z(1)+3])
        %axis square
        
        hold off
        plotno=plotno+1;
        
        %flush motifvector
        clear sepmotifs

end