%folder to process: all your data should be in folders with the name of the
%condition
rootfolder='Z:\Clustering_datasets\';


%list folders (and only folders)
cd(rootfolder);
listfol=dir;
listfol=listfol(3:end);
folderrank=[];
p=1;
for ifol=1:size(listfol,1)
    if listfol(ifol).isdir==1
        folderrank(p)=ifol;
        p=p+1;
    end
end
folderrank=folderrank';
listfol=listfol(folderrank);

    MMaster_summary=cell(1,9);
    MMaster_summary(1,1)={'folder'};
    MMaster_summary(1,3)={'nROI'};
    MMaster_summary(1,4)={'mean % colloc Ch1Ch2Ch3'}; 
    MMaster_summary(1,5)={'sem % colloc Ch1Ch2Ch3'}; 
    MMaster_summary(1,6)={'mean % colloc Ch1Ch2 unique'};
    MMaster_summary(1,7)={'sem % colloc Ch1Ch2 unique'};
    MMaster_summary(1,8)={'mean % colloc Ch1Ch3 unique'};   
    MMaster_summary(1,9)={'sem % colloc Ch1Ch3 unique'};
    ppp=2;


% 
 for ifol=1:size(listfol,1)  
    cd(cat(2,listfol(ifol).folder,'\',listfol(ifol).name));
    dirtif=dir("*_c1c2_summary.mat");
    
     Master_summary=cell(1,9);
     Master_summary(1,1)={'folder'};
     Master_summary(1,2)={'image'};  
     Master_summary(1,3)={'ROI'};
     Master_summary(1,4)={'nb spots Ch1'};
     Master_summary(1,5)={'nb spots colloc Ch1Ch2'};
     Master_summary(1,6)={'nb spots colloc Ch1Ch3'};
     Master_summary(1,7)={'nb spots colloc Ch1Ch2Ch3'}; 
     Master_summary(1,8)={'nb spots colloc Ch1Ch2 unique'};
     Master_summary(1,9)={'nb spots colloc Ch1Ch3 unique'};
     Master_summary(1,10)={'% colloc Ch1Ch2Ch3'}; 
     Master_summary(1,11)={'% colloc Ch1Ch2 unique'};
     Master_summary(1,12)={'% colloc Ch1Ch3 unique'};
      pp=2;
        temp=[];
    
      for itif=1:size(dirtif,1)            
        name=dirtif(itif).name(1:end-17);
        
        %get number of ROI
        load(cat(2,name,'_c1c2_summary.mat'));
        if isempty(CollocSpots)==0
        nROI=max(CollocSpots(:,20));
        clear cellPoolfoundspots3Dcleared CollocSpots CollocSpotsdiam1 CollocSpotsdiam2 
        
        for iROI=1:nROI
        
        [A] = extract_triple_colloc(name,iROI);
        A(1,7)=A(4)/A(1)*100;
        A(1,8)=A(5)/A(1)*100;
        A(1,9)=A(6)/A(1)*100;
        
        Master_summary(pp,1)={listfol(ifol).name};
        Master_summary(pp,2)={name};  
        Master_summary(pp,3)={num2str(iROI)};
        Master_summary(pp,4)={num2str(A(1))};
        Master_summary(pp,5)={num2str(A(2))};
        Master_summary(pp,6)={num2str(A(3))};
        Master_summary(pp,7)={num2str(A(4))}; 
        Master_summary(pp,8)={num2str(A(5))};
        Master_summary(pp,9)={num2str(A(6))};
        Master_summary(pp,10)={num2str(A(4)/A(1)*100)}; 
        Master_summary(pp,11)={num2str(A(5)/A(1)*100)};
        Master_summary(pp,12)={num2str(A(6)/A(1)*100)};
        pp=pp+1;
        temp=cat(1,temp,A);       
        end
        end
      end
    
    name=cat(2,cat(2,listfol(ifol).folder,'\',listfol(ifol).name),'\triplecolloc_summary.csv');
    writecell(Master_summary,name);
    
    MMaster_summary(ppp,1)={listfol(ifol).name};
    MMaster_summary(ppp,3)={num2str(size(temp,1))};
    MMaster_summary(ppp,4)={num2str(mean(temp(:,7)))}; 
    MMaster_summary(ppp,5)={num2str(std(temp(:,7))/sqrt(size(temp,1)))}; 
    MMaster_summary(ppp,6)={num2str(mean(temp(:,8)))}; 
    MMaster_summary(ppp,7)={num2str(std(temp(:,8))/sqrt(size(temp,1)))}; 
    MMaster_summary(ppp,8)={num2str(mean(temp(:,9)))};   
    MMaster_summary(ppp,9)={num2str(std(temp(:,9))/sqrt(size(temp,1)))}; 
    ppp=ppp+1;
    end
    
  name=cat(2,cat(2,listfol(ifol).folder,'\triplecolloc_summary.csv'));
    writecell(MMaster_summary,name);
