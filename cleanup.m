%this code will go in all the subfolders and clean up all the intermediate files.
%folder to process: all your data should be in folders with the name of the
%condition
rootfolder='Z:\Clustering_datasets\';


%------end of Global parameters

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

 
for ifol=1:size(listfol,1)  
    cd(cat(2,listfol(ifol).folder,'\',listfol(ifol).name));
    display(cat(2,'Processing ',listfol(ifol).folder,'\',listfol(ifol).name));
    system('del *.txt');
    system('del *.csv');
    system('del *.m');
    %system('del *.mat');
    system('del *.xls');
    system('del *.png');
    system('del *.fig');
end
