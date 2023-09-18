clc
clear

% Code for analysing colocalization in 3D.
% From Watson, J.L. & Kruger, L.K. et al. 2023
% Synthetic Par polarity induces cytoskeleton asymmetry in unpolarized mammalian cells

%------Global parameters

CS_detection=1; %if you want to remove the coverslip (CS) (movie has to be top to bottom, meaning CS at higher Z planes);
deltaCS=3; %number of planes to remove from bottom
CS_smooth=1;%apply a 3-plane moving average smoothing of the Z-detection when detecting the coverslip (better for cells with weird shapes data);

XYspacings=110; %nm
Zspacings=200; %nm
xyzratio=Zspacings/XYspacings;

doshift=0; %do the xyz shifting (this is very compute intensive to calculate)
Resolxy=292; %this is the coloc threshold in XY set at 80% of the resolution
Resolz=497;  %this is the coloc threshold in Z

% XY FWHM [nm]=292.8564 done on 60mer (see Ben-Sasson et al., Nature 2021)
% Z  FWHM [nm]=497.6494

%sigma for 3D gaussian fitting
% sigma_xy=Resolxy/(2*sqrt(2*log(2))*XYspacings);       %this is the sigma xy of the PSF of the scope for PSF fitting
% sigma_z=Resolz/(2*sqrt(2*log(2))*Zspacings);        %this is the sigma z  of the PSF of the scope for PSF fitting

sigma_xy=1;       %because the spacing in XY are 110nm and the one in z 200, and the psf is twice in z than xy, actually, the PSF is not asymmetrical in pixel space
sigma_z=1;        

%post_filtering

%Channel1 arrays

int_inf_threshold_ch1=0; %default 0
int_sup_threshold_ch1=68000; %default 68000
sigma_inf_threshold_ch1=0; %default 0
sigma_sup_threshold_ch1=4; %default 1000 %only keep spot below a sigma in XY of 4. this ensures we mostly work with diffraction limited spots this is in FWHM (so 2.3*sigma)
offset_inf_threshold_ch1=0; %default 0
offset_sup_threshold_ch1=68000; %default 68000
sda_a_thr_ch1=0.1;  %filter to keep only spots where the SD of the intensity divided by the intensity is lower than 0.1. Put very high if you don't want to filer (like at 1000000)
a_sigma_thr_ch1=30; %keep only spots for which the density is high enough put zero if you don't want to filter.

%Channel2 Par3

int_inf_threshold_ch2=200; %update from 100
int_sup_threshold_ch2=68000;
sigma_inf_threshold_ch2=0;
sigma_sup_threshold_ch2=4;
offset_inf_threshold_ch2=0;
offset_sup_threshold_ch2=68000; %could work on this
sda_a_thr_ch2=0.1; 
a_sigma_thr_ch2=330;

%Channel3 aPKC
int_inf_threshold_ch3=0;
int_sup_threshold_ch3=68000;
sigma_inf_threshold_ch3=0;
sigma_sup_threshold_ch3=4;
offset_inf_threshold_ch3=0;
offset_sup_threshold_ch3=68000;
sda_a_thr_ch3=0.1; 
a_sigma_thr_ch3=0;


%folder to process: all your data should be in subfolders with the name of the
%condition
rootfolder='Z:\Clustering_datasets\';

%Number of channels to process

%mastermode='c1c2'; %colloc C1c2 only
%mastermode='c1c3'; %colloc C1c3 only
mastermode='c1c2_c1c3'; %colloc C1c2 and also C1C3 

%------end of Global parameters

%processing of parameters

todo=[];
if strcmp(mastermode,'c1c2')==1
    todo{1,1}='c1c2' ;%mode
    todo{1,2}='2'   ; %scd
elseif strcmp(mastermode,'c1c3')==1
    todo{1,1}='c1c3' ;%mode
    todo{1,2}='3'   ; %scd
elseif strcmp(mastermode,'c1c2_c1c3')==1
    todo{1,1}='c1c2'; %mode
    todo{1,2}='2'   ; %scd
    todo{2,1}='c1c3' ;%mode
    todo{2,2}='3'    ; %scd
end

filter_ch1=[int_inf_threshold_ch1 int_sup_threshold_ch1 sigma_inf_threshold_ch1 sigma_sup_threshold_ch1 offset_inf_threshold_ch1 offset_sup_threshold_ch1 sda_a_thr_ch1 a_sigma_thr_ch1];
filter_ch2=[int_inf_threshold_ch2 int_sup_threshold_ch2 sigma_inf_threshold_ch2 sigma_sup_threshold_ch2 offset_inf_threshold_ch2 offset_sup_threshold_ch2 sda_a_thr_ch2 a_sigma_thr_ch2];
filter_ch3=[int_inf_threshold_ch3 int_sup_threshold_ch3 sigma_inf_threshold_ch3 sigma_sup_threshold_ch3 offset_inf_threshold_ch3 offset_sup_threshold_ch3 sda_a_thr_ch3 a_sigma_thr_ch3];
%system('taskkill /f /im excel.exe');

for ik=1:size(todo,1)

%read the parameters
    mode=todo{ik,1};
    scd=todo{ik,2};
     
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

    MMaster_summary=cell(1,25);
    MMaster_summary(1,1)={'folder'};
    MMaster_summary(1,2)={cat(2,'mean % colloc Ch1 with Ch',num2str(scd))};
    MMaster_summary(1,3)={cat(2,'median % colloc Ch1 with Ch',num2str(scd))};
    MMaster_summary(1,4)={cat(2,'sem % colloc Ch1 with Ch',num2str(scd))};
    MMaster_summary(1,5)={cat(2,'mean % colloc Ch',num2str(scd),' with Ch1')};
    MMaster_summary(1,6)={cat(2,'median % colloc Ch',num2str(scd),' with Ch1')};
    MMaster_summary(1,7)={cat(2,'sem % colloc Ch',num2str(scd),' with Ch1')};
    MMaster_summary(1,8)={cat(2,'mean % colloc Ch1 with Ch',num2str(scd),' Int')};
    MMaster_summary(1,9)={cat(2,'median % colloc Ch1 with Ch',num2str(scd),' Int')};
    MMaster_summary(1,10)={cat(2,'sem % colloc Ch1 with Ch',num2str(scd),' Int')};
    MMaster_summary(1,11)={cat(2,'mean % colloc Ch',num2str(scd),'with Ch1 int')};
    MMaster_summary(1,12)={cat(2,'median % colloc Ch',num2str(scd),'with Ch1 int')};
    MMaster_summary(1,13)={cat(2,'sem % colloc Ch',num2str(scd),'with Ch1 int')};
    MMaster_summary(1,14)={'number of ROIs'};
    MMaster_summary(1,15)={'mean intensity Ch1 spot [colloc only]'};
    MMaster_summary(1,16)={'median intensity Ch1 [colloc only]'};
    MMaster_summary(1,17)={'sem  intensity Ch1 [colloc only]'};
    MMaster_summary(1,18)={cat(2,'mean intensity Ch',num2str(scd),'[colloc only]')};
    MMaster_summary(1,19)={cat(2,'median intensity Ch',num2str(scd),'[colloc only]')};
    MMaster_summary(1,20)={cat(2,'sem intensity Ch',num2str(scd),'[colloc only]')};
    MMaster_summary(1,21)={cat(2,'total number of spots that colloc between Ch1 and Ch',num2str(scd))};
    MMaster_summary(1,22)={'total number of spots Ch1'};
    MMaster_summary(1,23)={cat(2,'total number of spots Ch',num2str(scd))};
    MMaster_summary(1,24)={'mean number of spot per FOV Ch1'};
    MMaster_summary(1,25)={cat(2,'mean number of spot per FOV Ch',num2str(scd))};
    ppp=2; 



for ifol=1:size(listfol,1)  
    cd(cat(2,listfol(ifol).folder,'\',listfol(ifol).name));
    dirtif=dir("C1-*.tif");
    CollocPool=[];
    nch1=[];
    nch2=[];
    zch1=[]; %this is the z of the spot that colocalize
    Master_summary=cell(1,23);
    Master_summary(1,1)={'folder'};
    Master_summary(1,2)={'image1'};
    Master_summary(1,3)={cat(2,'image',num2str(scd))};   
    Master_summary(1,4)={'ROI'};
    Master_summary(1,5)={'CollocThrehold'};
    Master_summary(1,6)={cat(2,'nb colloc spots Ch1Ch',num2str(scd))}; 
    Master_summary(1,7)={'nb spots Ch1'};
    Master_summary(1,8)={cat(2,'nb spots Ch',num2str(scd))};
    Master_summary(1,9)={cat(2,'% colloc Ch1Ch',num2str(scd))};
    Master_summary(1,10)={cat(2,'% colloc Ch',num2str(scd),'Ch1')};
    Master_summary(1,11)={cat(2,'% colloc Ch1Ch',num2str(scd),' intensity')};
    Master_summary(1,12)={cat(2,'% colloc Ch',num2str(scd),'Ch1 intensity')};
    Master_summary(1,13)={'Inclusion into Ch1'};
    Master_summary(1,14)={'nb colloc spots'};
    Master_summary(1,15)={cat(2,'nb spots Ch',num2str(scd))};
    Master_summary(1,16)={cat(2,'% colloc Ch',num2str(scd),'into Ch1')};
    Master_summary(1,17)={cat(2,'% colloc Ch',num2str(scd),'into Ch1 intensity')};
    Master_summary(1,18)={cat(2,'Inclusion into Ch',num2str(scd))};
    Master_summary(1,19)={'nb colloc spots'};
    Master_summary(1,20)={'nb spots Ch1'};
    Master_summary(1,21)={cat(2,'% colloc Ch1 into Ch',num2str(scd))};
    Master_summary(1,22)={cat(2,'% colloc Ch1 into Ch',num2str(scd),' intensity')};
    pp=2; 
    
    parfor itif=1:size(dirtif,1)  
        file=dirtif(itif).name(4:end);
        path3=cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\');
        if isfile(cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\',file,'_',mode,'_summary.csv'))~=1
            display(cat(2,'processing file ',listfol(ifol).folder,'\',listfol(ifol).name,'\',dirtif(itif).name));
            fit_func(file,path3,deltaCS,CS_detection,XYspacings,Zspacings,xyzratio,doshift,Resolxy,Resolz,sigma_xy,sigma_z,rootfolder,mode,CS_smooth,filter_ch1,filter_ch2,filter_ch3);             
        else
            display(cat(2,'file ',listfol(ifol).folder,'\',listfol(ifol).name,'\',dirtif(itif).name,' has previously been processed, skip'));
        end
   end
   % system('taskkill /f /im excel.exe');
    
    for itif=1:size(dirtif,1)            
        file=dirtif(itif).name(4:end);
        if exist(cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\',file,'_',mode,'_summary.mat'))>0
        load(cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\',file,'_',mode,'_summary.mat'));
            for q=1:size(summary,1)-1
                 if cell2mat(summary(q+1,6))>0
             Master_summary(pp,2:22)=summary(q+1,1:21);
             Master_summary{pp,1}=char(listfol(ifol).name);
             pp=pp+1;
                 end
            end

            if ~isempty(CollocSpots)
            CollocPool=cat(1,CollocPool,CollocSpots(find(CollocSpots(:,20)>0),:));
            end
            
            for iroi=1:max(cellPoolfoundspots3Dcleared{1,1}(:,20))
                 nch1=cat(1,nch1,size(cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,20)==iroi),1),1));
                 nch2=cat(1,nch2,size(cellPoolfoundspots3Dcleared{1,2}(find(cellPoolfoundspots3Dcleared{1,2}(:,20)==iroi),1),1));
                 if ~isempty(CollocSpots)
                temp=[];
                temp=cellPoolfoundspots3Dcleared{1,1}(:,:);
                temp(:,23)=0;
                temp(CollocSpots(:,2),23)=1;
                temp=temp(find(temp(:,20)==iroi),:);
                zch1=cat(1,zch1,temp(:,4));  %zposition and ROI for colloc spots
                 end
            end
             clear CollocSpots summary
        end
    end
    
    h=figure
        if isempty(CollocPool)~=1 & size(CollocPool,1)>1
        subplot(1,5,1)
        histfit(unique(CollocPool(find(CollocPool(:,20)>0),6).*XYspacings));
        title("FWHM channel 1 in nm for spots that colocalize");
        
        subplot(1,5,2)
        histfit(unique(CollocPool(find(CollocPool(:,20)>0),15).*XYspacings));
        title(cat(2,'FWHM channel ',num2str(scd),' in nm for spots that colocalize'));
        
        subplot(1,5,3)
        histfit(CollocPool(find(CollocPool(:,20)>0),7));
        title("Intensity channel 1 for spots that colocalize");
        
        subplot(1,5,4)
        histfit(CollocPool(find(CollocPool(:,20)>0),16));
        title(cat(2,'Intensity channel ',num2str(scd),'for spots that colocalize'));
        
        subplot(1,5,5)
        scatter(CollocPool(find(CollocPool(:,20)>0),7),CollocPool(find(CollocPool(:,20)>0),16));
        xlabel("Intensity channel 1");
        ylabel(cat(2,'Intensity channel ',num2str(scd)));
        title("Intensity ratio for spots that colocalize");
        else
            title(cat(2,"No spots were ever detected as colocalized between channel 1 and channel ",mode));
        end
        savefig(h,cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\',mode,'_master_colloc.fig'));
        set(gcf,'renderer','Painters');
        name=cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\',mode,'_master_colloc.png');
        print(name,'-dpng');
    close(h)
    
    save(cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\',mode,'_master_colloc.mat'), 'CollocPool', 'zch1');
        
    name=cat(2,listfol(ifol).folder,'\',listfol(ifol).name,'\',mode,'_master_summary.xls');
    %xlswrite(name,Master_summary);
	writecell(Master_summary,name);
    
    MMaster_summary{ppp,1}=char(listfol(ifol).name);
    MMaster_summary{ppp,2}=mean(cell2mat(Master_summary(2:end,9)));
    MMaster_summary{ppp,3}=median(cell2mat(Master_summary(2:end,9)));
    MMaster_summary{ppp,4}=std(cell2mat(Master_summary(2:end,9)))./sqrt(length(cell2mat(Master_summary(2:end,12))));
    MMaster_summary{ppp,5}=mean(cell2mat(Master_summary(2:end,10)));
    MMaster_summary{ppp,6}=median(cell2mat(Master_summary(2:end,10)));
    MMaster_summary{ppp,7}=std(cell2mat(Master_summary(2:end,10)))./sqrt(length(cell2mat(Master_summary(2:end,12))));
    MMaster_summary{ppp,8}=mean(cell2mat(Master_summary(2:end,11)));
    MMaster_summary{ppp,9}=median(cell2mat(Master_summary(2:end,11)));
    MMaster_summary{ppp,10}=std(cell2mat(Master_summary(2:end,11)))./sqrt(length(cell2mat(Master_summary(2:end,12))));
    MMaster_summary{ppp,11}=mean(cell2mat(Master_summary(2:end,12)));
    MMaster_summary{ppp,12}=median(cell2mat(Master_summary(2:end,12)));
    MMaster_summary{ppp,13}=std(cell2mat(Master_summary(2:end,12)))./sqrt(length(cell2mat(Master_summary(2:end,12))));
    MMaster_summary{ppp,14}=length(cell2mat(Master_summary(2:end,12)));
    if ~isempty(CollocPool)
    MMaster_summary{ppp,15}=mean(unique(CollocPool(:,7)));
    MMaster_summary{ppp,16}=median(unique(CollocPool(:,7)));
    MMaster_summary{ppp,17}=std(unique(CollocPool(:,7)))/sqrt( size(unique(CollocPool(:,7)),1));
    MMaster_summary{ppp,18}=mean(unique(CollocPool(:,16)));
    MMaster_summary{ppp,19}=median(unique(CollocPool(:,16)));
    MMaster_summary{ppp,20}=std(unique(CollocPool(:,16)))/sqrt( size(unique(CollocPool(:,16)),1));
    MMaster_summary{ppp,21}=size(unique(CollocPool(:,16)),1);
    else
    MMaster_summary{ppp,15}=NaN;
    MMaster_summary{ppp,16}=NaN;
    MMaster_summary{ppp,17}=NaN;
    MMaster_summary{ppp,18}=NaN;
    MMaster_summary{ppp,19}=NaN;
    MMaster_summary{ppp,20}=NaN;
    MMaster_summary{ppp,21}=0;   
    end
    
    MMaster_summary{ppp,22}=sum(nch1);
    MMaster_summary{ppp,23}=sum(nch2);
    MMaster_summary{ppp,24}=mean(nch1);
    MMaster_summary{ppp,25}=mean(nch2);  
    ppp=ppp+1; 
    
    
end
name=cat(2,listfol(ifol).folder,'\',mode,'_master_summary.xls');
%xlswrite(name,MMaster_summary);
writecell(MMaster_summary,name);
end
system('taskkill /f /im excel.exe');
