function [] = fit_func(file,path3,deltaCS,CS_detection,XYspacings,Zspacings,xyzratio,doshift,Resolxy,Resolz,sigma_xy,sigma_z,rootfolder,mode,CS_smooth,filter_ch1,filter_ch2,filter_ch3)  

%--------inputs---------        
        
        path1=cat(2,path3,'C1-',file);
        name1=cat(2,'C1-',file);
        if mode=='c1c2'
        path2=cat(2,path3,'C2-',file);
        name2=cat(2,'C2-',file);
        label='2';
        elseif mode=='c1c3'
        %cheap way of doing things we just swap channel 2 by channel3 everywhere
        path2=cat(2,path3,'C3-',file);
        name2=cat(2,'C3-',file);
        filter_ch2=filter_ch3;
        label='3';
        end

  
        %ROI loading either it's the whole image, or specific ROI) done

        ROIname=cat(2,path3,file,'.roi');
        ROIname2=cat(2,path1,'.roi');
        ROIsetname=cat(2,path3,file,'_Roiset.zip');
        ROIsetname2=cat(2,path1,'_Roiset.zip');
        
        if (isfile(ROIname)==1 || isfile(ROIname2)==1)%means only one cell
            %import ROI
            if isfile(ROIname)==1
        [sROI] = ReadImageJROI(ROIname);
            else
        [sROI] = ReadImageJROI(ROIname2);
            end
        
        skip=0;
        elseif (isfile(ROIsetname)==1 || isfile(ROIsetname2)==1) %means only multiple cells
            if isfile(ROIsetname)==1           
            [sROI] = ReadImageJROI(ROIsetname);  
            else
                [sROI] = ReadImageJROI(ROIsetname2);  
            end
            
        skip=0;
        else
            display('no ROI data, skipping dataset');
            skip=1;
        end

        if skip==0
%---------------



f=waitbar(0,'please wait');

if isfile(cat(2,path3,file,'_',mode,'_cellPoolfoundspots3D_unfiltered.mat'))~=1 %this means the channel needs to be processed
cellPoolfoundspots3Dcleared=cell(1,2);
    %process channel1
waitbar(0.2,f,'3D gaussian fitting channel 1 (1min for a 1000x700x76 stack)');

s = double(readtiff(path1));
s(s==0) = NaN;
[pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(s,[sigma_xy,sigma_z],'mode','xyzAcsr');


% [pstruct, mask] = pointSourceDetection3D(frame, sigma(mCh,:), 'Alpha', opts.Alpha,...
%         'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant,...
%         'RefineMaskLoG', false, 'WindowSize', opts.WindowSize); %#ok<PFBNS>
% 
% 

pstruct.zfloor=floor(pstruct.z); %add floorz for display
[sx sy sz]=size(s);

    %put stuff in cellPoolfoundspots3Dcleared matrix, 
    % 1)Channel
    % 2) X in pix
    % 3) Y in pix
    % 4) Z in XY pix
    % 5) width this is in FMWH in XY pixels
    % 6) integrated intensity 
    % 7) offset
    % 8-13) precision on these parameters
    % 8) precision X in pix
    % 9) precision Y in pix
    % 10) precision Z in XY pix
    % 11) precision width this is in FMWH in XY pixels
    % 12) precision integrated intensity 
    % 13) precision offset
    % 20)=ROI
    % 21)=cluster of spots connected in 3D
    % 22)= integrated intensity 

cellPoolfoundspots3Dcleared{1,1}=zeros(size(pstruct.zfloor,2),22);
cellPoolfoundspots3Dcleared{1,1}(:,1)=1;
cellPoolfoundspots3Dcleared{1,1}(:,2)=(pstruct.x)';
cellPoolfoundspots3Dcleared{1,1}(:,3)=(pstruct.y)';
cellPoolfoundspots3Dcleared{1,1}(:,4)=(pstruct.z.*xyzratio)';%needs conversion to pixel size
cellPoolfoundspots3Dcleared{1,1}(:,5)=(pstruct.s(1,:).*2*sqrt(2*log(2)))'; %this is FWHM
cellPoolfoundspots3Dcleared{1,1}(:,6)=pstruct.A';
cellPoolfoundspots3Dcleared{1,1}(:,7)=pstruct.c';
cellPoolfoundspots3Dcleared{1,1}(:,8)=(pstruct.x_pstd)';
cellPoolfoundspots3Dcleared{1,1}(:,9)=(pstruct.y_pstd)';
cellPoolfoundspots3Dcleared{1,1}(:,10)=(pstruct.z_pstd.*xyzratio)';%needs conversion to pixel size
cellPoolfoundspots3Dcleared{1,1}(:,11)=(pstruct.s_pstd(1,:).*2*sqrt(2*log(2)))'; %this is FWHM
cellPoolfoundspots3Dcleared{1,1}(:,12)=pstruct.A_pstd';
cellPoolfoundspots3Dcleared{1,1}(:,13)=pstruct.c_pstd;
cellPoolfoundspots3Dcleared{1,1}(:,22)=cellPoolfoundspots3Dcleared{1,1}(:,12)./cellPoolfoundspots3Dcleared{1,1}(:,6);
%add ROIs 
if isempty(cellPoolfoundspots3Dcleared{1,1})==0
cellPoolfoundspots3Dcleared{1,1}(:,20)=0; %default value outside the cell
    for i=1:size(sROI,2)
        bc=isinROI(cellPoolfoundspots3Dcleared{1,1}(:,2),cellPoolfoundspots3Dcleared{1,1}(:,3),sROI,i);
        cc=find(bc==1);
        cellPoolfoundspots3Dcleared{1,1}(cc,20)=i;                             
    end    
end


%process channel2
waitbar(0.4,f,'3D gaussian fitting channel 2 (1min for a 1000x700x76 stack)');

s = double(readtiff(path2));
s(s==0) = NaN;

[pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(s,[sigma_xy,sigma_z],'mode','xyzAcsr');
pstruct.zfloor=floor(pstruct.z);%add floorz for display
[sx sy sz]=size(s);

cellPoolfoundspots3Dcleared{1,2}=zeros(size(pstruct.zfloor,2),22);
cellPoolfoundspots3Dcleared{1,2}(:,1)=2;
cellPoolfoundspots3Dcleared{1,2}(:,2)=(pstruct.x)';
cellPoolfoundspots3Dcleared{1,2}(:,3)=(pstruct.y)';
cellPoolfoundspots3Dcleared{1,2}(:,4)=(pstruct.z.*xyzratio)';%needs conversion to pixel size
cellPoolfoundspots3Dcleared{1,2}(:,5)=(pstruct.s(1,:).*2*sqrt(2*log(2)))'; %this is FWHM
cellPoolfoundspots3Dcleared{1,2}(:,6)=pstruct.A';
cellPoolfoundspots3Dcleared{1,2}(:,7)=pstruct.c';
cellPoolfoundspots3Dcleared{1,2}(:,8)=(pstruct.x_pstd)';
cellPoolfoundspots3Dcleared{1,2}(:,9)=(pstruct.y_pstd)';
cellPoolfoundspots3Dcleared{1,2}(:,10)=(pstruct.z_pstd.*xyzratio)';%needs conversion to pixel size
cellPoolfoundspots3Dcleared{1,2}(:,11)=(pstruct.s_pstd(1,:).*2*sqrt(2*log(2)))'; %this is FWHM
cellPoolfoundspots3Dcleared{1,2}(:,12)=pstruct.A_pstd';
cellPoolfoundspots3Dcleared{1,2}(:,13)=pstruct.c_pstd';
cellPoolfoundspots3Dcleared{1,2}(:,22)=cellPoolfoundspots3Dcleared{1,2}(:,12)./cellPoolfoundspots3Dcleared{1,2}(:,6);

%add ROIs
if isempty(cellPoolfoundspots3Dcleared{1,2})==0
cellPoolfoundspots3Dcleared{1,2}(:,20)=0; %default value outside the cell
    for i=1:size(sROI,2)
        bc=isinROI(cellPoolfoundspots3Dcleared{1,2}(:,2),cellPoolfoundspots3Dcleared{1,2}(:,3),sROI,i);
        cc=find(bc==1);
        cellPoolfoundspots3Dcleared{1,2}(cc,20)=i;                             
    end    
end

save(cat(2,path3,file,'_',mode,'_cellPoolfoundspots3D_unfiltered.mat'), 'cellPoolfoundspots3Dcleared');


else
    display(cat(2,'bypassing 3D detection for ',path3,file,'_',mode));
    load(cat(2,path3,file,'_',mode,'_cellPoolfoundspots3D_unfiltered.mat'),'cellPoolfoundspots3Dcleared');
    
    s = double(readtiff(path1));
    [sx sy sz]=size(s); %this is just to get the size, we could use less RAM
end

%filtering 
waitbar(0.45,f,'filtering of both channels');
%ch1
int_inf_threshold_ch1=filter_ch1(1,1); 
int_sup_threshold_ch1=filter_ch1(1,2); 
sigma_inf_threshold_ch1=filter_ch1(1,3);  
sigma_sup_threshold_ch1=filter_ch1(1,4); 
offset_inf_threshold_ch1=filter_ch1(1,5); 
offset_sup_threshold_ch1=filter_ch1(1,6); 
sda_a_thr_ch1=filter_ch1(1,7);  
a_sigma_thr_ch1=filter_ch1(1,8)/(2*sqrt(2*log(2))); 

b=find(cellPoolfoundspots3Dcleared{1,1}(:,22)<sda_a_thr_ch1);
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,1}(:,5)<sigma_sup_threshold_ch1.*2*sqrt(2*log(2)));
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,1}(:,5)>sigma_inf_threshold_ch1.*2*sqrt(2*log(2)));
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,1}(:,6)<int_sup_threshold_ch1);
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,1}(:,6)>int_inf_threshold_ch1);
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,1}(:,7)>offset_inf_threshold_ch1);
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,1}(:,7)<offset_sup_threshold_ch1);
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,1}(:,6)./cellPoolfoundspots3Dcleared{1,1}(:,5)>a_sigma_thr_ch1);
cellPoolfoundspots3Dcleared{1,1}=cellPoolfoundspots3Dcleared{1,1}(b,:);

%ch2

int_inf_threshold_ch2=filter_ch2(1,1); 
int_sup_threshold_ch2=filter_ch2(1,2); 
sigma_inf_threshold_ch2=filter_ch2(1,3);  
sigma_sup_threshold_ch2=filter_ch2(1,4); 
offset_inf_threshold_ch2=filter_ch2(1,5); 
offset_sup_threshold_ch2=filter_ch2(1,6); 
sda_a_thr_ch2=filter_ch2(1,7);
a_sigma_thr_ch2=filter_ch2(1,8)/(2*sqrt(2*log(2)));

b=find(cellPoolfoundspots3Dcleared{1,2}(:,22)<sda_a_thr_ch2);
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,2}(:,5)<sigma_sup_threshold_ch2.*2*sqrt(2*log(2)));
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,2}(:,5)>sigma_inf_threshold_ch2.*2*sqrt(2*log(2)));
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,2}(:,6)<int_sup_threshold_ch2);
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,2}(:,6)>int_inf_threshold_ch2);
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,2}(:,7)>offset_inf_threshold_ch2);
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,2}(:,7)<offset_sup_threshold_ch2);
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

b=find(cellPoolfoundspots3Dcleared{1,2}(:,6)./cellPoolfoundspots3Dcleared{1,2}(:,5)>a_sigma_thr_ch2);
cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(b,:);

waitbar(0.5,f,'measuring colocalisation in 3D');
%----------Colocalization analysis---------
%For the colloc, We use 3 methods

% 1) distance between centroids is smaller to a given threshold : CollocSpots
% 2) inclusion of Channel2 centroid into Channel1 sphere [diameter is Gaussian width at max intensity] : CollocSpotsdiam1
% 3) inclusion of Channel1 centroid into Channel2 sphere [diameter is Gaussian width at max intensity]: CollocSpotsdiam2

%Syntax of these matrixes [XYZ coordinates are expressed in XY spacings in
%pixels]

% 1)Channel 1
% 2)number
% 3)X [px]
% 4)Y [px]
% 5)Z [px][XY spacings]
% 6)Width [ XY spacings]
% 7)Intensity (Max)
% 8)offset
% 9)Intensity (Integrated)
% 10)Channel 2
% 11)number
% 12)X [px]
% 13)Y [px]
% 14)Z [px]
% 15)Width [px]
% 16)Intensity (Max)
% 17)offset
% 18) Intensity (Integrated)
% 19) D
% 20) ROI

%-----------Distance threshold method--------
%In 2D, we say that 2 spot colocalize if the distance between them is below
%the lateral resolution of the scope: Resolxy. It's more complex in 3D,
%because the axial resolution Resolz is higher than the lateral one, we thus
%should not consider a sphere but an ovoid of big vertical axis Resolz and small
%horizontal axis Resolxy. The idea will thus to calculate the intersection of the
% theoretical ovoid PSF around the first point, then to calculate the distance 
% between the first point and this intersection, then to compare it to the 
%distance between the two centroids. If it's bellow, it's colocalizing
%CF "JACoP v2.0: improving the user experience with co-localization
%studies"

nspot1=size(cellPoolfoundspots3Dcleared{1,1},1);
nspot2=size(cellPoolfoundspots3Dcleared{1,2},1);
p=1;
CollocSpots=[];
Phi=zeros(nspot1,nspot2);
Teta=zeros(nspot1,nspot2);
Rref=zeros(nspot1,nspot2);
D=zeros(nspot1,nspot2);
for i=1:nspot1
    for j=1:nspot2
        x1=cellPoolfoundspots3Dcleared{1,1}(i,2)*XYspacings; %takes cycles to define new variables each time, but otherwise it's messy
        x2=cellPoolfoundspots3Dcleared{1,2}(j,2)*XYspacings;
        y1=cellPoolfoundspots3Dcleared{1,1}(i,3)*XYspacings;
        y2=cellPoolfoundspots3Dcleared{1,2}(j,3)*XYspacings;
        z1=cellPoolfoundspots3Dcleared{1,1}(i,4)*XYspacings; % Z is expressed in XY spacings, 
        z2=cellPoolfoundspots3Dcleared{1,2}(j,4)*XYspacings; % Z is expressed in XY spacings, 
        %this is for development purposes when you enter twice the same
        %image. there is a div by zero problem otherwise.
        if x1==x2 | y1==y2
            x1=x1+0.005;
            y1=y1+0.005;
        end
        Phi(i,j)=acos((x2-x1)/sqrt((x2-x1)^2+(y2-y1)^2));   
        Teta(i,j)=acos((z2-z1)/sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2));                     
        Rref(i,j)=sqrt((Resolxy*sin(Teta(i,j))*cos(Phi(i,j)))^2+(Resolxy*sin(Teta(i,j))*sin(Phi(i,j)))^2+(Resolz*cos(Teta(i,j)))^2);
        D(i,j)=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);           
        if D(i,j)<Rref(i,j)
        CollocSpots(p,1)=1;   
        CollocSpots(p,2)=i;
        CollocSpots(p,3:8)=cellPoolfoundspots3Dcleared{1,1}(i,2:7);
        CollocSpots(p,9)=cellPoolfoundspots3Dcleared{1,1}(i,22);
        CollocSpots(p,10)=2; 
        CollocSpots(p,11)=j;
        CollocSpots(p,12:17)=cellPoolfoundspots3Dcleared{1,2}(j,2:7);
        CollocSpots(p,18)=cellPoolfoundspots3Dcleared{1,2}(j,22);
        CollocSpots(p,19)=D(i,j); 
        CollocSpots(p,21)=Rref(i,j); 
        p=p+1;             
        end
    end
end

%look if spot in channel2 is included within gaussian width of spot in channel1
p=1;
CollocSpotsdiam1=[];
for i=1:nspot1
    for j=1:nspot2       
        %Rrefdiam1(i,j)=sqrt((cellPoolfoundspots3Dcleared{1,1}(i,5)*XYspacings*sin(Teta(i,j))*cos(Phi(i,j)))^2+(cellPoolfoundspots3Dcleared{1,1}(i,5)*XYspacings*sin(Teta(i,j))*sin(Phi(i,j)))^2+(cellPoolfoundspots3Dcleared{1,1}(i,5)*XYspacings*Resolz/Resolxy*cos(Teta(i,j)))^2);             
        %if D(i,j)<Rrefdiam1(i,j) %could put just the width here cellPoolfoundspots3Dcleared{1,1}(i,5)*XYspacings
        if D(i,j)<cellPoolfoundspots3Dcleared{1,1}(i,5)*XYspacings
        CollocSpotsdiam1(p,1)=1;   
        CollocSpotsdiam1(p,2)=i;
        CollocSpotsdiam1(p,3:8)=cellPoolfoundspots3Dcleared{1,1}(i,2:7);
        CollocSpotsdiam1(p,9)=cellPoolfoundspots3Dcleared{1,1}(i,22);
        CollocSpotsdiam1(p,10)=2; 
        CollocSpotsdiam1(p,11)=j;
        CollocSpotsdiam1(p,12:17)=cellPoolfoundspots3Dcleared{1,2}(j,2:7);
        CollocSpotsdiam1(p,18)=cellPoolfoundspots3Dcleared{1,2}(j,22);
        CollocSpotsdiam1(p,19)=D(i,j); 
        p=p+1;   
       end
    end
end

%look if spot in channel1 is included within gaussian width of spot in channel2
p=1;
CollocSpotsdiam2=[];
for i=1:nspot1
    for j=1:nspot2              
        %Rrefdiam2(i,j)=sqrt((cellPoolfoundspots3Dcleared{1,2}(j,5)*XYspacings*sin(Teta(i,j))*cos(Phi(i,j)))^2+(cellPoolfoundspots3Dcleared{1,2}(j,5)*XYspacings*sin(Teta(i,j))*sin(Phi(i,j)))^2+(cellPoolfoundspots3Dcleared{1,2}(j,5)*XYspacings*Resolz/Resolxy*cos(Teta(i,j)))^2);                  
        %if D(i,j)<Rrefdiam2(i,j)%could put just the width here cellPoolfoundspots3Dcleared{1,2}(j,5)*XYspacings
        if D(i,j)<cellPoolfoundspots3Dcleared{1,2}(j,5)*XYspacings
        CollocSpotsdiam2(p,1)=1;   
        CollocSpotsdiam2(p,2)=i;
        CollocSpotsdiam2(p,3:8)=cellPoolfoundspots3Dcleared{1,1}(i,2:7);
        CollocSpotsdiam2(p,9)=cellPoolfoundspots3Dcleared{1,1}(i,22);
        CollocSpotsdiam2(p,10)=2; 
        CollocSpotsdiam2(p,11)=j;
        CollocSpotsdiam2(p,12:17)=cellPoolfoundspots3Dcleared{1,2}(j,2:7);
        CollocSpotsdiam2(p,18)=cellPoolfoundspots3Dcleared{1,2}(j,22);
        CollocSpotsdiam2(p,19)=D(i,j); 
        p=p+1;   
       end
    end
end
%now we will add the ROI the spot is in in the 20th collumn of each matrix
%default value outside the cell: ROI=0
if isempty(CollocSpotsdiam2)==0
CollocSpotsdiam2(:,20)=0;
    for i=1:size(sROI,2)
        b=isinROI(CollocSpotsdiam2(:,3),CollocSpotsdiam2(:,4),sROI,i);
        c=find(b==1);
        CollocSpotsdiam2(c,20)=i; 
    end
end

if isempty(CollocSpotsdiam1)==0
CollocSpotsdiam1(:,20)=0; 
    for i=1:size(sROI,2)      
        b=isinROI(CollocSpotsdiam1(:,3),CollocSpotsdiam1(:,4),sROI,i);
        c=find(b==1);
        CollocSpotsdiam1(c,20)=i;
    end
end

if isempty(CollocSpots)==0
CollocSpots(:,20)=0;
 for i=1:size(sROI,2)
        b=isinROI(CollocSpots(:,3),CollocSpots(:,4),sROI,i);
        c=find(b==1);
        CollocSpots(c,20)=i;        
    end
end
%end of 3D colloc


waitbar(0.75,f,'Coverslip detection and removal');

for iroi=1:size(sROI,2)
Zplane=zeros(sz,3);
Zplane(1:sz,1)=1:sz';

temp=floor(cellPoolfoundspots3Dcleared{1,1}(:,4)./xyzratio);
temp=temp(find(cellPoolfoundspots3Dcleared{1,1}(:,20)==iroi),:);
for i=1:sz
Zplane(i,2)=size(find(temp(:)==i),1);
end
Zplane(:,3)=Zplane(:,2)./sum(Zplane(:,2));
Zplane(:,4)=movmean(Zplane(:,3),3);
h=figure;
plot(Zplane(1:sz,1),Zplane(1:sz,3),'-r');
hold on
plot(Zplane(1:sz,1),Zplane(1:sz,4),'--r');
hold on
if CS_detection==1
    [~,bottom_plane]=min(diff(Zplane(1:sz,3)));
    bottom_plane=bottom_plane-1-deltaCS;
    [~,bottom_plane_s]=min(diff(Zplane(1:sz,4))); %CS detection done on smoothed data
    bottom_plane_s=bottom_plane_s-1-deltaCS;
    plot(Zplane(1:bottom_plane,1),Zplane(1:bottom_plane,3),'-g');
    plot(Zplane(1:bottom_plane_s,1),Zplane(1:bottom_plane_s,3),'--g');
    
    if CS_smooth==1
        bottom_plane=bottom_plane_s; %this way the sorting is done on the smoothed value
    end
    
    if ~isempty(cellPoolfoundspots3Dcleared{1,1}) 
    cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,4)>(bottom_plane*xyzratio) & cellPoolfoundspots3Dcleared{1,1}(:,20)==iroi ),20)=0;
    end
     if ~isempty(cellPoolfoundspots3Dcleared{1,2}) 
    cellPoolfoundspots3Dcleared{1,2}(find(cellPoolfoundspots3Dcleared{1,2}(:,4)>(bottom_plane*xyzratio) & cellPoolfoundspots3Dcleared{1,2}(:,20)==iroi ),20)=0;
     end
    %cellPoolfoundspots3Dcleared{1,2}=cellPoolfoundspots3Dcleared{1,2}(find(cellPoolfoundspots3Dcleared{1,2}(:,4)<(bottom_plane*xyzratio) ),:);
    %CollocSpots=CollocSpots(find(CollocSpots(:,5)<(bottom_plane*xyzratio)),:);
    if ~isempty(CollocSpots)
    CollocSpots(find(CollocSpots(:,5)>(bottom_plane*xyzratio) & CollocSpots(:,20)==iroi),20)=0;
    end
    if ~isempty(CollocSpotsdiam1)
    CollocSpotsdiam1(find(CollocSpotsdiam1(:,5)>(bottom_plane*xyzratio) & CollocSpotsdiam1(:,20)==iroi),20)=0;
    end
     if ~isempty(CollocSpotsdiam2)
    CollocSpotsdiam2(find(CollocSpotsdiam2(:,5)>(bottom_plane*xyzratio) & CollocSpotsdiam2(:,20)==iroi),20)=0;
     end
    
    %CollocSpotsdiam1=CollocSpotsdiam1(find(CollocSpotsdiam1(:,5)<(bottom_plane*xyzratio)),:);
    %CollocSpotsdiam2=CollocSpotsdiam2(find(CollocSpotsdiam2(:,5)<(bottom_plane*xyzratio)),:);
        

end
set(gcf,'renderer','Painters');
print(cat(2,path1,'_roi',num2str(iroi),'plot_Zdetect.png'),'-dpng');
close(h);
   
end


%make a plot
h=figure;
if ~isempty(cellPoolfoundspots3Dcleared{1,1}) && ~isempty(cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,20)>0)))
scatter(cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,20)>0),2),cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,20)>0),3),'xg');
end
hold on
%scatter(cellPoolfoundspots3Dcleared{1,2}(find(cellPoolfoundspots3Dcleared{1,2}(:,20)>0),2),cellPoolfoundspots3Dcleared{1,2}(find(cellPoolfoundspots3Dcleared{1,2}(:,20)>0),3),'xr');
if ~isempty(CollocSpots) && ~isempty(CollocSpots(find(CollocSpots(:,20)>0)))
scatter(CollocSpots(find(CollocSpots(:,20)>0),3),CollocSpots(find(CollocSpots(:,20)>0),4),'xb');
end
if size(sROI,2)>1
    for iroi=1:size(sROI,2)
        plot(sROI{1,iroi}.mnCoordinates(:,1),sROI{1,iroi}.mnCoordinates(:,2));
    end
else
    try
    plot(sROI(1,iroi).mnCoordinates(:,1),sROI(1,iroi).mnCoordinates(:,2));
    catch ME
      if(strcmp(ME.identifier,'MATLAB:structRefFromNonStruct')==1)
     plot(sROI{1,iroi}.mnCoordinates(:,1),sROI{1,iroi}.mnCoordinates(:,2)); 
      end
    end
end

set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
set(gcf,'renderer','Painters');
title(cat(2,'detected spots ch1 in green, colloc spots with chan ',label,' in blue'));
print(cat(2,path3,file,'_',mode,'_plotcollocproj.png'),'-dpng');
close(h);

waitbar(0.8,f,'exporting');

%make IJ output (hijacks thunderstorm)
if ~isempty(cellPoolfoundspots3Dcleared{1,2})
Matrix=zeros(size(cellPoolfoundspots3Dcleared{1,2},1),10);
Matrix(:,1)=1:size(cellPoolfoundspots3Dcleared{1,2},1);
Matrix(:,3)=(cellPoolfoundspots3Dcleared{1,2}(:,2)-1); %we have to remove 1 for display between IJ and Matlab...
Matrix(:,4)=(cellPoolfoundspots3Dcleared{1,2}(:,3)-1); %we have to remove 1 for display between IJ and Matlab...
Matrix(:,2)=floor(cellPoolfoundspots3Dcleared{1,2}(:,4)./xyzratio);
Matrix(:,5)=cellPoolfoundspots3Dcleared{1,2}(:,5);
Matrix(:,6)=cellPoolfoundspots3Dcleared{1,2}(:,6);
Matrix(:,7)=cellPoolfoundspots3Dcleared{1,2}(:,7);
Matrix(:,10)=cellPoolfoundspots3Dcleared{1,2}(:,20);
Matrix(:,11)=cellPoolfoundspots3Dcleared{1,2}(:,22);
header = {'id','frame','x [pixel]','y [pixel]','sigma [pixel]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]','ROI','SDA/A'};
Matrix = [header; num2cell(Matrix)];
writecell(Matrix,cat(2,path2,'_coord_3Ddetec.csv'));
end
%make IJ output (hijacks thunderstorm)
if ~isempty(cellPoolfoundspots3Dcleared{1,1})
Matrix=zeros(size(cellPoolfoundspots3Dcleared{1,1},1),10);
Matrix(:,1)=1:size(cellPoolfoundspots3Dcleared{1,1},1);
Matrix(:,3)=(cellPoolfoundspots3Dcleared{1,1}(:,2)-1); %we have to remove 1 for display between IJ and Matlab...
Matrix(:,4)=(cellPoolfoundspots3Dcleared{1,1}(:,3)-1); %we have to remove 1 for display between IJ and Matlab...
Matrix(:,2)=floor(cellPoolfoundspots3Dcleared{1,1}(:,4)./xyzratio);
Matrix(:,5)=cellPoolfoundspots3Dcleared{1,1}(:,5);
Matrix(:,6)=cellPoolfoundspots3Dcleared{1,1}(:,6);
Matrix(:,7)=cellPoolfoundspots3Dcleared{1,1}(:,7);
Matrix(:,10)=cellPoolfoundspots3Dcleared{1,1}(:,20);
Matrix(:,11)=cellPoolfoundspots3Dcleared{1,1}(:,22);
header = {'id','frame','x [pixel]','y [pixel]','sigma [pixel]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]','ROI','SDA/A'};
Matrix = [header; num2cell(Matrix)];
writecell(Matrix,cat(2,path1,'_coord_3Ddetec.csv'));
end
%make IJ output (hijacks thunderstorm)
if ~isempty(CollocSpots)
Matrix=zeros(size(CollocSpots,1),10);
Matrix(:,1)=1:size(CollocSpots,1);
Matrix(:,3)=(CollocSpots(:,3)-1)'; %we have to remove 1 for display between IJ and Matlab...
Matrix(:,4)=(CollocSpots(:,4)-1)'; %we have to remove 1 for display between IJ and Matlab...
Matrix(:,2)=floor(CollocSpots(:,5)./xyzratio);
Matrix(:,5:7)=CollocSpots(:,6:8);
Matrix(:,10)=CollocSpots(:,20);
header = {'id','frame','x [pixel]','y [pixel]','sigma [pixel]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]','ROI'};
Matrix = [header; num2cell(Matrix)];
writecell(Matrix,cat(2,path3,file,'_',mode,'_coord_3Dcolloc.csv'));
end
%make IJ output (hijacks thunderstorm)
if ~isempty(CollocSpotsdiam1)
Matrix=zeros(size(CollocSpotsdiam1,1),10);
Matrix(:,1)=1:size(CollocSpotsdiam1,1);
Matrix(:,3)=(CollocSpotsdiam1(:,3)-1)'; %we have to remove 1 for display between IJ and Matlab...
Matrix(:,4)=(CollocSpotsdiam1(:,4)-1)'; %we have to remove 1 for display between IJ and Matlab...
Matrix(:,2)=floor(CollocSpotsdiam1(:,5)./xyzratio);
Matrix(:,5:7)=CollocSpotsdiam1(:,6:8);
Matrix(:,10)=CollocSpotsdiam1(:,20);
header = {'id','frame','x [pixel]','y [pixel]','sigma [pixel]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]','ROI'};
Matrix = [header; num2cell(Matrix)];
writecell(Matrix,cat(2,path3,file,'_',mode,'_coord_3Dcollocdiam1.csv'));
end
if ~isempty(CollocSpotsdiam2)
%make IJ output (hijacks thunderstorm)
Matrix=zeros(size(CollocSpotsdiam2,1),10);
Matrix(:,1)=1:size(CollocSpotsdiam2,1);
Matrix(:,3)=(CollocSpotsdiam2(:,12)-1)'; %we have to remove 1 for display between IJ and Matlab...
Matrix(:,4)=(CollocSpotsdiam2(:,13)-1)'; %we have to remove 1 for display between IJ and Matlab...
Matrix(:,2)=floor(CollocSpotsdiam2(:,14)./xyzratio);
Matrix(:,5:7)=CollocSpotsdiam2(:,6:8);
Matrix(:,10)=CollocSpotsdiam2(:,20);
header = {'id','frame','x [pixel]','y [pixel]','sigma [pixel]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]','ROI'};
Matrix = [header; num2cell(Matrix)];
writecell(Matrix,cat(2,path3,file,'_',mode,'_coord_3Dcollocdiam2.csv'));
end
%excel output (one extended and one summarized)

output=[];
summary=cell(1,21);
summary(1,1)={'image1'};
summary(1,2)={'image2'};   
summary(1,3)={'ROI'};
summary(1,4)={'CollocThrehold'};
summary(1,5)={'nb colloc spots'};
summary(1,6)={'nb spots Ch1'};
summary(1,7)={cat(2,'nb spots Ch',mode)};
summary(1,8)={cat(2,'% colloc Ch1Ch',mode)};
summary(1,9)={cat(2,'% colloc Ch',mode,'Ch1')};
summary(1,10)={cat(2,'% colloc Ch1Ch',mode,' intensity')};
summary(1,11)={cat(2,'% colloc Ch',mode,'Ch1 intensity')};
summary(1,12)={'Inclusion into Ch1'};
summary(1,13)={'nb colloc spots'};
summary(1,14)={cat(2,'nb spots Ch',mode)};
summary(1,15)={cat(2,'% spots Ch',mode,'into Ch1')};
summary(1,16)={cat(2,'% spots Ch',mode,'into Ch1 intensity')};
summary(1,17)={cat(2,'Inclusion into Ch',mode)};
summary(1,18)={'nb colloc spots'};
summary(1,19)={'nb spots Ch1'};
summary(1,20)={cat(2,'% spots Ch1 into Ch',mode)};
summary(1,21)={cat(2,'% spots Ch1 into Ch',mode,' intensity')};

summary2=[];%Summary2 is just summary without the header for easy sending to imageJ

for i=1:size(sROI,2)
    [outputtemp,summarytemp]=generateexceloutputROI(cat(1,cellPoolfoundspots3Dcleared{1,:}),CollocSpots,Resolxy,Resolz,CollocSpotsdiam1,CollocSpotsdiam2,i,name1,name2);
    output=cat(1,output,outputtemp);
    summary2=cat(1,summary2,summarytemp);
    summarytemp=num2cell(summarytemp);
    summarytemp(1,1)={name1};
    summarytemp(1,2)={name2}; 
    summarytemp(1,4)={cat(2,num2str(Resolxy),' um lateral and ',num2str(Resolz),' um axial')};    
    summarytemp{1,12}=[];
    summarytemp{1,17}=[]; 
    summary=cat(1,summary,summarytemp);  
end
name=cat(2,path3,file,'_',mode,'_extended_results.csv');
%xlswrite(name,output);
writecell(output,name);

name=cat(2,path3,file,'_',mode,'_summary.csv');
% xlswrite(name,summary);
writecell(summary,name);

if doshift==1
    waitbar(0.9,f,'doing the XYZ shift matrix');
    
    dX=-20:1:20;
    dY=-20:1:20;
    dZ=-20:1:20;
    Xshift=cell(size(dX,2),1);
    Yshift=cell(size(dY,2),1);
    Zshift=cell(size(dZ,2),1);
    parfor i=1:size(dX,2)
    Xshift{i,1} = colloc_shift(cellPoolfoundspots3Dcleared,dX(i),0,0,XYspacings,Resolxy,Resolz)
    end
    parfor i=1:size(dY,2)
    Yshift{i,1} = colloc_shift(cellPoolfoundspots3Dcleared,0,dY(i),0,XYspacings,Resolxy,Resolz)
    end
    parfor i=1:size(dZ,2)
    Zshift{i,1} = colloc_shift(cellPoolfoundspots3Dcleared,0,0,dZ(i),XYspacings,Resolxy,Resolz)
    end
    
    Matrix=zeros(size(dX,2),4);
    Matrix(:,1)=dX';
    Matrix(:,2)=cell2mat(Xshift(:,1));
    Matrix(:,3)=cell2mat(Yshift(:,1));
    Matrix(:,4)=cell2mat(Zshift(:,1));
    header = {'shift (px)','in X','in Y','in Z'};
    Matrix = [header; num2cell(Matrix)];
    writecell(Matrix,cat(2,path3,file,'_',mode,'_pixelshift.csv'));
end
close(f);

    save(cat(2,path3,file,'_',mode,'_summary.mat'), 'summary', 'cellPoolfoundspots3Dcleared', 'CollocSpots','CollocSpotsdiam2','CollocSpotsdiam1','filter_ch1','filter_ch2');


clear temp b bc imgLM imgLoG CollocSpots CollocSpotsdiam1 CollocSpotsdiam2 cellPoolfoundspots3Dcleared output outputtemp D Rref s summary summary2 summarytemp Teta
end

