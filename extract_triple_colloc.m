function [A] = extract_triple_colloc(name,roi)
%name='aPKC-wt-60min_p8.tif';
%roi=1;

load(cat(2,name,'_c1c2_summary.mat'));
Colloc_c1c2=CollocSpots;
load(cat(2,name,'_c1c3_summary.mat'));
Colloc_c1c3=CollocSpots;

if (isempty(Colloc_c1c2)==0 && isempty(Colloc_c1c3)==0)

x=Colloc_c1c2(find(Colloc_c1c2(:,20)==roi),2);
y=Colloc_c1c3(find(Colloc_c1c3(:,20)==roi),2);

[val,pos]=intersect(x,y);

nc1=size(cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,20)==roi)),1);
nc1c2=size(unique(x),1);
nc1c3=size(unique(y),1);
nc1c2c3=size(val,1);
nc1c2only=size(unique(x),1)-size(val,1);
nc1c3only=size(unique(y),1)-size(val,1);

A=[nc1 nc1c2 nc1c3 nc1c2c3 nc1c2only nc1c3only];
else
    
    if isempty(Colloc_c1c3)==1
        x=Colloc_c1c2(find(Colloc_c1c2(:,20)==roi),2);
        nc1=size(cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,20)==roi)),1);
        nc1c3=0;
        nc1c2=size(unique(x),1);
        nc1c2c3=0;
        nc1c2only=size(unique(x),1);
        nc1c3only=0;
        A=[nc1 nc1c2 nc1c3 nc1c2c3 nc1c2only nc1c3only];    
    elseif isempty(Colloc_c1c2)==1
        y=Colloc_c1c3(find(Colloc_c1c3(:,20)==roi),2);
        nc1=size(cellPoolfoundspots3Dcleared{1,1}(find(cellPoolfoundspots3Dcleared{1,1}(:,20)==roi)),1);
        nc1c2=0;
        nc1c3=size(unique(y),1);
        nc1c2c3=0;
        nc1c2only=0;
        nc1c3only=size(unique(y),1);
        A=[nc1 nc1c2 nc1c3 nc1c2c3 nc1c2only nc1c3only];
    end       
end

end

