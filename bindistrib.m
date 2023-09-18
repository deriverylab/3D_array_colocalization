function [Distrib,Databinned] = bindistrib(data,bins)
Distrib=zeros(length(bins)-1,1);
Databinned={};
    for i=1:length(bins)-1
    b=find(data(:,1)>=bins(i) & data(:,1)<bins(i+1));
    Distrib(i,1)=length(b);
    Databinned{i,1}=data(b);
    end
end

