function [ouput] = isinROI(x,y,ROI,n)

%input x,y coordinates of point, ROI=name of the ROI structure (from
%ReadImageJROI function), n=number of the ROI you want to consider
%output =1 of inside ROI, 0 otherwise

% The field '.strType' is guaranteed to exist, and defines the ROI type:
% {'Rectangle', 'Oval', Line', 'Polygon', 'Freehand', 'Traced', 'PolyLine',
% 'FreeLine', 'Angle', 'Point', 'NoROI'}.
if iscell(ROI)==1

if strcmp(char(ROI{1,n}.strType),'Freehand')
    Polyx=ROI{1,n}.mnCoordinates(:,1);
    Polyy=ROI{1,n}.mnCoordinates(:,2);
elseif strcmp(char(ROI{1,n}.strType),'Rectangle')
    B=ROI{1,n}.vnRectBounds;            
    Polyy=[B(1),B(1),B(3),B(3)]';
    Polyx=[B(2),B(4),B(4),B(2)]';
elseif strcmp(char(ROI{1,n}.strType),'Polygon')
    Polyx=ROI{1,n}.mnCoordinates(:,1);
    Polyy=ROI{1,n}.mnCoordinates(:,2);      
end


else
   if strcmp(char(ROI(1,n).strType),'Freehand')
    Polyx=ROI(1,n).mnCoordinates(:,1);
    Polyy=ROI(1,n).mnCoordinates(:,2);
elseif strcmp(char(ROI(1,n).strType),'Rectangle')
    B=ROI(1,n).vnRectBounds;            
    Polyy=[B(1),B(1),B(3),B(3)]';
    Polyx=[B(2),B(4),B(4),B(2)]';
elseif strcmp(char(ROI(1,n).strType),'Polygon')
    Polyx=ROI(1,n).mnCoordinates(:,1);
    Polyy=ROI(1,n).mnCoordinates(:,2);      
end 
end

ouput = inpolygon(x,y,Polyx,Polyy);

end

