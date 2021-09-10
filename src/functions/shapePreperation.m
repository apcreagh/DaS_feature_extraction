function [shape, x, y, x1, y1, t,xref, yref, t_interpol]=shapePreperation(shapeData, ref_shape, index, mtype, attempt)
%% Draw-a-Shape: Shape Preperation Function:
% Function to extract the x- and y-coordinate touch sceen points for a
% shape; the corresponding timestamps; the reference way-point coordinates;
% the ideal shape as interpolation between reference way-point coordinates,
% along with corresponding timestamps.
%--------------------------------------------------------------------------
% Input:
%       - shapeData: [N x 2] matrix of x- and y-coordinate touch screen data 
%                    (corresponding to shapes drawn).
%       - ref_shape: [N x 2] matrix of x- and y-coordinate reference way-point
%                    coordinates. 
%       - index: the index corresponding to the selected shape to pull.
%       - mtype: the interpolation method (string), denoting the number of
%                points to interpolate to:
%             - 'fixed': interpolate to a fixed number of touch screen pts
%             - 'length': interpolate to the length of drawn touch screen pts
%       - attempt: the drawing attempt to extract (1st: 1 /2nd: 2).
% =========================================================================
% Output: 
%       - shape: the shape name (string). e.g. 'CIRCLE' 
%       - x, y: the "ideal" (x,y) interpolated reference shape screen
%               coordinates.
%       - t_interpol: the correspdoning interpolated timestamp.
%       - x1, y1: the drawn (x,y) touch screen coordinates.
%       - t: the correspdoning timestamp.
%       - xref, yref: the (x,y) reference way-point coordinates.
%--------------------------------------------------------------------------
% Reference: 
% [1] Creagh, A.P., Simillion, C., Scotland, A., Lipsmeier, F., Bernasconi,
% C., Belachew, S., van Beek, J., Baker, M., Gossens, C., Lindemann, M. and
% De Vos, M., 2020. Smartphone-based remote assessment of upper extremity
% function for multiple sclerosis using the Draw a Shape Test.
% Physiological measurement, 41(5), p.054002.
%% Andrew Creagh. andrew.creagh@eng.ox.ac.uk
%  Last modified on Sept. 2017
%--------------------------------------------------------------------------
%% Initialise Variables
t_interpol = 0;
shape=shapeData{index,3};
%% Extract Shape Data
% Quality control,return if shape data empty or attempt doesn't exist
if isempty(shapeData{index,1}) && isempty(shapeData{index,2}) || ...
        isempty(shapeData{index,1}) &&  exist('attempt', 'var') && attempt==1 || ...
        isempty(shapeData{index,2}) &&  exist('attempt', 'var') && attempt==2

        x1=[]; y1=[]; t=[]; x=[]; y=[]; xref=[]; yref=[];
        return 
%otherwise extract touch-screen coordinates       
elseif exist('attempt', 'var') && attempt > 0
    %drawn touch screen x-coordinates   
    x1= getcolumn(shapeData{index, attempt}(:,1:2),1);
    %drawn touch screen y-coordinates
    y1=getcolumn(shapeData{index, attempt}(:,1:2),2);
    %corresponding timestamps
    t=getcolumn(shapeData{index, attempt}(:,3),1);
else
   %(if no explicit instructions given for which attempt, follow a series 
   % of rules to get the best attempt)
   
   %if all waypoints not completed (False) 
   if  strcmp(shapeData{index,4},'False') 
           %1) if number of attempts equals 1
           if shapeData{index,5}==1    
               %extract touch-screen coordinates
               x1= getcolumn(shapeData{index, 1}(:,1:2),1);
               y1=getcolumn(shapeData{index, 1}(:,1:2),2);
               t=getcolumn(shapeData{index, 1}(:,3),1);
               
           %2) if number of attempts equals 2
           elseif shapeData{index,5}==2
               %and the length of the first shape drawn is more than the
               %second shape.
               if length(shapeData{index, 1}) > length(shapeData{index, 2})
                    %extract touch-screen coordinates from the first shape
                    x1= getcolumn(shapeData{index, 1}(:,1:2),1);
                    y1=getcolumn(shapeData{index, 1}(:,1:2),2);
                    t=getcolumn(shapeData{index, 1}(:,3),1);
               else 
                    %extract touch-screen coordinates from the second shape
                    x1= getcolumn(shapeData{index, 2}(:,1:2),1);
                    y1=getcolumn(shapeData{index, 2}(:,1:2),2);
                    t=getcolumn(shapeData{index, 2}(:,3),1);
               end 
           %3) if number of attempts equals zero (0), return empty touch points 
           elseif shapeData{index,5}==0
                x1=[]; y1=[]; t=[]; x=[]; y=[]; xref=[]; yref=[];
               return
           end 
     
    else %if all waypoints completed (True) 
            %1) if number of attempts equals 1
            if shapeData{index,5}==1
                   %extract touch-screen coordinates
                   x1= getcolumn(shapeData{index, 1}(:,1:2),1);
                   y1=getcolumn(shapeData{index, 1}(:,1:2),2);
                   t=getcolumn(shapeData{index, 1}(:,3),1);
            %2) if number of attempts equals 2
            elseif shapeData{index,5}==2
                    %extract touch-screen coordinates
                    x1= getcolumn(shapeData{index, 2}(:,1:2),1);
                    y1=getcolumn(shapeData{index, 2}(:,1:2),2);
                    t=getcolumn(shapeData{index, 2}(:,3),1);
            end
   end
end

%take the timestamp point(s) from the first to have a monotonically
%increasing time vector
t=t-t(1);
t=t';

%% SHAPE INTERPOLATION
%--------------------------------------------------------------------------
%Interpolation Parameters:
%   m: the number of interpolated samples; 
%       - can be expressed as the number of samples: for m>) 
%       or 
%       - the step-size parameter: for m<1 (e.g. 1:0.5:10), where m=0.5
%--------------------------------------------------------------------------
if strcmp(mtype, 'length')  
    %interpolate to the length of drawn touch screen points
    m=length(x1); 
elseif strcmp(mtype, 'fixed')
    %interpolate to a fixed number of touch screen points
    m=0.05;    
    %m=100;
end 
% =========================================================================
%Perform interpolation:
if strcmp(shape, 'SPIRAL') || strcmp(shape, 'FIGURE_8') || strcmp(shape, 'CIRCLE')    
    %use spline interpolation in polar coordinates for circular-based shapes 
    interp_type='spline';
    [x, y] = interpolate_polar(ref_shape, interp_type,m, 0);
    xref = ref_shape(:,1); yref = ref_shape(:,2);
    xref(isnan(xref))=[];  yref(isnan(yref))=[];
    
elseif strcmp(shape, 'SQUARE') || strcmp(shape, 'LINE_TOP_TO_BOTTOM') || ...
       strcmp(shape, 'LINE_BOTTOM_TO_TOP')
    %use linear interpolation in cartesian coordinates for line-based shapes     
    interp_type='linear';  
    [x, y] = interpolate_cartesian(ref_shape, interp_type, m, 0);
    xref = ref_shape(:,1); yref = ref_shape(:,2);
    xref(isnan(xref))=[];  yref(isnan(yref))=[];
    
else %no interpolation
    x = ref_shape(:,1); y = ref_shape(:,2); % x-y ref coordinates 
    x(isnan(x))=[];
    y(isnan(y))=[];
    
end
%QC: remove nans
xref(isnan(xref))=[];
yref(isnan(yref))=[];
end
%EOF