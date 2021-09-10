function [REFERENCE_COORDINATE, x, y, x1, y1, t, xref, yref, extra_options]=extractReferenceCoordinates(choose_shape, shapeData, index, mtype, extra_options)
%% Draw-a-Shape: Function to Extract Drawing Touch Points
% Function to extract the drawing touch points (x- and y-coordinates) as
% well as the reference coordinates for that shape. 
%--------------------------------------------------------------------------
% Input:
%       - choose_shape: the shape name (string). e.g. 'CIRCLE'
%       - shapeData: [N x 2] matrix of x- and y-coordinate touch screen data 
%                    (corresponding to shapes drawn)
%       - shape: the shape name (string). e.g. 'CIRCLE'
%       - healthyData: [N x 2] matrix of refereence x- and y-coordinate touch
%                       screen data (corresponding to shapes drawn)
%       - index: the index of the shape to extract
%       - mtype: see shapePreperation.m file
% _________________________________________________________________________
%       - extra_options: structure containing optional inputs to be used in 
%                        each feature extraction function. 
%                        See specific extraction functions.
%       - extra_options.phone_type: the smartphone type e.g. Samsung 'S7',
%                        which can be used to call the correct reference 
%                        coordinates file. 
%       - extra_options.attempt: the drawing attempt to extract (1st: 1 /2nd: 2)
%       - extra_options.QC: bool (0/1) to perform drawing quality control (QC) 
%                       based on Hausdorff Distance (HD) threshold
%       - extra_options.HD_Tolerance: the HD threshold (in pixels)
% =========================================================================
% Output: 
%       - REFERENCE_COORDINATE: the reference coordinate filename used (string)
%       - x, y: the "ideal" (x,y) interpolated reference shape screen
%               coordinates.
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
%--------------------------------------------------------------------------%
%% Parameterisation
%set the reference coordinates as default SamsungS7
REFERENCE_COORDINATE='REFERENCE_COORDINATES_S7';
% REFERENCE_COORDINATE='REFERENCE_COORDINATES_S5';

%set defaults
attempt=0; % the shape attempt to pull;
HDT=450; %Hausdorf distance threshold for quality control check (hard-coded);
QC=false;

%load extra options parameters, if they've been called
if exist('extra_options', 'var')
    if isfield(extra_options, 'phone_type')  
        REFERENCE_COORDINATE=strcat('REFERENCE_COORDINATES_', extra_options.phone_type); end 
    if isfield(extra_options, 'attempt') % the shape attempt to pull
        attempt=extra_options.attempt; end
    if isfield(extra_options, 'HD_Tolerance') %Hausdorf Tolerance (QC)
        HDT=extra_options.HD_Tolerance; end 
    if isfield(extra_options, 'QC') %Hausdorf Tolerance (QC)
        QC=extra_options.QC; end 
end
%% Extract Shapes    
%load reference coordinates for specific model/OS of phone: e.g. SamsungS7
ref_coordinates = readtable([pwd, '/data/', REFERENCE_COORDINATE], 'Delimiter', '\t' ); 
%extract the reference coordinates corresponding to a specific shape
[ref_shape] = refShape(ref_coordinates, choose_shape);
%extract the x- and y-touch screen coordinates for ideal shape, shape
%drawn, and reference coordinates. See shapePreperation.m for more details.
[~, x, y, x1, y1, t,xref, yref, ~]=shapePreperation(shapeData, ref_shape, index,  mtype, attempt);
%% Quality Control 
% evaluate the hausdorff distance between ideal shape and draw shape.
% If hausdorff distance is empty, or below a threshold (HDT) discard the
% shape and return empty;
[hd, ~] = HausdorffDist([x1, y1],[x,y],[]);

if isempty(hd)
    return
end 

% QC: can be turned on or off with extra options
if hd > HDT && QC 
    x1=[]; y1=[]; t=[]; x=[]; y=[]; xref=[]; yref=[];
    return
end 

%% Visual Quality Control
%%plot the shape 
% reflectiony=1000;
% (refelection is to flip the shape the correct way it would have been seen
% by the user on the screen)
% figure
% plot(x, 2*reflectiony-y,'k:.','LineWidth',1 )
% hold on
% plot(xref, 2*reflectiony-yref,'ko','LineWidth',1 )
% hold on
% plot(x1,2*reflectiony-y1,'bo', 'LineWidth',0.5)
end 
%EOF