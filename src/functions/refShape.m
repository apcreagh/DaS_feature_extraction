function [ref_shape] = refShape(ref_coordinates, shape)
%% Draw-a-Shape: Function to match Shape Reference Coridnates 
% Function to match the shape reference coordinates 
%--------------------------------------------------------------------------
% Input:
%       - ref_coordinates: a table of reference coordinate way-points for 
%                          various shapes
%       - choose_shape: the shape name (string) to extract. e.g. 'CIRCLE'
% =========================================================================
% Output: 
%       - ref_shape:  [N x 2] matrix of reference coordinate way-points for
%                      the selected shape; stored as x- and y-coordinates 
%                                               (column 1 and column 2)
%--------------------------------------------------------------------------
% Reference: 
% [1] Creagh, A.P., Simillion, C., Scotland, A., Lipsmeier, F., Bernasconi,
% C., Belachew, S., van Beek, J., Baker, M., Gossens, C., Lindemann, M. and
% De Vos, M., 2020. Smartphone-based remote assessment of upper extremity
% function for multiple sclerosis using the Draw a Shape Test.
% Physiological measurement, 41(5), p.054002.
%
%% Andrew Creagh. andrew.creagh@eng.ox.ac.uk
%  Last modified on Sept. 2017
%--------------------------------------------------------------------------
% extract the reference coordinates corresponding to the inputted shape
if strcmp(shape,'CIRCLE')
    ref_shape(:,1) = ref_coordinates.CIRCLE_X;
    ref_shape(:,2) = ref_coordinates.CIRCLE_Y;
elseif strcmp(shape,'FIGURE_8') 
    ref_shape(:,1) = ref_coordinates.FIGURE_EIGHT_X;
    ref_shape(:,2) = ref_coordinates.FIGURE_EIGHT_Y;
elseif strcmp(shape,'LINE_TOP_TO_BOTTOM') 
    ref_shape(:,1) = ref_coordinates.LINE_TOP_TO_BOTTOM_X;
    ref_shape(:,2) = ref_coordinates.LINE_TOP_TO_BOTTOM_Y;
elseif strcmp(shape,'LINE_BOTTOM_TO_TOP')    
    ref_shape(:,1) = ref_coordinates.LINE_BOTTOM_TO_TOP_X;
    ref_shape(:,2) = ref_coordinates.LINE_BOTTOM_TO_TOP_Y;
elseif strcmp(shape,'SPIRAL')  
    ref_shape(:,1) = ref_coordinates.SPIRAL_X;
    ref_shape(:,2) = ref_coordinates.SPIRAL_Y;
elseif strcmp(shape,'SQUARE')   
    ref_shape(:,1) = ref_coordinates.SQUARE_X;
    ref_shape(:,2) = ref_coordinates.SQUARE_Y;
end
end
%EOF