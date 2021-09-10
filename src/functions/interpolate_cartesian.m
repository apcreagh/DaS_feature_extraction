function [X, Y] = interpolate_cartesian(ref_shape, interp_type, m, extra_options)
%% Draw-a-Shape: Interpolation (Cartesian)
%Perform interpolation of x- and y-coordinate reference way-points in
%cartesian coordinates.
%--------------------------------------------------------------------------
% Input:
%       - ref_shape: [N x 2] matrix of x- and y-coordinate reference way-point
%                    coordinates. 
%       - interp_type: interplotion method, e.g. 'linear' / 'spline'
%                      see interpn.m for more details.
%       - m: the number of interpolated samples; 
%            - can be expressed as the number of samples: for m>) 
%           or 
%            - the step-size parameter: for m<1 (e.g. 1:0.5:10), where m=0.5
%--------------------------------------------------------------------------
% Optional:
%       - extra_options.plot_transformation: bool (0/1) to plot the interpolation.
% =========================================================================
% Output: 
%       - X, X: the "ideal" (x,y) interpolated reference shape screen
%               coordinates.
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
%% Initialisation
x = ref_shape(:,1); y = ref_shape(:,2);
x(isnan(x))=[];
y(isnan(y))=[];
 
plot_transformation=false;
if exist('extra_options', 'var')
    if isfield(extra_options,'plot_transformation')
        plot_transformation=extra_options.plot_transformation;
    end 
end
%% Interpolation Parameterisation
pX=linspace(1,length(x),length(x));
pY=linspace(1,length(y),length(y));

if m > 1
    xq=linspace(1,length(x), m);
    yq=linspace(1,length(y), m);
else 
   xq=(1:m:length(x));
   yq=(1:m:length(y));
end 
%% Interpolation
X = interpn(pX,x, xq, interp_type);
Y = interpn(pY,y, yq, interp_type);
%% Visual Inspection
%plot interpolation
if plot_transformation
        figure
        plot(X,Y, '+')
        hold on
        plot(pX,x,'o',xq,X,'-');
        legend('Samples','Interpolation');

        figure
        plot(pX,x,'ko-')
        hold on

        figure
        plot(pY,y,'ko-')
        hold on
        plot(pY,y,'o',yq,Y,'-');
        legend('Samples','Interpolation');
end  
X=X'; Y=Y';
end
%EOF