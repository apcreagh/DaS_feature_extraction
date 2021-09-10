function [X, Y] = interpolate_polar(ref_shape, interp_type, m, plot_shape)
%% Draw-a-Shape: Interpolation (Polar)
%Perform interpolation of x- and y-coordinate reference way-points in polar
%coordinates. 
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
%       - plot_shape: bool (0/1) to plot the interpolation.
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
x(isnan(x))=[];  y(isnan(y))=[];

%convert (x,y) coordinates to polar coordinates
[theta,rho] = cart2pol(x,y);
%% Interpolation Parameterisation
pTheta=linspace(1,length(theta),length(theta));
pRho=linspace(1,length(rho),length(rho));
if m > 1
    tq=linspace(1,length(theta), m);
    rq=linspace(1,length(rho), m);
else 
    tq=(1:m:length(theta));
    rq=(1:m:length(rho));
end 
%% Interpolation
thetaq = interpn(pTheta,theta, tq, interp_type);
rhoq = interpn(pRho,rho, rq, interp_type);
%% Visual Inspection
%plot interpolation
if nargin>3 && plot_shape==1

    figure
    subplot(1,3,1)
    plot(pTheta,theta,'ko-')
    subplot(1,3,1)
    plot(pRho,rho,'ko-')
    subplot(1,3,3)
    plot(pRho,rho,'o',rq,rhoq,'-');
    legend('Samples','Spline Interpolation');

    figure
    plot(pTheta,theta,'o',tq,thetaq,'-');
    legend('Samples','Spline Interpolation');

    figure
    plot(pTheta,theta,'o',tq,thetaq,'-');
    legend('Samples','Spline Interpolation');

end 
%return to cart. coordinates
[X,Y] = pol2cart(thetaq,rhoq);
X=X'; Y=Y';
end
%EOF