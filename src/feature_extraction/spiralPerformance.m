function [features, feature_names]=spiralPerformance(xyref, xxyy, xy, extra_options)
%% Draw-a-Shape Spiral-Specific Feature Exraction
% Function to run the feature extraction for spiral shapes outlined in [1].
%--------------------------------------------------------------------------
% Input:
%       - xyref: [Nx2] matrix of reference shape x- and y-coordinates;
%       - xxyy:  [Nx3] matrix of interpolated reference shape x- and y-coordinates
%       - xy:    [Nx3] matrix of drawing x- and y-coordinates;
%                (x=[:,1],y=[:,2], t=[:,3])
% _________________________________________________________________________
%    extra_options: structure containing optional inputs to be used in each
%    feature extraction function. See specific extraction functions for
%    specific functionality.
%         - extra_options.sub_id: subject id (string, or integer)
% =========================================================================
% Output: 
%    features: a [1xN] vector of feature values 
%    feature_names:  a [1xN] string vector of corresponding feature labels.
%--------------------------------------------------------------------------
% Reference: 
% [1] Creagh, A.P., Simillion, C., Scotland, A., Lipsmeier, F., Bernasconi,
% C., Belachew, S., van Beek, J., Baker, M., Gossens, C., Lindemann, M. and
% De Vos, M., 2020. Smartphone-based remote assessment of upper extremity
% function for multiple sclerosis using the Draw a Shape Test.
% Physiological measurement, 41(5), p.054002.
%
%% Andrew Creagh. andrew.creagh@eng.ox.ac.uk
%  Last modified on Sept. 2018
%--------------------------------------------------------------------------%
%% Initialisation
feature_names={'AUC', 'GOF',  'd_bf', 'rmse', 'rmse_fit',  'rmse_ratio'};
features=NaN(1, length(feature_names));
if nargin < 1
    return
end 
%% Extract Shape Data
%Extract reference coordinates, interpolated reference coordinates and
%drawing touch-screen coordinates 
xref=xyref(:,1);    yref=xyref(:,2);  
x=xxyy(:,1);        y=xxyy(:,2);      t=xxyy(:,3);
x1=xy(:,1);         y1=xy(:,2) ;      t1=xy(:,3);
%--------------- Spiral Specific Quality Control -------------------------%
%If subject started at the wrong end  - drew from outside -> inside
%flip elements around
if CalcDistance(x1(1), y1(1), xref(1), yref(1))>...
        (CalcDistance(xref(1), yref(1), xref(end), yref(end)))*0.75
    %if the distance from the starting point in realtion to the distance
    %from the start of the ref to end of the ref is greater than arbitarty
    %level (chosen that if true the spiral probably is drawn opposite way
    %around), flip the elements around....
    x1=flipud(x1);
    y1=flipud(y1);
end 
%_________________________________________________________________________%

%% POLAR TRANSFORMATION
%-------------- polar coordinate transformation --------------------------%
x2=x1(:)-x1(1); %normalise
y2=y1(:)-y1(1); %normalise
% xref1=xyref(:,1)-xyref(1,1);
% yref1=xyref(:,2)-xyref(1,2);

%plot normalisation
% figure
% plot(x1(1:200), y1(1:200), 'r')
% hold on 
% plot(x1(200:end), y1(200:end), 'b')

% Transform Cartesian x- & y-coordinates to polar coordinates
[THETA,RAD] = cart2pol(x2,y2);
TH=THETA; R=RAD;

%plot transformation
% figure
% polarplot(TH,R, '.')
% hold on

%define cut-off value
cutoff=pi/6; % see unwrap.m
% unwrap phase angle
TH=unwrap(TH, cutoff);

%find the crossing points
[idx0,~] = crossing(TH);
%disregard...
TH(idx0(1):idx0(end))=[];
R(idx0(1):idx0(end))=[];

%normalise Radians
normR = R - min(R(:));
normR = normR ./ max(normR(:)); 

normTH = TH - min(TH(:));
normTH = normTH ./ max(normTH(:)); 

%% Feature Calculation
%-------------------------------------------------------------------------%
%              Drawing Error Metrics (Area-Based)
%-------------------------------------------------------------------------%
% create dummy vector...
refx=linspace(0,1,length(normTH));
refy=linspace(0,1,length(normTH));
%Trapezoidal numerical integration for each coordinate axis against itself;
%see calPerformance.m for further details...
dTH=trapz(normTH, refx);
dR=trapz(normR, refy);
AUC=dTH+dR;
%-------------------------------------------------------------------------%
%              Drawing Error Metrics (RMSE-Based)
%-------------------------------------------------------------------------%
%determine a linear fit between (unwraped) phase angle and radians
lm = fitlm(TH, R,'linear');
%get the RMSE of the fit (using matlabs inbuild function)
rmse=lm.RMSE;
%--> alternative:
%    f = fittype('a*x');
%    fit1 = fit(R,TH,f);

%plot the fit...
% figure
% plot(lm)
% xlabel('Theta (Radians)')
% ylabel('Radius')

%-------------- repeat for reference coordinates -------------------------%
xx=x(:)-x(1);
yy=y(:)-y(1);

[Tref,Rref] = cart2pol(xx, yy);
TH=Tref; R=Rref;
TH=unwrap(TH, cutoff);
[idx0,~] = crossing(TH);
TH(idx0(1):idx0(end))=[];
R(idx0(1):idx0(end))=[];
lm = fitlm(TH, R,'linear');
RMSE_ref=lm.RMSE;

% figure
% plot(lm)
% xlabel('Theta (Radians)')
% ylabel('Radius')

%determine "goodness of fit" by the ratio of the drawing rmse to reference
%shape rmse
GOF=rmse/RMSE_ref;

%-------------------------------------------------------------------------%
%              Emperical Archimedean Spiral Error Metrics (Fit-Based)
%-------------------------------------------------------------------------%
%Calculate the error metrics based on an emperical fit of the spiral
%drawing to a perfect (archimedean) spiral Note: this package will only
%work with spirals which have been constructed based on geometric spiral
%reference coordinates, and will not work with FLOODLIGHT PoC.

centre=NaN;
%reference data...
[best_center, rmse_refdata_fit]= fit_ArchimedeanSpiral(xref, yref, centre); 
%drawing touch-screen data...
[~, rmse_fit]= fit_ArchimedeanSpiral(x1, y1, centre);
%determine the ratio between drawing and reference archimedean sprial rmse
rmse_ratio=rmse_fit/rmse_refdata_fit;

%determine the absolute distance between the drawing best center and the first
%coordinate (which should be start at the "centre")...
[d_bf]= CalcDistance(best_center(1), best_center(2), x1(1), y1(1));
d_bf=abs(d_bf);

%% Save Features 
features=[AUC, GOF,  d_bf, rmse, rmse_fit,  rmse_ratio];
if length(features)~=length(feature_names)
    error('feature mismatch...'); end 
end 
%EOF