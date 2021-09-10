function [features, feature_names]=figure8Performance(xyref, xxyy, xy, sid, extra_options)
%% Draw-a-Shape Figure-8-Specific Feature Exraction
% Function to run the feature extraction for figure-8-shapes outlined in [1].
%--------------------------------------------------------------------------
% Input:
%       - xyref: [Nx3] matrix of reference shape x- and y-coordinates;
%       - xxyy:  [Nx3] matrix of interpolated reference shape x- and y-coordinates
%       - xy:    [Nx3] matrix of drawing x- and y-coordinates;
%                (x=[:,1],y=[:,2], t=[:,3])
%       - sid (legacy, redundant): subject id (string)
% _________________________________________________________________________
%    extra_options: structure containing optional inputs to be used in each
%    feature extraction function. See specific extraction functions for
%    specific functionality.
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
feature_names={'asymX','asymY', 'centre_offshoot','TOTSYMM', 'areaMismatch'};
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

%% Calculate Features
%determine the minimum peak heights and distances
mpd=length(x1)*0.15; %arbitary 15%
mph=mean(x1)+std(x1)*0.5;
%get inverse signal...
sigxInv = 1.01*max(x1) - x1;
mph_inv=mean(x1)-std(x1);
MPP=std(x1)*0.75;

%determine the peak locations for all 4 sides of the square (i.e. positive
%and negative peaks
[~, locs]=findpeaks(x1, 'minpeakheight',  mph, 'minpeakdistance' ,mpd, 'MinPeakProminence',MPP);
[~, locs1]=findpeaks(sigxInv, 'minpeakheight',  mph_inv, 'minpeakdistance' ,mpd);

%create interpolation grids
xx=linspace(x1(1),x1(end), length(x1));xx=xx';
xxref=linspace(x1(1),x1(end), length(xref));xxref=xxref';
[X, Y]=interpolate_cartesian([xxref, xref], 'cubic', 0.01, extra_options);

%--------------- Figure-8 Specific Quality Control -------------------------%
%if the person started in the wrong directon, and is shifted upon the origin,flip signal
if locs1(1)<locs(1)
    x1=(-x1);
    x1=x1-mean(x1);
    xref=xref-mean(xref);
end 
 
%-------------------------------------------------------------------------%
%              Error Feature(s)
%-------------------------------------------------------------------------%
%the square-specific area between drawing and reference coordinates: 
%normalised area mismatch
areaMismatch = abs(trapz(X,Y)-trapz(xx,x1))/length(X);
%-------------------------------------------------------------------------%
%              Asymmetry Feature(s)
%-------------------------------------------------------------------------%

% Get distanxces from the centre 
xX=[xx,x1]; XY=[xxref(4),xref(4)];
[~, Dcentre]=dsearchn(xX,XY);
XY=[xxref(10),xref(10)];
[~, Dcentre(2)]=dsearchn(xX,XY);

%calculate center offshoot by the sum of the offshoots...
centre_offshoot=sum(Dcentre);

%---------- calculate x-coordinate (horrizontal) asymmetry %-------------%
asymX=zeros(2,2);
asymX(:,2)=abs(xref(4)-x1(locs)); %right side
asymX(:,1)=abs(xref(4)-x1(locs1)); %left side

%refactoring  x-/y- : left-/right-: horrizontal/vertical to match up...
if asymX(1,2)>asymX(1,1)
    asym_XTop=asymX(1,2)/asymX(1,1);
else 
     asym_XTop=asymX(1,1)/asymX(1,2);
end 
if asymX(2,2)>asymX(2,1) 
    asym_XBot=asymX(2,2)/asymX(2,1);
else 
    asym_XBot=asymX(2,1)/asymX(2,2);
end 

%---------- calculate y-coordinate (vertical) asymmetry %-------------%
[~,iMin] = min(y1);
[~,iMax] = max(y1); 
asymY(1)=abs(yref(4)-y1(iMin)); %bottom
asymY(2)=abs(yref(4)-y1(iMax)); %top

%refactoring  x-/y- : left-/right-: horrizontal/vertical to match up...
if asymY(1)>asymY(2) 
    asymY=asymY(1)/asymY(2);
else 
    asymY=asymY(2)/asymY(1); 
end 

%add to determine "total" asymmetry
asymX=asym_XTop + asym_XBot;
TOTSYMM= asymY + asymX;

%% Save Features
features=[asymX,asymY,centre_offshoot,TOTSYMM,areaMismatch];
if length(features)~=length(feature_names)
    error('feature mismatch...'); end 
end