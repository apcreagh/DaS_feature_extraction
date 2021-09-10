function [features, feature_names]=squarePerformance(xyref, xxyy, xy, speed, corners, sub_id, extra_options)
%Draw-a-Shape Spiral-Specific Feature Exraction
% Function to run the feature extraction for spiral shapes outlined in [1].
%--------------------------------------------------------------------------
% Input:
%       - xyref: [Nx2] matrix of reference shape x- and y-coordinates;
%       - xxyy:  [Nx3] matrix of interpolated reference shape x- and y-coordinates
%       - xy:    [Nx3] matrix of drawing x- and y-coordinates;
%                (x=[:,1],y=[:,2], t=[:,3])
%       - corners:  the corner points of square, see calPerformance.m
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
feature_names={'polar_mag_error','polar_mag_error_time', 'dwellTime', 'dwellTime_mean', 'dwellTime_sd', 'dwell_to_move', 'overShoot','overShoot_mean','overShoot_sd'};
features=NaN(1, length(feature_names));
if nargin < 1
    return
end 
%% Extract Shape Data
xref=xyref(:,1);    yref=xyref(:,2);  
x=xxyy(:,1);        y=xxyy(:,2);      t=xxyy(:,3);
x1=xy(:,1);         y1=xy(:,2) ;      t1=xy(:,3);

extra_options.shape='SQUARE'; %legacy issues [22-08-19]

%% Calculate Features
%-------------------------------------------------------------------------%
%              Error Feature(s)
%-------------------------------------------------------------------------%
[~, ~, ~, ~, polar_mag_error]=calculatePolarVelocity(x1, y1, t, xref, yref);
%time-normalised magnitude error
time=abs(t(end)-t(1));
polar_mag_error_time=polar_mag_error/time; 
%-------------------------------------------------------------------------%
%              Temporal Metrics (Hesitation-Based)
%-------------------------------------------------------------------------%
%---- determine speed minimum and maximum, i.e. stopping points-----------%
speed=[0; speed];
speed=speed-mean(speed);
%speed=detrend(speed, 'linear');
mpd=length(speed)*0.15; %arbitary 15%, see calPerformance.m
mph=mean(speed)+std(speed); 
mph_inv=mean(speed)+std(speed);

%update [26-02-19]:
[mx,mini,  w_max, p_max] = findpeaks(speed,'MinPeakProminence',std(speed), 'minpeakheight', mph ,'minpeakdistance', mpd, 'WidthReference', 'halfprom');
speed_inv = 1.1*max(speed) - speed;
[mi,maxi, w_min, p_min] = findpeaks(speed_inv,'MinPeakProminence',std(speed),  'minpeakheight', mph_inv ,'minpeakdistance', mpd, 'WidthReference', 'halfheight');
mi = speed(mini);
w_max=floor(w_max/2); w_min=floor(w_min/2);

%-------------- determine dwell (hesitation) time ------------------------%
for k=1:min([length(mini), length(w_min)])
    thres(k,1)=mini(k)-w_min(k); %low threshold
    if thres(k,1) <1
       thres(k,1)=1;
    end
    thres(k,2)=mini(k)+w_min(k); %high threshold
    if  thres(k,1) > length(t1) || thres(k,2) > length(t1) 
    %--------------- Square Specific Quality Control ---------------------%
    %denote as corrupted signal
    [dwellTime, dwellTime_mean, dwellTime_sd, dwell_to_move, overShoot,overShoot_mean,overShoot_sd]=deal(NaN);
    else
    dwell_time(k)=abs(t1(thres(k,2))-t1(thres(k,1)));
    end 
end

%----------- moments of the dwell (hesitation) time ----------------------%
dwellTime=sum(dwell_time);
dwellTime_mean=mean(dwell_time);
dwellTime_sd=std(dwell_time);
dwell_to_move=dwellTime/t(end); %length of time stationary 

%-------------------------------------------------------------------------%
%              Error Metrics (Overshoot-Based)
%-------------------------------------------------------------------------%
%-------------------- overshoot at the corners ---------------------------%
x1y1=[x1(mini),y1(mini)];
[knn, overshoot]=dsearchn(x1y1, corners);
% knn=sort(knn);
[~,ia,~]=unique(knn, 'stable'); 
overshoot=overshoot(ia);
%remove outliers that are most likely errors...
overshoot((overshoot(:)>(mean(overshoot)+(std(overshoot)*3))))=[]; 

%overshoot metrics...
overShoot=sum(abs(overshoot)); 
overShoot_mean=mean(abs(overshoot));
overShoot_sd=std(abs(overshoot)); 

%% Feature Plotting
% tX=linspace(1,length(T),length(T));
% tX=tX';
% xq=1:1:(length(x1));
% xq=xq';
% tq=interpn(tX,T, xq, 'linear');
% 
% yq=1:1:(length(error));
% yq=yq';
% yq=interpn(tX,T, yq, 'linear');
% min_buffer_lead=round(w_min/2);
% min_buffer_trail=round(w_min/2);
% 
% max_buffer_lead=2;
% % max_buffer_lead=round(fs*0.25);
% max_buffer_trail=1;
% 
% reflectiony=1000;
% 
% Sq_fig=figure;
% plot(x, 2*reflectiony-y, 'k:', 'LineWidth',2)
% hold on
% plot(xref, 2*reflectiony-yref, 'ko', 'LineWidth',2)
% hold on
% scatter(x1,2*reflectiony-y1, 'bo','LineWidth',1)
% hold on
% for i=1:length(MaxIdx)
%         scatter(x1(max(MaxIdx(i)-max_buffer_lead, 1):MaxIdx(i)+max_buffer_trail), 2*reflectiony-y1(max(MaxIdx(i)-max_buffer_lead, 1):MaxIdx(i)+max_buffer_trail), 'go','filled', 'LineWidth',2 ); end
% hold on
% for j=1:length(MinIdx)
%         scatter(x1(max(MinIdx(j)-min_buffer_lead, 1):MinIdx(j)+min_buffer_trail), 2*reflectiony-y1(max(MinIdx(j)-min_buffer_lead, 1):MinIdx(j)+min_buffer_trail), 'ro','filled', 'LineWidth',2 ); end
% % % hold on
% % % plot(x1(1:round(fs*1)),2*reflectiony-y1(1:round(fs*1)), 'ko')
%  axis equal tight
% axis ([-100 1500 -200 1500])
% hold on
% set(gca,'XTick',[])
% hold on
% set(gca,'YTick',[])
% axis off;
% 
% filename=strcat(extra_options.save_pathname, 'sqplot_', extra_options.group_name,'_',extra_options.subject_name, '_', strrep(extra_options.test_date, ' ', '_'), '_',extra_options.shape); 
% saveas(Sq_fig,strcat(filename, '.png'))
% saveas(Sq_fig,strcat(filename, '.fig'))

%% Save Featues
features=[polar_mag_error, polar_mag_error_time, dwellTime, dwellTime_mean, dwellTime_sd, dwell_to_move, overShoot,overShoot_mean,overShoot_sd];
if length(features)~=length(feature_names)
    error('feature mismatch...'); end 
end
%EOF