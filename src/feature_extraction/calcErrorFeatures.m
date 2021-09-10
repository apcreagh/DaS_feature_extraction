function [features, feature_names]=calcErrorFeatures(xyref, xxyy, xy, extra_options)
%% Draw-a-Shape Error Feature Exraction
% Function to run the error (distance-based and area-based) feature
% extraction for outlined in [1].
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
%       - extra_options.trgt_ptr_threshold: the overlap threshold (in
%           pixels) between a touch screen coordinate and a reference
%           way-point prompted to the user, in order to be considered a 
%           "hit" (e.g. for celerity)
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
feature_names={'err_mean', 'err_sd','err_diff_cuml', 'error_time',...
'err25', 'err75', 'errIQR', 'areaError', 'mag_areaError', 'sqrtError',...
'sqrtErrorTime', 'drawPathRMSE', 'num_hits', 'acc_percent', 'celerity',...
'acc_celerity'};

features=NaN(1, length(feature_names));
if nargin < 1
    return
end 


% % % input the threshold size (radius) of pointer & target
trgt_ptr_threshold = [37.5 25];

if isfield(extra_options, 'trgt_ptr_threshold')
    trgt_ptr_threshold=extra_options.trgt_ptr_threshold;end

%% Extract Shape Data
xref=xyref(:,1);    yref=xyref(:,2);  
x=xxyy(:,1);        y=xxyy(:,2);      t=xxyy(:,3);
x1=xy(:,1);         y1=xy(:,2) ;      t1=xy(:,3);

time=abs(t1(end)-t1(1));

nx=length(x);
nx1=length(x1);
nxref=length(xref); 

%% Calculate Features
%-------------------------------------------------------------------------%
%              Drawing Error Metrics (Distance-Based)
%-------------------------------------------------------------------------%
% Caculate drawing error using N-D nearest point search:
% 1) calculate w.r.t. interpolated waypoints
[kNN, D1]=dsearchn([x1 y1], [x, y]);
% kNN=sort(kNN);
D1=D1/length(D1); %to normalise the distances

err_sum = nansum(D1); %error sum
% err_sum=err_sum/length(x1); %normalised version (removed [08-05-18])
err_mean = nanmean(D1); %error mean
err_sd = nanstd(D1); %error standard deviation
err_cuml = cumsum(D1)*length(x1); %cumulateive error
err_diff_cuml=mean(diff(err_cuml)); %mean difference in error between consecuative waypoints 
error_time=err_sum*length(x1); %cumulateive error w.r.t. time

% 2) calculate w.r.t. reference points
[kref, Dref]=dsearchn([x1 y1], [xref yref]); 
% kref=sort(kref);

%Drawing error over first 25% of shape vs. last 25% vs. overall
wp25=round(nx*0.25); % first 25% number of waypoints 
wp75=round(nx*0.75); % last 75% number of waypoints 
err25=sum(D1(1:wp25))/wp25;
err75=sum(D1(wp75:end))/wp25;

%interquartile error range
errIQR=iqr(D1)/length(D1);
%-------------------------------------------------------------------------%
%              Drawing Error Metrics (Area-Based)
%-------------------------------------------------------------------------%
%Trapezoidal numerical integration for each coordinate axis against itself
%zero-baseline
Ax=trapz(x1); Ay=trapz(y1); 
AUCshape=Ax+Ay; %simple AUC combination of the x- and y-coordinate 
%normalised AUC
AUCshape=AUCshape/length(x1);

%normalised AUC for each x- and y-coordinate individually
AxNorm=Ax/length(x1);
AyNorm=Ay/length(y1);

%reference coordinates AUC
Axref=trapz(x); Ayref=trapz(y); AUCref=Axref+Ayref;
AUCref=AUCref/length(x);

%absolute AUC, of x- w.r.t to y-coordinates 
AUCxy=abs(trapz(x1, y1));

%calculate the (area-based) error as the difference between the reference
%and drawn coordinate AUCs:
%1) releative difference
areaError=AUCshape-AUCref; 
%2) absolute difference
mag_areaError=abs(AUCshape-AUCref);
%3) root squared error difference
sqrtError=sqrt((AUCshape-AUCref).^2)/length(x1);
%4) root squared error difference w.r.t. drawing time
sqrtErrorTime=trapz(t1, sqrt(x1.^2 + y1.^2)/length(x1));
%5) RMSE between the raw drawing paths (interpolated coordinates)
%legacy....
[~, xinterp, yinterp, ~, ~, ~, ~, ~]=extractReferenceCoordinates(...
                        extra_options.shape, extra_options.shapeData, extra_options.index, 'length', extra_options);

remove_index=extra_options.remove_index;
xinterp(remove_index)=[]; yinterp(remove_index)=[];
drawPathRMSE=rmse(sqrt(xinterp.^2 + yinterp.^2), sqrt(x1.^2 + y1.^2));  
drawPathRMSE=drawPathRMSE/length(x1);

%-------------------------------------------------------------------------%
%              Drawing Error Metrics (Accuracy-Based)
%-------------------------------------------------------------------------%
hit = zeros(1,nxref);
%cycle through the shape and register a hit if the distance between a
%drawing touch point and the reference coordinate is within a target
%pointer threshold:
hit=Dref < max(trgt_ptr_threshold);

%number of successful hits
num_hits = sum(nonzeros(hit(1,:)));
% percent(%) accuracy  
acc_percent = (num_hits / nxref); 

hit=hit'; hit=logical(hit);
khit=kref(hit);

% %number of successful baseline shape hits(for comparison)
% num_hits_B = sum(nonzeros(hitBaseline(1,:)));
% acc_percent_B = (num_hits_B / nxref)*100; 
% accuracy_B = {hitBaseline, num_hits_B, acc_percent_B};
% hitBaseline=hitBaseline'; hitBaseline=logical(hitBaseline);
% khitH=krefH(hitBaseline);

%celerity
celerity=num_hits/time;
%accuracy based celerity
acc_celerity=acc_percent/time; 

% ------------------ plotting accuracy features -----------------------%
% figure;
% plot(x, y,':k','LineWidth',1 )
% hold on
% plot(xref, yref,'ko','LineWidth',1 )
% hold on
% plot(x1,y1,'ro', 'MarkerSize', 5)
% hold on
% plot(x1(khit), y1(khit), 'g*', 'MarkerSize', 8)
% hold on
% hline(yref(4), 'k:')
% vline(xref(4),'k:')
% title(['UEFG,  Subject # ', sub_id])
% axis equal tight
% axis ([-100 1100 -200 1800])
% set(gca,'XTick',[])
% hold on
% set(gca,'YTick',[])
% hold off;

% reflectiony=1000;
% raw_fig=figure;
% plot(x, 2*reflectiony-y,':k','LineWidth',3 )
% hold on
% plot(xref, 2*reflectiony-yref,'ko','LineWidth',3)
% hold on
% scatter(x1,2*reflectiony-y1,60, 'r', 'LineWidth',1)
% hold on
% scatter(x1(1:50),2*reflectiony-y1(1:50),60, 'b', 'LineWidth',1)
% hold on
% plot(x1(khit), 2*reflectiony-y1(khit), 'g.', 'MarkerSize', 15)
% % plot(x1(kref), 2*reflectiony-y1(kref), 'g*', 'MarkerSize', 7)
% %  plot(x1(kNN), 2*reflectiony-y1(kNN), 'g*', 'MarkerSize', 7)
% axis equal tight
% axis ([-100 1500 -200 1500])
% hold on
% set(gca,'XTick',[])
% hold on
% set(gca,'YTick',[])
% axis off;
% 
% filename=strcat(extra_options.pathname, 'rawplot_',...
% extra_options.group_name,'_',extra_options.subject_name, '_',
% strrep(extra_options.test_date, ' ', '_'), '_',extra_options.shape);
% saveas(raw_fig,strcat(filename, '.png'))
% saveas(raw_fig,strcat(filename,...
% '.fig')) close(raw_fig)

% xaxis=linspace(1,100,length(err_cuml));
% xaxis=xaxis';
% 
% figure
% plot(xaxis, err_cuml,'r','LineWidth', 2)
% hold on
% plot(xaxis, cum_error_BL, 'b','LineWidth', 2)
% hold on
% a=[cellstr(num2str(get(gca,'xtick')'))]; 
% % Create a vector of '%' signs
%      pct = char(ones(size(a,1),1)*'%'); 
% % Append the '%' signs after the percentage values
%      new_xticks = [char(a),pct];
% % 'Reflect the changes on the plot
%      set(gca,'xticklabel',new_xticks)
% xlabel('Percent of Shape Length [%]')
% ylabel('Cumulative Drawing Error')
% 
% tX=linspace(0,length(t),length(t));
% % tX=linspace(t(1),t(end),length(t));
% 
% time_D1 = linspace(t(1),t(end),length(D1));
% time_DHealthy = linspace(tBaseline(1),tBaseline(end),length(DBaseline));
% % 
% tX=t;
% 
% tX=tX';
% xq=1:1:(length(err_cuml)); %query points 
% xq=xq';
% % tq=interp1(t,tX, xq, 'linear');
% 
% tX=tX';
% xq=1:1:(length(err_cuml));
% xq=xq';
% tq=interpn(tX,t, xq, 'spline');
% 
% tB=linspace(tBaseline(1),tBaseline(end),length(tBaseline));
% tB=tB';
% xqb=1:1:(length(cum_error_BL));
% xqb=xqb';
% tq_BL=interp1(tBaseline,tB, xqb, 'linear');
% 
% taxis=linspace(1,100,length(err_cuml));
% xaxis=xaxis';
% 
% figure
% plot(time_D1, err_cuml,'r', 'LineWidth', 2)
% hold on 
% plot(time_DHealthy, cum_error_BL, 'b', 'LineWidth', 2)
% hold on
% xlabel('Time to Complete Shape [s]')
% ylabel('Cumulative Drawing Error')
% set(gca,'YTick',[])
% title(shape)
% 
% figure
% plot(tq, err_cuml,'r', 'LineWidth', 2)
% hold on 
% plot(tq_BL, cum_error_BL, 'b', 'LineWidth', 2)
% % plot(cum_error_BL, 'b', 'LineWidth', 2)
% % hold on
% xlabel('Time to Complete Shape [s]')
% ylabel('Cumulative Drawing Error')
% title(shape)
% 
% figure
% plot(err_cuml,'r', 'LineWidth', 2)
% hold on 
% plot(cum_error_BL, 'b', 'LineWidth', 2)
% % plot(cum_error_BL, 'b', 'LineWidth', 2)
% % hold on
% xlabel('Waypoint')
% ylabel('Cumulative Drawing Error')
% title(shape)
%________________________ end plotting features __________________________%
 
%% Save Features
features=[err_mean, err_sd, err_diff_cuml, error_time, err25, err75,...
errIQR, areaError, mag_areaError, sqrtError, sqrtErrorTime, drawPathRMSE,...
num_hits, acc_percent, celerity, acc_celerity];

if length(features)~=length(feature_names)
    error('feature mismatch...'); end
end 
%EOF