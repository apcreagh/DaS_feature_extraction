function [features, feature_names] = featureExtraction(shapeData, shape, extra_options)
%% Draw-a-Shape Temporal and Spatial Feature Exraction
% Shell function to run the feature extraction outlined in [1]. Extraction
% is performed on groups of functions characterising various feature
% domains.
%--------------------------------------------------------------------------
% Input:
%       - shapeData: a cell containing the shape data and information:
%          - shapeData{:, 1:2}: [Nx3] matrix of drawing x- and y-coordinates, 
%                               along with corresponding monotonically
%                               increasing timestamp vector. 
%               * shapeData{:, 1}: first attempt
%               * shapeData{:, 2}: second attempt
%                                  (if applicable, empty otherwise)
%          - shapeData{:, 3}: shape name (string, CAPS). e.g. 'CIRCLE'
%          - shapeData{:, 4}: bool (true/flase) flag to determine if the 
%                             shape was successfuly completed
%                             (as determined by the floodlight app.)
%          - shapeData{:, 5}: (int) the number of attempts
%       - shape: the shape name (string, CAPS). e.g. 'CIRCLE'
%       - referenceData (legacy, redundant): [N x 2] matrix of reference x- 
%                   and y-coordinate touch screen data (corresponding to 
%                   shapes drawn) taken from healthy subject (for 
%                   development of features);
%       - trgt_ptr_threshold:
%       - index: the shape index 
% _________________________________________________________________________
%    extra_options: structure containing optional inputs to be used in each
%    feature extraction function. See specific extraction functions for
%    specific functionality.
%         - extra_options.sub_id: subject id (string, or integer)
%         - extra_options.test_id: test id (string or datetime)
%         - extra_options.phoneType (legacy, redundant): the phone type,
%           'string' e.g. iPhone, SamsungsS6
%         - extra_options.mtype: see extractReferenceCoordinates.m
%         - extra_options.calculate_image_features: bool (true/false),
%           calculate image-based features
% =========================================================================
% Output: 
%    features: a [1xN] vector of feature values 
%    feature_names:  a [1xN] string vector of corresponding feature labels.
% ------------------------------------------------------------------------%
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
%% Example Usage
% (Best practice for creating a scaleable feature extraction)
% It is best practice when budiling modular feature extraction code to
% embed feature extraction functions in try/catch statements. If the file
% fails to execute or if there is an error associated with the feature
% extraction, the whole codebase won't fail, just the file you were
% extracting. We can then execute the same function in the catch statement
% to output a NaN feature file by calling the function with no argument.
%Example:
%--------------------------------------------------------------------------%
% try
%     [features, feature_names] = featureExtraction(shapeData, shape, extra_options);
% catch 
%     [features, feature_names] = featureExtraction();
% end
%_________________________________________________________________________%
%% Initialise Features
%Initialse all shape specific features (as NaNs). Calling with no argument
%outputs NaN features and the corresponding feature names. 
[speed_features, speed_feature_names]=calcTemporalFeatures();
speed_feature_names=strcat('speed_', speed_feature_names);
[aDS_features, aDS_feature_names]=calcTemporalFeatures();
aDS_feature_names=strcat('aDS_', aDS_feature_names);
[error_features, error_feature_names]=calcErrorFeatures();
[frequency_features, frequency_feature_names]=calcSpectralFeatures();
[circular_features, circular_feature_names]=circularPerformance();
[square_features, square_feature_names]= squarePerformance();
[spiral_features, spiral_feature_names]=spiralPerformance();
[figure8_features, figure8_feature_names]= figure8Performance();
[image_features, image_feature_names] = calcImageFeatures();

temporal_feature_names=[speed_feature_names, aDS_feature_names];
temporal_features=[speed_features, aDS_features];

feature_names=[temporal_feature_names,frequency_feature_names,...
error_feature_names, circular_feature_names, spiral_feature_names,...
square_feature_names, figure8_feature_names, image_feature_names];

features=[...
temporal_features, frequency_features, error_features,...
circular_features, spiral_features, square_features, figure8_features, ...
image_features];

if nargin < 1
    return
end
%% Options Initialisation
save_pathname=strcat(pwd, '/features/');
sub_id='MS0000';
test_id='30-09-2016';
mtype='length';%or: mtype='ground truth'; see extractReferenceCoordinates.m

calculate_image_features=true;

if isfield(extra_options, 'save_pathname')
    save_pathname=extra_options.save_pathname;end
if isfield(extra_options, 'sub_id')
    sub_id=extra_options.sub_id;end
if isfield(extra_options, 'test_id')
    test_id=extra_options.test_id;end
if isfield(extra_options, 'mtype')
    mtype=options.mtype;end
if isfield(extra_options, 'calculate_image_features')
    calculate_image_features=extra_options.calculate_image_features;end

if ~isfolder(save_pathname)
    mkdir(save_pathname)
    cd(save_pathname) 
end

warning('off','all');  % Turn off warning (be careful!)

%% Extract Shape Data
index=find(strcmp(shapeData(:, 3), shape));

%extract x- and y-touch coordinates & shape reference coordinates
[~, x, y, x1, y1, t, xref, yref]=extractReferenceCoordinates(shape, shapeData, index, mtype, extra_options);

%Pre-processing, filtering of touch-screen coordinates
[x1, y1, t, remove_index]=removeErroneousTouchPoints(x1, y1, t, 0.5); 
x(remove_index)=[];y(remove_index)=[];

%add to extra option (legacy issues)
extra_options.shape=shape;
extra_options.shapeData=shapeData;
extra_options.index=index;
extra_options.remove_index=remove_index;
%%
%get sampling freqency
fs=1/nanmedian(diff(t));
extra_options.fs=fs;

%drawing time
time=t(end); 

 
%get number of touch points in each shape
nx=length(x);
nx1=length(x1);
nxref=length(xref); 

% --- aside ---:
% I = location(1,1); J = location(1,2); % key to access desired shape
% x1 = getcolumn(shapeData{I, J}(:,1:2),1); %access desired shape x-,y-,t- coordinates
% y1 = getcolumn(shapeData{I, J}(:,1:2),2);

%% TEMPORAL FEATURES
%-------------------------------------------------------------------------%
%                          Temporal Metrics
%-------------------------------------------------------------------------%
t=t';t1=t;
%instantaneous velocity / speed:
vel=sqrt(diff(x1).^2 + diff(y1).^2)./(t(2:end));
speed=sqrt(diff(x1).^2 + diff(y1).^2)./diff(t);

%QC: remove NaNs and INFs
%[vel, t_vel]=removeNaNs_DS(vel, t);
[speed, t1]=removeNaNs_DS(speed, t1);
     
%compute moments of the speed signal:
[speed_features, speed_feature_names]=calcTemporalFeatures(speed, t1, extra_options);
speed_feature_names=strcat('speed_', speed_feature_names);
 
aDS=abs(speed);
[aDS_features, aDS_feature_names]=calcTemporalFeatures(aDS, t1, extra_options);
aDS_feature_names=strcat('aDS_', aDS_feature_names);

temporal_features=[speed_features, aDS_features];
temporal_feature_names=[speed_feature_names, aDS_feature_names];

% -------------------- plotting temporal features ------------------------%
% Plot instantaneous velocity / speed
% RGB=[166,206,227; 31, 120, 180; 123, 192, 62]./255;
% idx=1;
% raw_fig=figure;
% refPlot=plot(x, y,':k','LineWidth',2.5);
% hold on
% refMarkers=plot(xref, yref,'ko','LineWidth',2.5);
% hold on
% pDrawing=plot(x1,y1,'o', 'LineWidth',1.5, 'MarkerSize', 10);
% pDrawing.Color=RGB(idx, :);
% hold on
% axis equal tight
% set(gca,'XTick',[])
% hold on
% set(gca,'YTick',[])
% box off
% axis off
% axis square
% hold off;

% filename=strcat(extra_options.save_pathname, 'raw_Plot_',...
% extra_options.group_name,'_',extra_options.subject_name, '_',
% strrep(extra_options.test_date, ' ', '_'), '_',extra_options.shape);
% saveas(raw_fig,filename, 'fig') saveas(raw_fig,filename, 'png')
% saveas(raw_fig,filename, 'svg')


% fc=5; n=4; % cut off freq. @ 5Hz.; 4th order butterworth filter
% speed=filterData([speed, t1(2:end)], fc, n);   
% distance=(sqrt(diff(x1).^2 + diff(y1).^2));
% 
% fs=1/median(diff(t(:,end)));
% fc=5; n=4;
% Wn=fc/fs;
% [B,A] = butter(n, Wn);
% distance = filtfilt(B, A, distance);
% 
% 
% figure;
% plot(t(2:end), distance, 'k', 'LineWidth', 0.75)
% hold on
% xlabel ('Time [s]')
% ylabel ('Absolute Drawing Speed  (Pixels / s)')
% set(gca,'ytick',[])
% hold on

      
% fs=round(1/nanmedian(diff(t)));        
% fc=10; n=4;
% Wn=fc/fs;
% [B,A] = butter(n, Wn);
% speed = filtfilt(B, A, speed);
% speed=speed(:, 1);

% DS_fig=figure;
% plot(t1(2:end), ([speed(:)]), 'k', 'LineWidth', 1)
% hold on
% xlabel ('Time [s]')
% ylabel ({'Drawing Speed (Pixels / s)'})
% ylim([0, 1.1])
% ax=gca;
% ax.YAxis.Exponent =3;
% set(gca,'XMinorTick','on','YMinorTick','on')
% hold on
% set(gca, 'fontsize', 16)

% filename=strcat(extra_options.save_pathname, 'dsplot_',
% extra_options.group_name,'_',extra_options.subject_name, '_',
% strrep(extra_options.test_date, ' ', '_'), '_',extra_options.shape);
% saveas(DS_fig,strcat(filename, '.png')) saveas(DS_fig,strcat(filename,
% '.fig')) close(DS_fig)
        
% speed=detrend(speed, 'constant');  
% speed=rescale(speed, 0, 1);
% 
% aDS_fig=figure;
% plot(t1(2:end), abs([speed(:)]), 'k', 'LineWidth', 1)
% hold on
% xlabel ('Time [s]')
% ylabel ('| Drawing Speed |')
% ylim([0, 1.1])
% set(gca,'XMinorTick','on','YMinorTick','on')
% hold on
% set(gca, 'fontsize', 16)
% 
% filename=strcat(extra_options.pathname, 'aDSPlot_',
% extra_options.group_name,'_',extra_options.subject_name, '_',
% strrep(extra_options.test_date, ' ', '_'), '_',extra_options.shape);
% saveas(aDS_fig,strcat(filename, '.png')) saveas(aDS_fig,strcat(filename,
% '.fig')) close(aDS_fig)
%________________________ end plotting features __________________________%

%-------------------------------------------------------------------------%
%                          Jerk Metrics (removed [08-05-18])
%-------------------------------------------------------------------------%
% [b,a] = butter(4,0.6);
% % freqz(b,a)
% dataIn = vel;
% dataOut = filter(b,a,dataIn);
% vel=dataOut;
% % accel= diff(vel(1:end))./dt(2:end-1);
% % dt=diff(t_vel)';
% accel= diff(vel)./dt(2:end);
% 
% [b,a] = butter(4,0.6);
% % freqz(b,a)
% dataIn = accel;
% dataOut = filter(b,a,dataIn);
% accel=dataOut;
% 
% jerk= diff(accel)./dt(3:end);
% 
% [b,a] = butter(4,0.6);
% % freqz(b,a)
% dataIn = jerk;
% dataOut = filter(b,a,dataIn);
% jerk=dataOut;
% 
% %removed [08-05-18]
% % jerk_m=nanmean(jerk);
% % jerk_sd=nanstd(jerk);
% % jerk_max=max(jerk,'omitnan');
%% Spetral Freatures
%Frequency-based features
[frequency_features, frequency_feature_names]=calcSpectralFeatures(speed, t1, extra_options);
%% Error Features
% Error-based features
[error_features, error_feature_names]=calcErrorFeatures([xref, yref], [x,y, t], [x1,y1,t1], extra_options);
%% CIRCULAR SHAPE-SPECIFIC FEATURES
% calculate a sub-set of features based on a circular shapes
if strcmp(shape, 'CIRCLE') || strcmp(shape, 'FIGURE_8') || strcmp(shape, 'SPIRAL')
    %calculate circular-based featueres (seperate function circularPerformance.m
    [circular_features, circular_feature_names]=circularPerformance([xref, yref], [x,y, t], [x1,y1,t1], extra_options);
end
%% SPIRAL SHAPE-SPECIFIC FEATURES
if strcmp(shape, 'SPIRAL')
%calculate spiral-based featueres (seperate function spiralPerformance.m)
[spiral_features, spiral_feature_names]=spiralPerformance([xref, yref], [x,y, t], [x1,y1,t1], sub_id);
end 
%% SQUARE SHAPE-SPECIFIC FEATURES
if strcmp(shape, 'SQUARE')

%Determine the corner points of square
XY_MM=[min(x), max(y); min(x), min(y); ...
    max(x), min(y); max(x) max(y)]; 

%plot the corner points...
% figure
% plot(x, y, 'k.')
% hold on
% plot(XY_MM(:,1), XY_MM(:,2), 'ro')

% see squarePerformance.m
corners=XY_MM;
[square_features, square_feature_names]= squarePerformance([xref, yref], [x,y, t], [x1,y1,t1], speed, corners, [], extra_options);

end 
%% FIGURE-8 SHAPE-SPECIFIC FEATURES
if strcmp(shape, 'FIGURE_8')
[figure8_features, figure8_feature_names]= figure8Performance([xref, yref], [x,y, t], [x1,y1,t1], [], extra_options);
end 

%% Image-Specific Features
if calculate_image_features
[image_features, image_feature_names] = calcImageFeatures([xref, yref],[x,y, t], [x1,y1,t1], extra_options);
end 

%% Save Features
features=[...
temporal_features, frequency_features, error_features,...
circular_features, spiral_features, square_features, figure8_features, ...
image_features];

if length(features)~=length(feature_names)
    error('feature mismatch...'); end 
end
%EOF