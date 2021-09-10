function [heatmap, max_pixel_time, n]=calcHeatmapFeatures(imageData, time_stamps, extra_options)
%% Draw-a-Shape Spatio-Temporal Heatmap Generation
% Function to run the spatio-temporal (heatmap-based) feature generation
% outlined in [1].
%--------------------------------------------------------------------------
% Input:
%       - imageData: [Nx2] matrix of reference shape x- and y-coordinates;
%                (x=[:,1],y=[:,2])
%       - time_stamps: the corresponding monitonically increasing timestamps
% _________________________________________________________________________
%    extra_options: structure containing optional inputs to be used in each
%    feature extraction function. See specific extraction functions for
%    specific functionality.
%      - extra_options.sub_id: subject id (string, or integer)
% =========================================================================
% Output: 
%      - heatmap: the raw heatmap image values
%      - max_pixel_time:  the maximum time spent in one region [ms]
%      - n: the pixel bin counts
%--------------------------------------------------------------------------
%% Initialisation
figpos=[3.7604    5.9375    4.1458    3.6667];
hist_bin=40; % set default is to 20x20 bins
grid_method='fixed';
ymin=600;  ymax=2200; xmin=0;  xmax=1450; 
fs=0;

if isfield(extra_options, 'grid_method')
    grid_method=options.grid_method; end
if isfield(extra_options, 'hist_bin')
    hist_bin=options.hist_bin; end
if isfield(extra_options, 'figpos')
    figpos=options.figpos; end

if ~isempty(time_stamps)
    time=abs(time_stamps(1)-time_stamps(end));
    fs=1/median(diff(time_stamps));
    pixel=1/fs; %the value each pixel represents in time
    pixel_max=max(diff(time_stamps)); %the maximum value each pixel represents in time
end
%% Pre-Proceesing
%fix image scale, pad image
padding=[xmin, ymin; xmin, ymax; xmax, ymin; xmax, ymax];
imageData=[padding; imageData];
%% Heatmap Calculation

%plot drawn shape
heatmap_fig=figure;
n = hist3(imageData, [hist_bin, hist_bin]); 
%set hist padding
n1 = n'; n1(size(n,1) + 1, size(n,2) + 1) = 0;
%set the pcolor grid either: 
% (1) a fixed shape/size 
% (2) set the the drawn touch coordinates
if strcmp(grid_method, 'fixed') 
    yb = linspace(ymin,ymax,size(n,1)+1);
    xb = linspace(xmin,xmax,size(n,1)+1);
else
    xb = linspace(min(imageData(:,1)),max(imageData(:,1)),size(n,1)+1);
    yb = linspace(min(imageData(:,2)),max(imageData(:,2)),size(n,1)+1);
end 

%convert to 2D color blocks based on bincounts
h = pcolor(xb,yb,n1);
colormap(hot) % heat map
caxis([0,max(max(n))]) %fix the color intensities 
grid on
%view(3);
%set the color interpolation parameters
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto'); 
set(gca,'XTickLabel','', 'YTickLabel','','XTick', [],'YTick', [])
set(gca,'visible','off');
heatmap_fig.Units='inches'; %denote units and keep fixed 
heatmap_fig.WindowState='normal'; %'pixels'
heatmap_fig.Position=figpos;

%get the raw image values
F = getframe(heatmap_fig);
imData=F.cdata;
heatmap=rgb2gray(imData);
close(heatmap_fig)

%% Feature Calculation
%determine the maximum time spent in one region 
max_pixel_time=max(max(n))/fs;
end 
%EOF