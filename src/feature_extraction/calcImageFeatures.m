function [features, feature_names] = calcImageFeatures(xyref, xxyy, xy, extra_options)
%Draw-a-Shape Image-Based Feature Exraction
% Function to run the image-based feature extraction for outlined in [1].
%--------------------------------------------------------------------------
% Input:
%       - xyref: [Nx2] matrix of reference shape x- and y-coordinates,
%       - xxyy:  [Nx3] matrix of interpolated reference shape x- and
%                y-coordinates, and corresponding time-stamps;
%       - xy:    [Nx3] matrix of drawing x- and y-coordinates; and 
%                corresponding time-stamps;
%                (x=[:,1],y=[:,2],time=[:,3])
% _________________________________________________________________________
%    extra_options: structure containing optional inputs to be used in each
%    feature extraction function.
%       - extra_options.sub_id: subject id (string, or integer)
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
feature_names={'mIm','sdIM','imCorr','peaksnrIM','snrIM','ssimvalIM',...
'MSEerrIM','nMSE','nRMSE','MI_im','max_pixel_time','mHM',...
'sdHM', 'hmCorr', 'MSEerrHM', 'peaksnrHM', 'snrHM', 'ssimvalHM',...
'hmEntropy', 'hmEntropy_ratio', 'MI_hm', 'imRatio',...
'isdRatio', 'pDchisq', 'errorValue', 'iqrHausD', 'stdHausD', 'varHausD',...
'maxHausD', 'minHausD', 'skewHausD', 'kurtHausD', 'hausD_t',...
'randomError',  'HausD_25', 'HausD_75', 'hausDmiddle_t'};

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

%% Pre-Processing

%------------------Shape conversion to Images ----------------------------%
%Convert drawn touchpoints to image
figpos=[3.7604    5.9375    4.1458    3.6667];
fig1=figure;
plot(x1,y1, 'k.', 'MarkerSize',20,'MarkerFaceColor','auto')
% axis([-200, 1600, 0, 2500])
set(gca,'XTickLabel','', 'YTickLabel','','XTick', [],'YTick', [])
set(gca,'visible','off');
colormap gray
caxis([0, max([y1;x1])])
%set axis and position constant
fig1.Units='inches';
fig1.WindowState='normal'; %'pixels'
fig1.Position=figpos;
F = getframe(fig1);
imData1=F.cdata;
im1=rgb2gray(imData1);

%Convert reference touchpoints to image
fig2=figure;
plot(x,y, 'k.', 'MarkerSize',20,'MarkerFaceColor','auto')
% axis([-200, 1600, 0, 2500])
set(gca,'XTickLabel','', 'YTickLabel','','XTick', [],'YTick', []);
set(gca,'visible','off');
colormap gray
caxis([0, max([y1;x1])])
%set axis and position constant
fig2.Units='inches';
fig2.WindowState='normal'; %'pixels'
fig2.Position=figpos;
F = getframe(gcf);
imData2=F.cdata;
im2=rgb2gray(imData2);

close(fig1);close(fig2);

%----------------- Shape conversion to Heatmaps---------------------------%
% hm1=figure;
% h=dscatter(x1,y1, 'PLOTTYPE', 'surf', 'BINS',[800,800], 'SMOOTHING',10)
% view(0, 90)
% set(gca,'visible','off');
% F = getframe(gcf);
% imData=F.cdata;
% hm1=rgb2gray(imData);
% imshow(hm1)

% hm2=figure;
% h=dscatter(x,y, 'PLOTTYPE', 'surf', 'BINS',[800,800], 'SMOOTHING',10)
% view(0, 90)
% set(gca,'visible','off');
% F = getframe(gcf);
% imData=F.cdata;
% hm2=rgb2gray(imData);
% imshow(hm2)

%% Feature Calculation

%-------------------------------------------------------------------------%
%               Image-Based Similarity Metrics & Moments
%-------------------------------------------------------------------------%
%Compute measures between reference image and drawn image
%e.g. moments of the image
mIm=mean2(im1);
sdIM=std2(im1);
%2-D correlation
imCorr=corr2(im1, im2);
%Peak Signal-to-Noise Ratio
[peaksnrIM, snrIM] = psnr(im1, im2);
%Self-similarity metrics
[ssimvalIM,ssimmapIM] = ssim(im1,im2); 
%MSE, RMSE
MSEerrIM = immse(im1,im2);
nMSE  = (immse([x1, y1],[x, y]));
nRMSE = sqrt(immse([x1, y1],[x, y]));
% Compute the mutual information of two images
MI_im = MI(im1, im2); 
%-------------------------------------------------------------------------%
%               Heatmap-Based Similarity Metrics & Moments
%-------------------------------------------------------------------------%

%----------------- Generate Pixel Density Heat Maps ----------------------%
%downsample the touchpoints
yds = downsample(y,2);  xds = downsample(x,2); 
[hm1, max_pixel_time, n]=calcHeatmapFeatures([x1, y1], t, extra_options);
[hmref, ~, nref]=calcHeatmapFeatures([xds,yds], [], extra_options);

%--------- Calculate Drawn and Image Comparison Metrics ------------------%
%moments of the heatmap image
mHM=mean2(hm1);
sdHM=std2(hm1);
%2-D correlation
hmCorr=corr2(hm1, hmref);
%2-D mean square error
MSEerrHM = immse(hm1,hmref);
%peak signal-to-noise ratio
[peaksnrHM, snrHM] = psnr(hm1, hmref);
%Structural Similarity Index for measuring image quality
[ssimvalHM,ssimmapHM] = ssim(hm1, hmref); 
%image entropy 
hmEntropy = entropy(hm1);   
%entropy ratio of the drawn heatmap in comparison with the reference image
hmEntropy_ratio=hmEntropy/entropy(hmref);
% Compute the mutual information of two images
MI_hm = MI(hm1, hmref);

%Histogram of image data
[pixelCount1,grayLevels1]=imhist(hm1);
[pixelCount2,grayLevels2]=imhist(hmref);

imRatio=mean2(hm1)/mean2(hmref);
isdRatio=std2(hm1)/std2(hmref);

% Histogram-Similarity Distance Measure
%The smaller the value, the more similar the histograms of the two images are....
% So, this value can be used as an indicator of how similar the images are.
% Note that this is a primitive indicator.
% - use chi-squared distance (chisq):
pDchisq=pdist2(pixelCount1', pixelCount2', 'chisq');

% ------------------- removed features [08-05-18]-------------------------%
% mHistIntensity1= sum(pixelCount1 .* grayLevels1) / sum(pixelCount1);
% mHistIntensity2= sum(pixelCount2 .* grayLevels2) / sum(pixelCount2);
% pixelCount1=pixelCount1/size(hm1,1)/size(hm1,2);
% pixelCount2=pixelCount2/size(hm2,1)/size(hm2,2);
% meanBinHeight1=mean(pixelCount1);
% sdBinHeight1 =std(pixelCount1);
% meanBinHeight2=mean(pixelCount2);
% sdBinHeight2=std(pixelCount2);
% overlap = hm1 & hm2; %which discrete pixels match between reference and 
% overlap=~overlap;
% pixels=sum(sum(overlap));
% S=size(hm1); sS=S(1)*S(2);
% pixel_not_matched=double(pixels/sS)*100;
%_________________________________________________________________________%

%-------------------------------------------------------------------------%
%                   Hausdorff Distance Metrics
%-------------------------------------------------------------------------%
%determine the hausdorff distances between interpolated and drawn shape
% - hausD: the hausdorff distance
% - D: the matrix of distances where D(n,m) is the distance of the nth 
%      point in P from the mth point in Q.

[hausD, D, QueryDistances] = HausdorffDist([x1, y1],[x,y],[]);

% % determine the hausdorff distances between interpolated and refernce shape
% [hausDREF, ~] = HausdorffDist([x1,y1],[xref, yref],[]);

nL=size(QueryDistances, 1);
L=length(D);

% % %normalise errors; scaled to number of data points 
% % D=D/L;
% % QueryDistances=QueryDistances/nL;


%Implementation Options:
%%(1) all point distances
%errors=D(:); 
%%(2) all miniumum-point distances
errors=QueryDistances(:); 

%----------------------- Hausdorff Moments ------------------------------%
errorValue=mean(errors);
stdHausD=std(errors);
varHausD=var(errors);
maxHausD=max(errors);
minHausD=min(errors);
skewHausD=skewness(errors);
kurtHausD=kurtosis(errors);
hausD_t=hausD/t(end);
randomError=hmEntropy*hausD;
iqrHausD=iqr(errors);

wp25=round(L*0.25); %   25(%)
wp75=round(L*0.75); %   75(%)
wp15=round(L*0.15); %   15(%)
wp85=round(L*0.85); %   85(%)

HausD_25=sum(errors(1:wp25))/wp25;
HausD_75=sum(errors(wp75:end))/wp25;

% Hausdorff distanes for the middle of the shape (i.e. without
% starting/ending drawing effects
middle_index=wp15:wp85;
[hausDmiddle, ~] = HausdorffDist([x1(middle_index), y1(middle_index)],[x,y],[]);
hausDmiddle_t=hausDmiddle/t(end);

%% Save Features
features=[...
mIm,sdIM,imCorr,peaksnrIM,snrIM,ssimvalIM,...
MSEerrIM,nMSE,nRMSE,MI_im,max_pixel_time,mHM,...
sdHM, hmCorr, MSEerrHM, peaksnrHM, snrHM, ssimvalHM,...
hmEntropy, hmEntropy_ratio, MI_hm, imRatio,...
isdRatio, pDchisq, errorValue, iqrHausD, stdHausD, varHausD,...
maxHausD, minHausD, skewHausD, kurtHausD, hausD_t,...
randomError, HausD_25, HausD_75, hausDmiddle_t];

if length(features)~=length(feature_names)
    error('feature mismatch...'); end 
end
%EOF