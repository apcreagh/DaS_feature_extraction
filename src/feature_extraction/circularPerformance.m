function [features, feature_names]=circularPerformance(xyref, xxyy, xy, extra_options)
%% Draw-a-Shape Circular-Specific Shape Feature Exraction
% Function to run the feature extraction for circular shapes outlined in [1].
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
feature_names={...
'mag_error','mag_error_time','RHO_mean','RHO_sd','RHO_cov','RHO_max',...
'RHO_min','RHO_kurt','RHO_skew','RHO_nPeaks','RHO_nPeaks_norm',...
'RHO_mPeak','RHO_sdPeak','RHO_ps_RHO','RHO_sp_RHO','RHO_apen','RHO_hf_m',...
'RHO_hf_cov','RHO_hf_sd','RHO_hf_min','RHO_hf_max','RV_mean','RV_sd',...
'RV_cov','RV_max','RV_min','RV_kurt','RV_skew','RV_nPeaks',...
'RV_nPeaks_norm','RV_mPeak','RV_sdPeak','RV_ps_RHO','RV_sp_RHO',...
'RV_apen','RV_hf_m','RV_hf_cov','RV_hf_sd','RV_hf_min','RV_hf_max'};

features=NaN(1, length(feature_names));
if nargin < 1
    return
end 

%ApEn Measures
M=2; %window length
R=0.2;%measure of self-simialarity 

%% Extract Shape Data
xref=xyref(:,1);    yref=xyref(:,2);  
x=xxyy(:,1);        y=xxyy(:,2);      t=xxyy(:,3);
x1=xy(:,1);         y1=xy(:,2) ;      t1=xy(:,3);

%% Calculate Features
time=abs(t1(end)-t1(1));

%------------------convert to polar coordinates---------------------------% 
%radial velocity (RV) and angular velocity (RHO)
[AV, timeAV, RV, timeRV, mag_error]=calculatePolarVelocity(x1, y1, t, xref, yref);

%-------------------------------------------------------------------------%
%              Error Feature(s)
%-------------------------------------------------------------------------%
% mag_error;
%time-normalised magnitude error
mag_error_time=mag_error/time; 

%QC: remove NaNs and INFs
[AV, timeAV]=removeNaNs_DS(AV, timeAV);
[RV, timeRV]=removeNaNs_DS(RV, timeRV);

%-------------------------------------------------------------------------%
%              Temporal Features
%-------------------------------------------------------------------------%
%radial velocity...
[RV_features, RV_feature_names]=calcTemporalFeatures(RV, timeRV);
RV_feature_names=strcat('RV_', RV_feature_names);

%angular velocity...
[RHO_features, RHO_feature_names]=calcTemporalFeatures(AV, timeAV);
RHO_feature_names=strcat('RHO_', RHO_feature_names);

%-------------------------------------------------------------------------%
%              Spectral Features
%-------------------------------------------------------------------------%
[RV_wavelet_features, RV_wavelet_feature_names] = wavelet_filtering(RV, 'db10');
RV_wavelet_feature_names=strcat('RV_', RV_wavelet_feature_names);

[RHO_wavelet_features, RHO_wavelet_feature_names] = wavelet_filtering(AV, 'db10');
RHO_wavelet_feature_names=strcat('RHO_', RHO_wavelet_feature_names);

%-------------------------------------------------------------------------%
%              Entropy Features
%-------------------------------------------------------------------------%
% M=2; %window length
% R=0.2;%measure of self-simialarity 
[RV_apen] = approx_entropy(M,R,RV);
[RHO_apen] = approx_entropy(M,R,AV);

RHO_feature_names=[RHO_feature_names, {'RHO_apen'}]; RHO_features=[RHO_features,RHO_apen];
RV_feature_names=[RV_feature_names, {'RV_apen'}]; RV_features=[RV_features, RV_apen];

%% Save Featues
features=[mag_error, mag_error_time, RHO_features,...
RHO_wavelet_features, RV_features, RV_wavelet_features];

feature_names=[{'mag_error', 'mag_error_time'},...
RHO_feature_names,RHO_wavelet_feature_names, RV_feature_names,...
RV_wavelet_feature_names];

if length(features)~=length(feature_names)
    error('feature mismatch...'); end
end 
%EOF