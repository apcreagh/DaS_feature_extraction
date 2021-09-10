function [features, feature_names]=calcTemporalFeatures(signal, t, extra_options)
%% Draw-a-Shape Temporal Feature Exraction
% Function to run the temporal-based feature extraction for outlined in [1].
%--------------------------------------------------------------------------
% Input:
%       - signal: [1xN] signal vector (e.g. speed);
%       - t:  [1xN] corresponding monitonically increasing time vector
% _________________________________________________________________________
%    extra_options: structure containing optional inputs to be used in each
%    feature extraction function.
%       - extra_options.temporal_filter_window: the median filter window
%         size in samples 
%       - extra_options.correlation_filter_window: the correlation filter 
%         window size in samples
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
feature_names={'mean', 'sd', 'cov','max', 'min', 'kurt', 'skew', 'nPeaks', 'nPeaks_norm', 'mPeak', 'sdPeak', 'ps_RHO', 'sp_RHO'};
features=NaN(1, length(feature_names));
if nargin < 1
    return
end 

temporal_filter_window=5; %set the filter to 5 samples either side
correlation_filter_window=30;

if exist('extra_options', 'var')
    if isfield(extra_options, 'temporal_filter_window')
        temporal_filter_window=options.temporal_filter_window; end 
    if isfield(extra_options, 'correlation_filter_window')
        correlation_filter_window=options.correlation_filter_window; end 
end 
%% Pre-Proccessing
%--------------- Temoporal Specific Quality Control ----------------------%
% Filter the signal with a median filter to smooth over erronous touchpoints 
signal = medfilt1(signal,temporal_filter_window); %median filter 

%discard the first and last 5% of touch points as these are usually very
%error-filled. [07-05-18]: May need a more robust QC aproach in future iterations
discard_pc=0.05;
discard=floor(length(signal)*discard_pc);
n=length(signal)-discard;
signal(n:end)=[];       t(n:end)=[];
signal(1:discard)=[];   t(1:discard)=[];
%% Feature Calculation
%-------------------------------------------------------------------------%
%                   compute moments of the signal
%-------------------------------------------------------------------------%
sig_mean=mean(signal);
sig_sd=std(signal);
sig_kurt=kurtosis(signal);
sig_skew=skewness(signal);
sig_max=max(signal);
sig_min=min(signal);
sig_cov = sig_sd/sig_mean; 

%-------------------------------------------------------------------------%
%                   Signal Amplitude Metrics
%-------------------------------------------------------------------------%
%-------------- normalise the signal for next steps-----------------------%
signal = (signal - min(signal(:))) ./ max(signal(:));

%min peak distance
mpd=length(signal)*0.05; %arbitary 5%

if min(signal) == 0
    mph=std(signal)*1.5;
else
    mph=std(signal)*2; 
end

%get the peaks in the normalised signal
[peaks,idx]=findpeaks(signal,'minpeakheight',  mph, 'minpeakdistance' ,mpd);

%get the time points between peaks (delta-peaks)
deltaPeak=diff(t(idx));
%--------------- compute moments of the delta-peaks ----------------------%
nPeaks=length(deltaPeak);
mPeak=mean(deltaPeak);
sdPeak=std(deltaPeak); 
nPeaks_norm=nPeaks/(length(signal));
%-------------------------------------------------------------------------%
%                   Self-Similarity Metrics
%-------------------------------------------------------------------------%
%moving standard deviation of the signal with a window size of 30
 move_sd=movstd(signal, correlation_filter_window);
%compute the correlation with the variability over time
[ps_RHO,ps_PVAL] = corr(move_sd,t, 'type', 'Pearson'); %Pearson Correlation
[sp_RHO,sp_PVAL]= corr(move_sd,t, 'type', 'Spearman'); %Spearman Correlation

%remove non-significant correlations
if ps_PVAL > 0.05   
    ps_RHO=0; end 
if sp_PVAL > 0.05
    sp_RHO=0; end 

%% Save Featues
features=[sig_mean, sig_sd, sig_cov,sig_max,sig_min,  sig_kurt,...
sig_skew,nPeaks, nPeaks_norm, mPeak, sdPeak,ps_RHO, sp_RHO];
if length(features)~=length(feature_names)
    error('feature mismatch...'); end 
end
%EOF