function [features, feature_names]=calcSpectralFeatures(signal, t, extra_options)
%% Draw-a-Shape Spectral Feature Exraction
% Function to run the spectral (frequency-based) feature extraction for
% outlined in [1].
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
feature_names={'max_pxx','domFreq','aDS_hf_m','aDS_hf_cov','aDS_hf_sd',...
    'aDS_hf_min','aDS_hf_max','speed_apen'};
features=NaN(1, length(feature_names));
if nargin < 1
    return
end 

fs=1/nanmedian(diff(t));

if exist('extra_options', 'var')
    if isfield(extra_options, 'fs')
        fs=extra_options.fs;
    end 
end 
%% SPECTRAL / FREQUENCY FEATURES
%-------------------------------------------------------------------------%
%                          Spectrogram Metrics
%-------------------------------------------------------------------------%
%detrend speed before calculating the spectrogram
%speed=detrend(speed, 'constant'); 
signal=detrend(signal, 'linear'); 
%calculate the periodogram power spectral density (PSD); use hamming window
[pxx, f]=periodogram(signal, hamming(length(signal)), 1024, fs, 'psd');
[max_pxx, mi]=max(pxx);
domFreq=f(mi); %in Hz
% ------------------ plotting spectrogram features -----------------------%
% freq_fig=figure;
% NOTE: smooth(pxx); % for visual asthestic
% plot(f,smooth(pxx), 'k', 'LineWidth', 1.5)
% xlim([0.01, 8])
% xlabel('Frequency (Hz)')
% ylabel('Power (Pixels / s)^2')
% %set(gca,'yticklabels',[])
% set(gca, 'fontsize', 16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% 
% filename=strcat(extra_options.save_pathname, 'FreqPlot_',
% extra_options.group_name,'_',extra_options.subject_name, '_',
% strrep(extra_options.test_date, ' ', '_'), '_',extra_options.shape);
% saveas(freq_fig,strcat(filename, '.png'))
% saveas(freq_fig,strcat(filename, '.fig')) close(freq_fig)
%________________________ end plotting features __________________________%
              
%-------------------------------------------------------------------------%
%                          Wavelet Metrics
%-------------------------------------------------------------------------%
%Calculate wavelet filtering, using a db10 wavelet
aDS=abs(signal);
[aDS_wavelet_features, aDS_wavelet_feature_names] = wavelet_filtering(aDS, 'db10');
aDS_wavelet_feature_names=strcat('aDS_', aDS_wavelet_feature_names);
%-------------------------------------------------------------------------%
%                          Entropy Metrics
%-------------------------------------------------------------------------%
% calculate approximate entropy.... 
M=2; %window length
R=0.2;%measure of self-simialarity 
nsamples=round(fs*5); %5 seconds
%(rescaling added:[04-07-19])
%handy function to perform cublic interpolation to ensure all vectors are
%the same length
A = imresize(signal, [nsamples, 1]);
[speed_apen] = approx_entropy(M,R,A);

%% Save Featues
features=[max_pxx, domFreq, aDS_wavelet_features speed_apen];
feature_names=[{'max_pxx', 'domFreq'}, aDS_wavelet_feature_names, {'speed_apen'}];


if length(features)~=length(feature_names)
    error('feature mismatch...'); end
end 
%EOF