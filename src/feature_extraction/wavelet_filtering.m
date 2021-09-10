function [features, feature_names] = wavelet_filtering(signal, wname, extra_options)
%Draw-a-Shape Spiral-Specific Feature Exraction
% Function to run the feature extraction for spiral shapes outlined in [1].
%--------------------------------------------------------------------------
% Input:
%       - signal: [1xN] signal vector (e.g. speed);
%       - wname: 'string' corresponding to wavelet name e.g. 'db4'; see
%                  wfilters.m for further details
% _________________________________________________________________________
%    extra_options: structure containing optional inputs to be used in each
%    feature extraction function.
%       - (blank)
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
feature_names={'hf_m', 'hf_cov', 'hf_sd','hf_min', 'hf_max'};
features=NaN(1, length(feature_names));
if nargin < 1
    return
end 
%% Feature Calculation
% Compute the four filters associated with wavelet name given 
% by the input character vector wname. 
%---------------------- plot the filters ---------------------------------%
% [Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname); 
% figure
% subplot(221); stem(Lo_D); 
% title('Decomposition low-pass filter'); 
% subplot(222); stem(Hi_D); 
% title('Decomposition high-pass filter'); 
% subplot(223); stem(Lo_R); 
% title('Reconstruction low-pass filter'); 
% subplot(224); stem(Hi_R); 
% title('Reconstruction high-pass filter'); 
% xlabel(['The four filters for ', wname])

%cD are the approx. coefficients; %cD are the detail coefficients 
[cA,cD] = dwt(signal,wname);

%-------------------------------------------------------------------------%
%                   Moment-based Metrics
%-------------------------------------------------------------------------%
hf_sd=std(cD);
hf_m=mean(cD);
hf_cov=abs(hf_sd)/abs(hf_m);
hf_min = min(cD);
hf_max = max(cD);
%% Save Features
features=[hf_m, hf_cov, hf_sd,hf_min, hf_max];
if length(features)~=length(feature_names)
    error('feature mismatch...'); end 
end
%EOF