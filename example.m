%% Draw-a-Shape Feature Extraction Tutorial
% Tutorial for Draw-a-Shape feature extraction software. This example is 
% based off the work presented by Creagh et al. (2020) [1].
% Reference:
%[1] Creagh, A.P., Simillion, C., Scotland, A., Lipsmeier, F., Bernasconi,
%    C., Belachew, S., van Beek, J., Baker, M., Gossens, C., Lindemann, M. and
%    De Vos, M., 2020. Smartphone-based remote assessment of upper extremity
%    function for multiple sclerosis using the Draw a Shape Test.
%    Physiological measurement, 41(5), p.054002.
%% Andrew Creagh. andrew.creagh@eng.ox.ac.uk
% Last modified in August 2020
%--------------------------------------------------------------------------%
clear
addpath([pwd, '/data/'])
addpath([pwd, '/features/'])
addpath([pwd, '/src/']);addpath([pwd, '/src/feature_extraction/']);addpath([pwd, '/src/functions/'])

%% Load the data
%load example data
load('example_data.mat')

% print the shape data...
shapeData
%arrange into a cell for use in featureExtraction.m
shapeData=table2cell(shapeData);

%Generate dummy patient information
sub_id='PID001'; %participant/subject id
test_id='31-Dec-1999';

%possible shapes:
%{'CIRCLE';'FIGURE_8';'LINE_BOTTOM_TO_TOP';'LINE_TOP_TO_BOTTOM';'SPIRAL';'SQUARE'}
shape='SPIRAL'; 

%generate a information table for this feature file
FILE_INFO=table({sub_id}, {test_id}, {shape}, 'VariableNames', {'sub_id', 'test_id', 'shape'});
%% Run the Feature extraction
%add optional arguments to the extra_options structure; see
%featureExtraction.m for more details
extra_options.sub_id=sub_id;
extra_options.test_id=test_id;

%Perform feature extraction... 
[features, feature_names] = featureExtraction(shapeData, shape, extra_options);

%% Save Features

%Arrange the file informationa and features into a sigle table
FEATURE_TABLE=[FILE_INFO, array2table(features, 'VariableNames', feature_names)];

%generate a feature filename based on the file and shape information
N = matlab.lang.makeValidName([FILE_INFO.sub_id,  FILE_INFO.test_id, {''},{shape}]);
feature_filename=strcat(N{:});

%save the feature file as a delimited .txt file
writetable(FEATURE_TABLE,[pwd, '/features/', 'DaS_feat_' feature_filename, '.txt'],'Delimiter',' ')

%EOF