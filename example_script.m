% This script provides an example of how to find outliers within a transect
% of exposure ages based on the stratigraphic relationship using the 
% function 'find_outliers_transect.m'.
%
% Written to be compatible with the iceTEA tools suite by Richard Selwyn 
% Jones, Monash University (richard.s.jones@monash.edu).
%
% The dataset used here is a vertical transect (NEG1) from the North-East 
% Greenland Ice Stream (Roberts et al., in prep).


% SPECIFICATIONS
input_dataset = load('NEG1_ages.mat'); % Load data - file created using iceTEA script Import_Plot_age.m 
mask = [4:14]; % Mask to specify samples to include (here, samples are
    % excluded that are clearly outliers with substantial cosmogenic inheritance)
strat_level = 2; % Assess each age against how many samples either side (e.g. above/below)?
exclude_ends = 1; % Exclude the ages at either end of the transect? - 0=no, 1=yes

% RUN OUTLIER ASSESSMENT
input_dataset.new_ageska = find_outliers_transect(input_dataset.ages_ka,mask,strat_level,exclude_ends);

% PLOT AGES
plot_transect_outliers(input_dataset.ages_ka,input_dataset.new_ageska,'vert',0,mask,[],[]);

