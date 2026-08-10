function [ver, expectedMatlabVersion,expectedToolboxNames,expectedToolboxIncrements] = mrLoadRetVersion
%     $Id$	

ver = 4.7;

% Change this after testing Matlab upgrades
expectedMatlabVersion = [7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 7.10 7.11 7.12 7.13 7.14 8.0 8.1 8.2 8.4 8.5 9.0 9.1 9.2 9.4 9.5 9.7 9.8 9.13 9.14 9.15];

% expected toolbox
if verLessThan('matlab','8.5')
  expectedToolboxNames = {'Statistics Toolbox','Image Processing Toolbox','Optimization Toolbox'};
else
  expectedToolboxNames = {'Statistics and Machine Learning Toolbox','Image Processing Toolbox','Optimization Toolbox'};
end

expectedToolboxIncrements = {'statistics_toolbox','image_toolbox','optimization_toolbox'};

