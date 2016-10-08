% testPRF.m
%
%     usage: testPRF(    )
%        by: akshay jagadeesh
%      date: 10/07/16
%   purpose: Computes the probability of seeing the model response, 
%            given the observed data, and the computed covariance (noise).
%
%
%  Input variables:
%       observe:   true voxel response for the left-out scan at a single timepoint
%       modelResp:   expected voxel response trained on other n-1 scans
%       covarMat:   Covariance matrix (size: m x m, where m is # of voxels)
%       roiVoxels:   list of voxels in this ROI (length m)

function testPRF(voxelList, tSeries, residual, modelResponse, covarMat)

% Y = MVNPDF(X,MU,SIGMA) %returns the density of the multivariate normal
%   distribution with mean MU and covariance SIGMA, evaluated at each row
%   of X.
mvnpdf(testResponse, modelResponse, covarMat)


%compute log likelihood of the data and sum



