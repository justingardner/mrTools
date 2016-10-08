% testPRF.m
%
%     usage: testPRF(    )
%        by: akshay jagadeesh
%      date: 10/07/16
%   purpose: Computes the probability of seeing the model response, 
%            given the observed data, and the computed covariance (noise).
%

function testPRF(voxelList, tSeries, residual, modelResponse, covarMat)

%use mvnpdf to compute the probability density function


%compute log likelihood of the data and sum



