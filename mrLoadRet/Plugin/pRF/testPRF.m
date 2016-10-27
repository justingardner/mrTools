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
%       testTSeries:   true voxel response for the left-out scan at a single timepoint
%       modelResp:   expected voxel response trained on other n-1 scans
%       covarMat:   Covariance matrix (size: m x m, where m is # of voxels)
%
%  Output variables:
%       logLikelihood: total likelihood of all the observed tSeries given model response
%       probTable: value at (i, j) is the probability that i'th timepoint in tseries is 
%                  explained by j'th timepoint in model response.

function probTable = testPRF(testTSeries, modelResponse, covarMat)

numTimePoints = size(modelResponse, 2);
probTable = zeros(numTimePoints, numTimePoints);
logprobTable = zeros(numTimePoints, numTimePoints);

disppercent(-inf,'(testPRF) Calculating likelihoods');
for i = 1:numTimePoints
  model_i = modelResponse(:, i);
  tSeries_i = testTSeries(:, i);

  for j = 1:numTimePoints
    tSeries_j = testTSeries(:, j);
    probTable(j, i) = mvnpdf(tSeries_j, model_i, covarMat);
    if probTable(j, i) == 0
      disp(sprintf('\n(testPRF) mvnpdf at index (%i, %i) returning 0; Exiting...', j, i));
      return;
    end
    logProbTable(j, i) = log(probTable(j,i));
  end
  
  %modelLikelihood(i) = mvnpdf(tSeries_i, model_i, covarMat);
  %logLikelihood = logLikelihood + log(modelLikelihood(i));
  disppercent(i/numTimePoints);
end
disppercent(inf);
