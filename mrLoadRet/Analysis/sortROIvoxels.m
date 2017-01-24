
%sortROIvoxels.m
%
%
% author: Steeve Laquitaine
%   date: 150907
%purpose: Sort or keep voxels based on certain parameters
%
%  usage:
%       
%       %run mrLoadRet    
%       mrLoadRet;
%       v = getMLRView;
%
%       %load roi and pre-computed event-related analysis
%       roi = loadROITSeries(v,'V1');
%       v   = loadAnalysis(v);       
%
%       %keep roi voxels with r2>0.2
%       [newroi,roi] = sortROIvoxels(roi,'r2cutoff=0.2');


function [newroi,roi] = sortROIvoxels(roi,varargin)

%check varargin
if ~isempty(varargin)
    getArgs(varargin);
end

%if r2cutoff
if ~isempty(whos('r2cutoff'))
    
    %get view
    v = getMLRView;
    
    %get scan info
    scanNum = viewGet(v,'curScan');
    
    %check that an r2 map is loaded
    overlay = viewGet(v,'overlayType');
    if ~strcmp(overlay,'r2')
        fprintf('%s \n','(sortROIvoxels) No r2 analysis is loaded. Load one.')
        keyboard
    end
      
    %get voxel r2
    r2all = viewGet(v,'overlayData',scanNum);
    scanDims = viewGet(v,'scanDims');

    %get the linearized index of the matrix "scanCoords"
    %each col contains the 3D coordinates of a voxel (e.g., 10 10 10). We get 
    %the linear index (count elements until coordinate point [10 10 10]) 
    %in a matrix with the scan dimensions. We can then get the r2 values
    %corresponding to those indices in the 3D r2 maps with same scan
    %dimensions.
    newroi = roi;
    linCoor = sub2ind(scanDims,roi.scanCoords(1,:),roi.scanCoords(2,:),roi.scanCoords(3,:));
    r2all = r2all(linCoor);
    
    %update roi with r2cutoff voxels only
    r2cutoffVox         = r2all>r2cutoff;   
    newroi.r2           = r2all(r2cutoffVox);
    newroi.scanCoords   = roi.scanCoords(:,r2cutoffVox);
    newroi.n            = size(newroi.scanCoords,2);   
    newroi.tSeries      = roi.tSeries(r2cutoffVox,:);
    fprintf('\n %s %.2f % s \n','(slfMRIclassifAnalSensoryAreas2) ', r2cutoff, ' has been applied.')
    fprintf('%s %i % s \n','(slfMRIclassifAnalSensoryAreas2)',newroi.n,'voxels survived the r2 cutoff')
end

