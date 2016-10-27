function [residual, covMat, tSeries, modelResponse] = pRFNoise(v,scanNum,coordinates, analysisFileName)

% coordinates = [53,61,25; 53,59,24; 51,57,24; 51,56,22; 52,55,20];
% scanNum = 1;

if ieNotDefined('v')
    v = newView;
end

v = viewSet(v, 'currentGroup', 'Concatenation');
analysisFile = dir(sprintf('Concatenation/pRFAnal/%s', analysisFileName));

% get the analysis structure
analysis = viewGet(v,'analysis');
if ~isfield(analysis,'d') || (length(analysis.d) < scanNum) || isempty(analysis.d)
   v = loadAnalysis(v, ['pRFAnal/' analysisFile.name]);
  analysis = viewGet(v,'analysis');

end

d = viewGet(v,'d',scanNum);
if isempty(d),disp(sprintf('Could not find d structure for this scan'));return,end

numVoxels = size(coordinates, 1);

%figure;
disppercent(-inf, '(pRFNoise) Calculating residuals');
for voxel = 1:numVoxels
    [residual(voxel,:), tSeries(voxel,:), modelResponse(voxel,:)] = getResidual(v,analysis, d, scanNum, coordinates(voxel,1), coordinates(voxel,2), coordinates(voxel,3));
    %subplot(numVoxels,1,voxel)
    %plot(tSeries(voxel,:), 'k.-');
    %hold on;
    %plot(modelResponse(voxel,:), 'r-');
    %plot(residual(voxel,:), 'b-');
    %title(sprintf('voxel %i %i %i', coordinates(voxel,1), coordinates(voxel,2), coordinates(voxel,3)));
    disppercent(voxel/numVoxels);
end
disppercent(inf);
covMat = residual*residual';
    
function [residual, tSeries,modelResponse] = getResidual(v,a,d, scanNum,x,y,z)

% get the params that have been run
scanDims = viewGet(v,'scanDims',scanNum);
whichVoxel = find(d.linearCoords == sub2ind(scanDims,x,y,z));
r = d.r(whichVoxel,:);

params = d.params(:,whichVoxel);
if isfield(d,'paramsInfo')
  paramsInfo = d.paramsInfo;
else
  paramsInfo = [];
end

% get params
m = pRFFit(v,scanNum,x,y,z,'stim',d.stim,'getModelResponse=1','params',params,'concatInfo',d.concatInfo,'fitTypeParams',a.params.pRFFit,'paramsInfo',paramsInfo);

residual = m.tSeries - m.modelResponse;
tSeries = m.tSeries;
modelResponse = m.modelResponse;
