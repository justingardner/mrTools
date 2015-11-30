% pRFRegisterModel.m
%
%        $Id:$ 
%      usage: pRFRegisterModel(modelName,initFit,getStim,initScan,initVoxel,getModel,endScan,endFit)
%         by: justin gardner
%       date: 11/21/15
%    purpose: Register a model to use with pRF Fit code. See pRFModel
%             for an example of how to call. 
%             modelName: String that contains the (arbitrary) name of the model
%             initFit: Function that gets called at the beggining of the fit.
%             getStim:
%             initScan:
%             initVoxel:
%             getModel:
%             endScan:
%             endFit;
%
function retval = pRFRegisterModel(modelName,initFit,getStim,initScan,initVoxel,getModel,endScan,endFit)

% check arguments
if ~any(nargin == [8])
  help pRFRegisterModel
  if nargin >= 1
    mrWarnDlg('(prfRegisterModel) Could not register model %s. Imporper number of input arguments to pRFRegisterModel',modelName);
  end
  return
end

% get the global that contains the models.
global gPRFModels;
modelNum = length(gPRFModels)+1;

% and put the model name and function handles there
gPRFModels(modelNum).modelName = modelName;
gPRFModels(modelNum).initFit = initFit;
gPRFModels(modelNum).getStim = getStim;
gPRFModels(modelNum).initScan = initScan;
gPRFModels(modelNum).initVoxel = initVoxel;
gPRFModels(modelNum).getModel = getModel;
gPRFModels(modelNum).endScan = endScan;
gPRFModels(modelNum).endFit = endFit;

