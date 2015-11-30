% pRFRegisterModel.m
%
%        $Id:$ 
%      usage: pRFRegisterModel()
%         by: justin gardner
%       date: 11/21/15
%    purpose: Register a model to use with pRF Fit code.
%
function retval = pRFRegisterModel(modelName,initFit,getStim,initScan,initVoxel,getModel,dispFit)

% check arguments
if ~any(nargin == [7])
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
gPRFModels(modelNum).dispFit = dispFit;

