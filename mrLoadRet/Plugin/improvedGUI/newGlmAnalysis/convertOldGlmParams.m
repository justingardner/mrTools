% convertOldGlmParams.m
%
%      usage: params = convertOldGlmParams(params)
%         by: julien besle 
%       date: 04/12/2010
%    purpose: re-organizes old GLM params for new GLM code
%              $Id$

function params = convertOldGlmParams(params)

if isfield(params,'testParams')
  params = copyFields(params.testParams,params);
  params = rmfield(params,'testParams');
end

if isfield(params,'contrast') 
  if ~isfield(params,'contrasts')
    params.contrasts = params.contrast;
  end
  params = rmfield(params,'contrast');
end

if isfield(params,'f_tests')
  if ischar(params.f_tests)
    params.fTests = eval(params.f_tests);
  else
    params.fTests = params.f_tests;
  end
  params.numberFtests = size(params.fTests,1);
  params = rmfield(params,'f_tests');
end
if isfield(params,'params')  
  if isfield(params,'fTests')
    for iFtest = 1:size(params.fTests,1)
      params.restrictions{iFtest} = diag(params.fTests(iFtest,:));
    end
    params = rmfield(params,'fTests');
  end
end

if isfield(params,'hrfModel') && strcmp(params.hrfModel,'hrfDiffGamma')
  params.hrfModel='hrfDoubleGamma';
end
if isfield(params,'inplaceConcat') || isfield(params,'applyFiltering')
  params.hrfModel='hrfDeconvolution';
end
  
if isfield(params,'hrfParams')  
  if isfield(params.hrfParams,'incDeriv')
    params.hrfParams.includeDerivative = params.hrfParams.incDeriv;
    params.hrfParams = rmfield(params.hrfParams,'incDeriv');
  end
end

if isfield(params,'trSupersampling')
  params = rmfield(params,'trSupersampling');
end
if isfield(params,'correctionType')
  params.covCorrectionMethod = params.correctionType;
  params = rmfield(params,'correctionType');
end
if isfield(params,'n_rand')
  params.nResamples = params.n_rand;
  params = rmfield(params,'n_rand');
end
if isfield(params,'outputZStatistic')
  params = rmfield(params,'outputZStatistic');
end
if isfield(params,'outputPValue')
  params = rmfield(params,'outputPValue');
end
if isfield(params,'parametricTestOutput')
  if strcmp(params.parametricTestOutput,'T/F value')
    params.outputParametricStatistic =1;
  else
    params.testOutput = params.parametricTestOutput;
  end
  params = rmfield(params,'parametricTestOutput');
end
if isfield(params,'randomizationTestOutput')
  params = rmfield(params,'randomizationTestOutput');
end
if isfield(params,'bootstrapTestOutput')
  params = rmfield(params,'bootstrapTestOutput');
end
