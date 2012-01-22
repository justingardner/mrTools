% newGlmAnalysisPlugin.m
%
%        $Id: newGlmAnalysisPlugin.m 1929 2010-12-15 11:55:25Z julien $ 
%      usage: newGlmAnalysisPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/14/10
%
function retval = newGlmAnalysisPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help newGlmAnalysisPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(newGlmAnalysisPlugin) Need a valid view to install plugin'));
  else
    mlrAdjustGUI(thisView,'remove','menu','eventRelatedMenuItem');
    mlrAdjustGUI(thisView,'set','glmMenuItem','Callback',@glmAnalysisCallback);
    mlrAdjustGUI(thisView,'set','recomputeAnalysisMenuItem','Callback',@recomputeAnalysisCallback);

retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Subsumes event-related and "GLM" analyses in a common GLM framework with statistical inference tests (see http://www.psychology.nottingham.ac.uk/staff/ds1/lab/doku.php?id=analysis:glm_statistics_in_mrtools)';
 otherwise
   disp(sprintf('(newGlmAnalysisPlugin) Unknown command %s',action));
end

%------------------------- glmAnalysisCallback Function ------------------------------%
function glmAnalysisCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
glmAnalysis(thisView);

%------------------------- RecomputeAnalysisCallback Function ------------------------------%
% this is to replace the GUI function of old event-related analyses
function recomputeAnalysisCallback(hObject,dump)

viewNum = getfield(guidata(hObject),'viewNum');
thisView = viewGet(viewNum,'view');

if viewGet(thisView,'numAnalyses') > 0
  n = viewGet(thisView,'currentAnalysis');
  groupName = viewGet(thisView,'analysisGroupName',n);
  analysisFunction = viewGet(thisView,'analysisFunction',n);
  guiFunction = viewGet(thisView,'analysisGuiFunction',n);
  switch(guiFunction)
    case {'eventRelatedGUI','eventRelatedGlmGUI'}
      guiFunction = 'glmAnalysisGUI';
      analysisFunction = 'glmAnalysis';
  end
  params = viewGet(thisView,'analysisParams',n);
  % params = guiFunction('groupName',groupName,'dummy',0,'params',params,'thisView',thisView);
  evalstring = ['params = ',guiFunction,'(','''','groupName','''',',groupName,','''','params','''',',params,','''','thisView','''',',thisView);'];
  eval(evalstring);
  % params is empty if GUI cancelled
  if isempty(params)
    disp('(glmAnalysis) Analysis cancelled');
    return
  else
    %check if the group has changed, in which case we need to remove the tseriesFile field so that reconcileParams doesn't get confused
    if isfield(params,'groupName') && ~strcmp(groupName,params.groupName) && isfield(params,'tseriesFile')
       params = rmfield(params,'tseriesFile');
    end
    % thisView = analysisFunction(thisView,params);
    evalstring = ['thisView = ',analysisFunction,'(thisView,params);'];
    eval(evalstring);
    refreshMLRDisplay(viewNum);
  end
else
  mrWarnDlg(sprintf('(mrLoadRetGUI) No analyses loaded'));
end
