function analysesList = selectAnalyses(thisView,title,preselected)
% analysesList = selectAnalyses(thisView,[title],[preselected]);
%
%   Gather a list of analyses available in current group/analysis
%   and query the user for a sub-selection.
%
% Output:
%  analysesList: list of selected analyses.
%
% 14/07/10 jb adapted from selectScans
%
% $Id$ 

analysesList = [];

if ieNotDefined('title')
  title = 'Choose analyses';
end
if ieNotDefined('preselected')
  preselected = [];
end

analysisNames = viewGet(thisView,'analysisNames');

%Check for zero:
if isempty(analysisNames)
  mrWarnDlg('(selectAnalysis) No analysis found!');
  return
end

preselection = zeros(1,length(analysisNames));
preselection(preselected) = 1;

% Which overlays to analyze?
iSel = buttondlg(title, analysisNames,preselection);
if isempty(iSel)
  analysesList = iSel; %if cancel has been pressed, this will be a 0*0 matrix, 
  %but if the top close button has been pressed, it will be a 0*1 matrix
else
  analysesList = find(iSel); %if OK is pressed but nothing has been selected, this will output a 1*0 array
end

return;
