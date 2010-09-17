function analysesList = selectAnalyses(thisView,title)
% analysesList = selectAnalyses(thisView,[title]);
%
%   Gather a list of analyses available in current group/analysis
%   and query the user for a sub-selection.
%
% Output:
%  analysesList: list of selected analyses.
%
% 14/07/10 jb adapted from selectScans
%
% $Id: 

if ieNotDefined('title')
  title = 'Choose analyses';
end

analysisNames = viewGet(thisView,'analysisNames');

%Check for zero:
if isempty(analysisNames)
  mrErrorDlg('No analyses found!');
  return
end

% Which overlays to analyze?
iSel = buttondlg(title, analysisNames);
analysesList = find(iSel);

return;
