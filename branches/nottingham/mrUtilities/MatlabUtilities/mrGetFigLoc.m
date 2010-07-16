function pos = mrGetFigLoc(figname)
%
% pos = mrGetFigLoc(figname)
%
% Gets a field in the global variable mrDEFAULTS.figloc, which is a
% structure with fields for each figure name.
%
% figname is a string that specifies each type of figure window
% pos is a 4-vector specifying lowerleft corner and size
%
% Examples:
%   pos = mrGetFigLoc('mrLoadRetGUI');
%   pos = mrGetFigLoc('buttondlg');
%   pos = mrGetFigLoc('graphFigure');
%   pos = mrGetFigLoc('mrParamsDialog');
%   pos = mrGetFigLoc('mrParamsDialogHelp');
%
% djh, 5/2007

global mrDEFAULTS

% if mrDEFAULTS is empty, then we should try to load it
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

if ~isempty(mrDEFAULTS) && isfield(mrDEFAULTS.figloc,figname)
    % pos = getfield(mrDEFAULTS.figloc,figname);
    pos = mrDEFAULTS.figloc.(figname);
else
    pos = [];
end
