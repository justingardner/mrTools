% mrLoadRet, version 4.0 
%
% AUTHOR:   DJH
% PURPOSE:  Everything.
% DATE:     June, 2004
%
% This is the default mrLoadRet script that simply opens an
% inplaneView window.  You are encouraged to put a copy of this
% file in each of your data directories and modify the copied
% version of this script to customize it for that data set.
% Examples of a number of possible customizations are commented
% below.
function mrLoadRet

if ~isfile('mrSession.mat')
  disp('(mrLoadRet) No mrSession.mat found in current directory');
  return
end

% Define and initialize global variable MLR.
mrGlobals

% Open inplane window
v = mrOpenWindow('Volume');

