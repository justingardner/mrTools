% importTSeries.m
%
%        $Id:$ 
%      usage: importTSeries(v)
%         by: justin gardner
%       date: 07/03/09
%    purpose: 
%
function retval = importTSeries(v)

% check arguments
if ~any(nargin == [1])
  help importTSeries
  return
end

% go find the group that user wants to load here
[filename pathname] = uigetfile({sprintf('*%s',mrGetPref('niftiFileExtension')),'Nifti files'},'Select nifti tSeries that you want to import');

if (filename==0)
  return
end

if ~isempty(strfind(stripext(filename),'.'))
  mrWarnDlg(sprintf('(importTSeries) Ignoring file %s because it has a . in the filename that does not mark the file extension. If you want to use this file, consider renaming to %s',filename,setext(fixBadChars(stripext(filename),{'.','_'}),'hdr')));
  return
end

% get the full file name
fullFilename = fullfile(pathname,filename);

% read the nifti header
hdr = cbiReadNiftiHeader(fullFilename);
if isempty(hdr),return,end
% make sure it is 4 dimensional and then get the number of frames
if hdr.dim(1) ~= 4
  mrWarnDlg(sprintf('(importTSeries) Could not import tSeries because it is a %i dimensional and not 4 dimensional file',hdr.dim(1)));
  return
end
nFrames = hdr.dim(5);

paramsInfo = {{'filename',filename,'The name of the nifti file that you are importing','editable=0'}};
paramsInfo{end+1} = {'description','','A description for the nifti tSeries you are imporint'};
paramsInfo{end+1} = {'nFrames',nFrames,'incdec=[-1 1]',sprintf('minmax=[0 %i]',nFrames),'Number of total frames in your nfiti tSeries'};
paramsInfo{end+1} = {'junkFrames',0,'incdec=[-1 1]',sprintf('minmax=[0 %i]',nFrames),'How many frames should be junked at the beginning'};

params = mrParamsDialog(paramsInfo);
if isempty(params),return,end

scanParams.fileName = params.filename;
scanParams.description = params.description;
scanParams.nFrames = params.nFrames;
scanParams.junkFrames = params.junkFrames;

% now read the file
%tSeries = cbiReadNifti(fullFilename);

v = saveNewTSeries(v,fullFilename,scanParams);



