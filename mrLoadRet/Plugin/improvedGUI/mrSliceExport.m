% mrSliceExport - output series of images across all slices
%
%      usage: [  ] = mrSliceExportTest( <v>, <imrange>, <sliceList>, <fname>, <nRows>, <nColumns>, <scanList>) )
%         by: denis schluppeck, modified by julien besle for integration in mrLoadRet
%       date: 2007-11-22
%        $Id: mrSliceExportCor.m 583 2010-10-15 09:54:53Z lpzjb $:
%     inputs:       v: mrLoadRet View structure
%             imrange: relative dimensions of the display to export between 0 and 1 ([minX maxX minY maxY] form left to right and top to bottom)
%           sliceList:
%               fname: base name of the tif file to export to
%               nRows: number of rows in the image
%            scanList:
%    outputs: 
%
%    purpose: exports currently displayed base+overlay as a tif image, (goes across required scans and slices)
%
%        e.g:
%
function [  ]=mrSliceExport( v, imrange, sliceList, fname, nRows, scanList)

mrGlobals

if ieNotDefined('imrange')
  imrange = [0 1 0 1];
end
% clip between 0 and 1
imrange = min(imrange,1);
imrange = max(imrange,0);
  
if ieNotDefined('v')
  v = viewGet([],'view',viewGet([],'viewnums'));
  disp('... using current view')
  %disp(v)
end

if ieNotDefined('fname')
   fname = 'montage00.tif';
end

if ieNotDefined('sliceList')
   sliceList = viewGet(v,'nslices'):-1:1;% 1:nSlices;
end

if ieNotDefined('scanList')
   scanList = viewGet(v,'curscan');
else
   nScans = viewGet(v,'nscans');
   if any(scanList>nScans)
      mrWarnDlg(['There are only ' num2str(nScans) ' scans']);
      return;
   end
end

nSlices = length(sliceList);
if ieNotDefined('nRows')
  nRows = 2;
end
nColumns= ceil(nSlices/nRows);




if isempty(viewGet(v,'curanalysis'))
  % load in current analaysis
  v = loadAnalysis(v, 'corAnal');
end

if isempty(viewGet(v,'curbase'))
  % load anatomy
  v = loadAnat(v);
end

% don't want to hardcode any parameters here, but you could roll your own
% e.g.
% v = viewSet(v, 'curoverlay', 1); % co is 1
% v = viewSet(v,'overlaymin', 0.35);
% v = viewSet(v,'overlaymax', 1.0); 

viewNum = viewGet(v,'viewnum');
curGroup = viewGet(v,'curgroup');
curOverlay = viewGet(v,'currentoverlay');

% for alpha transparency, have to go via mlrGuiSet -- is this bad?
% mlrGuiSet(v,'alpha', 0.75);

% loop over slices and get images

basedims = viewGet(v,'basedims');
newx = ceil(imrange(1:2)*(basedims(1)-1))+1; % indeces are from 1 to basedims
newy =  ceil(imrange(3:4)*(basedims(2)-1))+1;

fprintf(1,'Image coordinates: X: %d -> %d - Y: %d -> %d \n',newx(1),newx(2),newy(1),newy(2));

[pathname, fname, extension] = fileparts(fname);
fname = [fname extension];

% check that the directory exists
if ~exist(pathname, 'dir')
  mkdir(pathname)
  disp(['(mrSliceExport) Created directory ' path]);
else
  % rm all the im*.tif files in Etc
  unix(['rm ' pathname '/im*.tif']);
end

% loop through scans and slices.
for iScan = scanList
   v=viewSet(v, 'curscan', iScan);
   for iSlice = sliceList
      v=viewSet(v, 'curslice', iSlice);

      disp(sprintf('rendering scan %d slice %d .',iScan,iSlice));

      img = refreshMLRDisplay(viewNum);
      reducedim = img(newy(1):newy(2),newx(1):newx(2),:);

      imwrite(reducedim, sprintf('%s/im_%02d_%02d.tif',pathname,iScan, iSlice), 'tif');
   end
end


% and create a montage for a final tif image (this is taken from mrCreateMontage by D. Schluppeck)
dirOutput = dir([pathname '/im*.tif']);
filenames = {dirOutput.name};
for iFile = 1:length(filenames)
  filenames{iFile} = [pathname '/' filenames{iFile}];
end

smartfig('Exported Image');
h_ = montage(filenames, 'size', [nRows nColumns]);
imwrite(get(h_,'cdata'), [pathname '/' fname], 'tiff');
disp(sprintf('wrote: %s', [pathname '/' fname]));

return ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function range = findRange(data)

ampMin = realmax;
ampMax = 0;
nScans = length(data);
for scan=1:nScans
  if ~isempty(data{scan})
    ampMin = min([ampMin min(data{scan}(:))]);
    ampMax = max([ampMax max(data{scan}(:))]);
  end
end
if (ampMin <= ampMax)
  range = [ampMin ampMax];
else
  % if amp data is empty, need to make sure min < max
  range = [0 1];
end




% mrCreateMontage - create montage from tif images in etc folder 
%
%      usage: [  ] = mrCreateMontage(  )
%         by: denis schluppeck
%       date: 2008-09-03
%        $Id$:
%     inputs: 
%    outputs: 
%
%    purpose: get tif files form Etc folder and make a montage
%
%        e.g:
%
function [  ]=mrCreateMontage( fname, nRows, nColumns  )














