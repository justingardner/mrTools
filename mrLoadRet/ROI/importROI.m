% importROI.m
%
%      usage: importROI(view,pathStr)
%         by: justin gardner
%       date: 03/16/07
%    purpose: import roi from mrLoadRet 3.1 version to 4.5 or from Nifti volume
%
function [thisView,params] = importROI(thisView,params,varargin)

% check arguments
if nargin < 1
  help importROI
  return
end

if ieNotDefined('params')
  params = struct;
end
if ischar(params)
  pathStr = params; %for backwards compatibility
  params = struct;
end

if fieldIsNotDefined(params,'from')
  params.from='mrLoadRet';
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end


switch(lower(params.from))
  case 'mrloadret'
    filterSpec={'*.mat','MAT files'; '*.*','All files'};
  case 'nifti'
    filterSpec={'*.nii;*.nii.gz;*.img','NIFTI files'; '*.*','All files'};
end

mrGlobals;
% Complete pathStr
if ieNotDefined('pathStr') && ~justGetParams && ~defaultParams
  % start in an roi directory
  %startPathStr = fullfile(viewGet(view,'viewType'),'ROIs');
  startPathStr = mrGetPref('importROIPath');
  if isempty(startPathStr)
    startPathStr = 'Inplane/ROIs';
  end
  if ~isdir(startPathStr),startPathStr='';,end
  % get the user defined path
  pathStr = mlrGetPathStrDialog(startPathStr,'Choose roi files to import',filterSpec,'on');
elseif justGetParams
  pathStr = 'none';
end
if isempty(pathStr),disp('No ROI selected');,return,end
mrSetPref('importROIPath',fileparts(pathStr{1}));

% get some info
baseNum = viewGet(thisView,'currentBase');
baseXform = viewGet(thisView,'basexform',baseNum);
baseVoxelSize = viewGet(thisView,'baseVoxelSize',baseNum);
baseDims = viewGet(thisView,'baseDims',baseNum);

switch(lower(params.from))
  case 'mrloadret'

    mrMsgBox('(importROI) Note that to import an ROI from the old mrLoadRet, you MUST have the same anatomy image loaded that the ROI was defined on. For example, if it is an inplane ROI, you must have that inplane as the current anatomy. If it is a 3D volume, you will need that 3D volume. Also, make sure that the anatomy has been correctly registered to your base anatomy by mrAlign and has its sform set appropriately. Old mrLoadRet ROIs do not have any alignment information and are just a list of voxels for a particular anatomy image.');


    for roinum = 1:length(pathStr)
      % try to load the roi
      l = load(pathStr{roinum});
      if isfield(l,'ROI')
        clear ROI;
        ROI.name = l.ROI.name;
        ROI.viewType = thisView.viewType;
        ROI.color = l.ROI.color;
        if isfield(l.ROI,'viewType') && ~strcmp(l.ROI.viewType,'Inplane')
          % not sure why gray rois are different from inplane but
          % this seems to work in conversion
          ROI.coords(1,:) = l.ROI.coords(3,:);
          ROI.coords(2,:) = baseDims(2)-l.ROI.coords(2,:)+1;
          ROI.coords(3,:) = baseDims(3)-l.ROI.coords(1,:)+1;
        else
          % there is just an x/y flip for the inplane ROIs
          ROI.coords(1,:) = l.ROI.coords(2,:);
          ROI.coords(2,:) = l.ROI.coords(1,:);
          ROI.coords(3,:) = l.ROI.coords(3,:);
        end
        ROI.coords(4,:) = 1;
        ROI.xform = baseXform;
        ROI.voxelSize = baseVoxelSize;
        ROI.date = datestr(now);
        % Add it to the view
        thisView = viewSet(thisView,'newROI',ROI);
        %ROI.coords
      else
        disp(sprintf('(importROI) No ROI variable found in mat file'));
      end
    end

  case 'nifti'
    
    scanDims = viewGet(thisView,'scanDims');
    scanXform = viewGet(thisView,'scanXform');
    scanVoxelSize = viewGet(thisView,'scanVoxelSize');
    
    for roinum = 1:length(pathStr)  % NB: if scripting, only one ROI can be imported at a time
      if ~justGetParams
        [data,hdr] = mlrImageLoad(pathStr{roinum});
        [~,name] = fileparts(pathStr{roinum});

        % make sure it has only 1 frame
        if hdr.nDim == 3
           hdr.nDim = 4;
           hdr.dim(4) = 1;
        end

        if hdr.dim(4) ~= 1
          mrWarnDlg(sprintf('(importOverlay) Could not import image because it has %d frames',hdr.dim(4)));
          return
        end

        importXformOptions = cell(0);
        if isequal(scanDims,hdr.dim(1:3))
          importXformOptions = putOnTopOfList('scanXform',importXformOptions);
        end
        if isequal(baseDims,hdr.dim(1:3))
          importXformOptions = putOnTopOfList('baseXform',importXformOptions);
        end
        if ~fieldIsNotDefined(hdr,'qform')
          if ~isequal(hdr.sform,eye(4))  % not sure whether sform is always defined
            importXformOptions = putOnTopOfList('niftiQform',importXformOptions);
          end
        end
        if ~fieldIsNotDefined(hdr,'sform')
          if ~isequal(hdr.qform,eye(4)) % not sure whether qform is always defined
            importXformOptions = putOnTopOfList('niftiSform',importXformOptions);
          end
        end
      else
        importXformOptions = {'niftiSform'};
      end
      colors = putOnTopOfList('black',color2RGB);

      paramsInfo = {{'name',name,'The name of the nifti file that you are importing'}};
      paramsInfo{end+1} = {'importXform',importXformOptions,'type=popupmenu','What transformation matrix to use when importing'};
      paramsInfo{end+1} = {'notes','','A description for the ROI you are importing (optional).'};
      paramsInfo{end+1} = {'color',colors,'type=popupmenu','Color of the ROI'};
      
      if defaultParams
        tempParams = mrParamsDefault(paramsInfo);
      else
        tempParams = mrParamsDialog(paramsInfo);
      end
      if isempty(tempParams),return,end
      params = copyFields(tempParams,params);
      params = rmfield(params,{'paramInfo'});
      if justGetParams,return,end
      
      ROI.name = params.name;
      ROI.viewType = 'Volume';
      ROI.color = params.color;
      switch(params.importXform)
        case 'scanXform'
          if ~isequal(scanDims,hdr.dim(1:3))
            mrWarnDlg(sprintf('(importROI) Could not import ROI because its dimensions differ from the current scan dimensions (%s vs %s)', num2str(hdr.dim([1 2 3])'),num2str(scanDims)) );
            return
          else
            ROI.xform = scanXform;
            ROI.voxelSize = scanVoxelSize;
          end
        case 'baseXform'
          if ~isequal(baseDims,hdr.dim(1:3))
            mrWarnDlg(sprintf('(importROI) Could not import ROI because its dimensions differ from the current base dimensions (%s vs %s)', num2str(hdr.dim([1 2 3])'),num2str(baseDims)) );
            return
          else
            ROI.xform = baseXform;
            ROI.voxelSize = baseVoxelSize;
          end
        case 'niftiSform'
          ROI.xform = hdr.sform;
          ROI.voxelSize = hdr.pixdim(1:3);
        case 'niftiQform'
          ROI.xform = hdr.qform;
          ROI.voxelSize = hdr.pixdim(1:3);
      end
      [ROI.coords(:,1),ROI.coords(:,2),ROI.coords(:,3)] = ind2sub(hdr.dim(1:3),find(data));
      ROI.coords = ROI.coords';
      ROI.coords(4,:) = 1;
      ROI.date = datestr(now);
      ROI.notes = params.notes;
      % Add it to the view
      thisView = viewSet(thisView,'newROI',ROI);

    end
end


if exist('ROI','var')
  ROInum = viewGet(thisView,'ROInum',ROI.name);
  if (ROInum > 0)
    thisView = viewSet(thisView,'currentROI',ROInum);
    thisView = viewSet(thisView,'prevROIcoords',[]);
  end
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
