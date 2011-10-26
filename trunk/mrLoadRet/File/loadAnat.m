function [view anatFilePath] = loadAnat(view,anatFileName,anatFilePath)
%
% view = loadAnat(view,[anatFileName],[anatFilePath])
%
% Loads an anatomy array and sets view.baseVolumes to include.
%
% anatFileName: Anatomy file. If not specified or empty, prompts user to
% select the file. Anatomy file must be nifti format. anatFileName must be
% in the <homedir>/Anatomy subdirectory. Use symbolic links to point to a
% shared directory for volume anatomy files.
%
% anatFilePath is the path of where to open up the dialog. This
% will bey returned so that the GUI can open up in the same
% place each time.
%
% djh, 1/9/98
% 5/2005, djh, update to mrLoadRet-4.0

%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the anatomy file %
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0,help loadAnat,return,end
  
if ieNotDefined('anatFilePath')
    anatFilePath = ''; 
end

% Open dialog box to have user choose the file
if ieNotDefined('anatFileName')
    if ieNotDefined('anatFilePath')
        startPathStr = viewGet(view,'anatomyDir');
    else
        startPathStr = anatFilePath;
    end
    filterspec = {'*.img;*.nii','Nifti/Analyze files'};
    title = 'Choose anatomy file';
    pathStr = getPathStrDialog(startPathStr,title,filterspec,'on');
elseif ieNotDefined('anatFilePath')
    pathStr = fullfile(viewGet(view,'anatomyDir'),anatFileName);
else
    pathStr = fullfile(anatFilePath,anatFileName);
end

% Aborted
if ieNotDefined('pathStr')
  return
end

% make sure we have a cell array of paths to load
pathStr = cellArray(pathStr);

% cycle through all paths and load each one
for pathNum = 1:length(pathStr)
  % Check whether extension of the file is img or nii
  % matlab function "fileparts" takes the last .foo as extension!
  [path,name,extension] = fileparts(pathStr{pathNum});
  % set default extension, if extension not specified
  if isempty(extension)
    pathStr{pathNum} = setext(pathStr{pathNum},mrGetPref('niftiFileExtension'),0);
 % extension is not .nii or .img
  elseif ~any(strcmp(extension,{'.nii', '.img', '.hdr'}))
    mrWarnDlg(['(loadAnat) File type ',extension,' is not a valid anatomy file format']);
    return
  end

  % File does not exist
  if ~exist(pathStr{pathNum},'file')
    mrWarnDlg(['(loadAnat) File ',pathStr{pathNum},' not found']);
    return
  end

  anatFilePath = path;

  % Load nifti file and reorder the axes and flip (site specific).
  h = mrMsgBox(['(loadAnat) Loading volume: ',pathStr{pathNum},'. Please wait']);
  hdr = cbiReadNiftiHeader(pathStr{pathNum});
  if ishandle(h)
    close(h)
  end

  % get volume dimension...
  % hdr.dim(1) should probably be the number of dimensions
  % but some files don't seem to have this set properly
  % so this check seems to work
  volumeDimension = sum((hdr.dim(2:end)~=1)&(hdr.dim(2:end)~=0));

  % Error if it dimension is greater than 4D.
  if (volumeDimension > 4)
    mrErrorDlg(['(loadAnat) Volume must be 3D or 4D. This file contains a ',num2str(volumeDimension),'D array.']);
  end

  % Handle 4D file
  if (volumeDimension == 4)
    paramsInfo = {{'frameNum',0,'incdec=[-1 1]',sprintf('minmax=[0 %i]',hdr.dim(5)),'This volume is a 4D file, to display it as an anatomy you need to choose a particular time point or take the mean over all time points. Setting this value to 0 will compute the mean, otherwise you can select a particular timepoint to display'}};
    params = mrParamsDialog(paramsInfo,'Choose which frame of 4D file. 0 for mean');
    drawnow
    if isempty(params)
      return
    end
    % if frameNum is set to 0, take the mean
    if params.frameNum == 0
      % either load the whole thing, or if it is too large,
      % then load volume by volume
      if 8*prod(hdr.dim(2:5)) < mrGetPref('maxBlocksize')
	[vol hdr] = cbiReadNifti(pathStr{pathNum});
	vol = nanmean(vol,4);
      else
	% too large, load volume by volume
	vol = zeros(hdr.dim(2:4)');
	disppercent(-inf,sprintf('(loadAnat) Reading %s',getLastDir(pathStr{pathNum})));
	for i = 1:hdr.dim(5)
	  vol = vol + (1/hdr.dim(5))*cbiReadNifti(pathStr{pathNum},{[],[],[],i});
	  disppercent(i/hdr.dim(5));
	end
	disppercent(inf);
      end
      % other wise single time slice
    else
      [vol hdr] = cbiReadNifti(pathStr{pathNum},{[] [] [] params.frameNum});
    end
  else
    % if 3D file just load
    [vol hdr] = cbiReadNifti(pathStr{pathNum});
  end

  % Warning if no alignment information in the header.
  if ~hdr.sform_code
    mrWarnDlg('(loadAnat) No base coordinate frame in the volume header.');
  end

  % Extract permutation matrix to keep track of slice orientation.
  % This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
  % make changes here, please update that function as well.
  [permutationMatrix sliceIndex] = getPermutationMatrix(hdr);

  % keeping this code commented for now
  if 0
    
    % check whether it is left/right reversed, first find which dimension goes from
    % left to right. That is the sliceIndex for a sagittal image (i.e. the sagittal
    % images slice directions goes from left to right).
    lrIndex = sliceIndex(1);
    coronalVector = [0 0 0 1]';
    coronalVector(lrIndex) = 1;
    % now multiply by the qform and see which way we go in the magnet space when we 
    % move one coordinate in this dimension.
    magShift = hdr.qform44*coronalVector-hdr.qform44*[0 0 0 1]';
    if magShift(1) < 0 
      disp(sprintf('(loadAnat) left/right ok'));
    else
      % swap left and right
      disp(sprintf('(loadAnat) left/right reversed, flipping anatomy'));
      vol = flipdim(vol,lrIndex);
      % flip the qform and sform to track this flip
      flipMatrix = eye(4);
      flipMatrix(lrIndex,lrIndex) = -1;
      flipMatrix(lrIndex,4) = size(vol,lrIndex)-1;
      hdr = cbiSetNiftiQform(hdr,hdr.qform44*flipMatrix);
      if hdr.sform_code
	hdr = cbiSetNiftiSform(hdr,hdr.sform44*flipMatrix);
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%
  % Add it to the view %
  %%%%%%%%%%%%%%%%%%%%%%

  % now check for an assoicated .mat file, this is created by
  % mrLoadRet and contains parameters for the base anatomy.
  % (it is the base structure minus the data and hdr) This is
  % essential for flat files

  matFilename = sprintf('%s.mat',stripext(pathStr{pathNum}));
  if isfile(matFilename)
    l = load(matFilename);
    if isfield(l,'base')
      base = l.base;
    end
  end

  % Set required structure fields (additional fields are set to default
  % values when viewSet calls isbase).
  base.name = name;
  base.data = vol;
  base.hdr = hdr;
  base.permutationMatrix = permutationMatrix;

  % Add it to the list of base volumes and select it
  view = viewSet(view,'newBase',base);
end
