% mrFixHdr
%
%      usage: mrFixHdr(<useIdentity>)
%         by: eli merriam
%       date: 03/20/07
%    purpose: Update all nifti headers in groups with
%             the nifti header from the files.
%
%             mrFixHdr;
% 
%             if useIdentity is set, sets all qform and sform
%             matrices to the identity matrix (does not change
%             the actual qform and sform in the nifti hdr, just
%             in the groups, so you can always go back)
%
%             mrUpdateNiftihdr(1);
%
function retval = mrFixHdr(useIdentity)

% check arguments
if ~any(nargin == [0 1])
  help mrFixHdr
  return
end

% default is to update nifti hdrs not set them to identity
if ~exist('useIdentity','var'),useIdentity = 0;,end

% set variables
view = newView('Volume');
mrGlobals
updateHdr = 0;

% go through this loop twice, first time is just to tell
% user what is going to change
for passNum = 1:2
  for iGroup = 1:viewGet(view, 'numberofGroups')
    for iScan = 1:viewGet(view, 'nScans', iGroup);
      % load the nifti header from the mrSession file
      curhdr = viewGet(view, 'niftiHdr', iScan, iGroup);

      % load the acutal Nifti header
      filename = viewGet(view, 'tseriesPath', iScan, iGroup);
      hdr = cbiReadNiftiHeader(filename);
      % get dimensions of voxels from headr
      voxdim = hdr.pixdim(2:4)';
      % get dimensions of voxels from qform44
      xdim = sqrt(sum(hdr.qform44(:,1).^2));
      ydim = sqrt(sum(hdr.qform44(:,2).^2));
      sdim = sqrt(sum(hdr.qform44(:,3).^2));
      voxdimFromQform = [xdim ydim sdim];
      % round to nearest 10000
      voxdim = round(voxdim*1000)/1000;
      voxdimFromQform = round(voxdimFromQform*1000)/1000;
      if ~isequal(voxdim,voxdimFromQform)
	updateHdr = 1;
	% print out info about what is going on
	[os og] = viewGet(view,'originalScanNum',iScan,iGroup);
	if ~isempty(os)
	  disp(sprintf('(%s:%i) pixdim [%s] mismatch with qform44: copy header from (%s:%i) to (%s:%i)',viewGet(view,'groupName',iGroup),iScan,num2str(voxdim),viewGet(view,'groupName',og(1)),os(1),viewGet(view,'groupName',iGroup),iScan));
	else
	  disp(sprintf('(%s:%i) pixdim [%s] mismatch with qform44: No original header to copy from',viewGet(view,'groupName',iGroup),iScan,num2str(voxdim)));
	  % user will have to specify which group/scan to copy from
	  if passNum == 2
	    og = getnum('Which group do you want to copy from? ',1:viewGet(v,'numberOfGroups'));
	    os = getnum('Which scan do you want to copy from? ',1:viewGet(v,'numberOfScans',og));
	  end
	end
	% on second pass, actualy due it
	if passNum == 2
	  % get the source header
	  srchdr = viewGet(view,'niftiHdr',os(1),og(1));
	  % and write the header out
	  cbiWriteNiftiHeader(srchdr,sprintf('%s.hdr',stripext(filename)));
	  % and change it in the view
	  view = viewSet(view,'niftiHdr',iScan,iGroup);
	end
      else
	disp(sprintf('(%s:%i) pixdim [%s] matches with qform44',viewGet(view,'groupName',iGroup),iScan,num2str(voxdim)));
      end
    end
  end

  % if there is nothing to do, or the user refuses, then return 
  if ~updateHdr
    disp('Headers do not need to be fixed');
    return
  elseif ((passNum==1) && ~askuser('Ok to change nifti headers in mrSession?'))
    return
  end
end


% save the mrSession file
saveSession;


