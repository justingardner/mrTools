% mrUpdateNiftiHdr.m
%
%      usage: mrUpdateNiftiHdr()
%         by: eli merriam
%       date: 03/20/07
%    purpose: 
%	$Id$	
%
function retval = mrUpdateNiftiHdr()

% check arguments
if ~any(nargin == [0])
  help mrUpdateNiftiHdr
  return
end

view = newView('Volume');
mrGlobals
updateHdr = 0;

% go through this loop twice, first time is just to tel
% user what is going to change
for j = 1:2
  for iGroup = 1:viewGet(view, 'numberofGroups')
    for iScan = 1:viewGet(view, 'nScans', iGroup);
    
      % load the sform44 from the mrSession file
      curhdr = viewGet(view, 'niftiHdr', iScan, iGroup);

      % load the sform44 from the Nifti header
      filename = viewGet(view, 'tseriesPath', iScan, iGroup);
      hdr = cbiReadNiftiHeader(filename);

      % compare them and update if mismatch
      if ~isequal(curhdr,hdr)
	updateHdr = 1;
	if j==1
	  disp(sprintf('Nifti hdr for scan %i in %s group has been modified', iScan, viewGet(view, 'groupName',iGroup)));
	  view = viewSet(view,'niftiHdr',hdr,iScan,iGroup);
	else
	  disp(sprintf('Updating nifti hdr for scan %i in %s group', iScan, viewGet(view, 'groupName',iGroup)));
	  view = viewSet(view,'niftiHdr',hdr,iScan,iGroup);
	end
      end
    end
  end
  % if there is nothing to do, or the user refuses, then return 
  if ~updateHdr || ((j==1) && ~askuser('Ok to change nifti headers in mrSession?'))
    disp('Headers appear up-to-date');
    return
  end
end


% save the mrSession file
saveSession;


