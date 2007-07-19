function xform = computeAlignment(inp, vol, xform, contrastReverse, crop, maxIter);
% xform = computeAlignment()
%
% Oscar Nestares - 5/99

% ON 5/00     - NEW PARAMETER NzLIMIT that controls if slices are replicated or
%               not in the edges. I've put a high value so that almost ALWAYS the
%               slices are replicated before doing the alignment.
%             - Added an alarm when maximum number of iterations is reached, with
%               the possibility of continuing the iterations.
% SL 7/25/02 -  Created regVolInp4 for use in mrAlign4 (based entirely on regVolInp)
% SL 7/31/02 -  Changed instances of A*inv(B), where A and B are both matrices,
%               to be A/B.
% SL 8/01/02 -  Added coarseFlag and fineFlag options
% djh 9/2003 -  Cleaned up various things, fixed a major bug (although I'm
%               entirely sure what it was).
% djh 5/2004 -  Added crop.
% djh 1/2007 -  Update to mrAlign-4.5
% djh 7/2007 -  Added contrastReverse param

inplaneSize = size(inp);
volumeSize = size(vol);
mindisp = 0.1;  % displacement used to end the iterations
CB = 2.5;       % cutoff parameter for robust estimation
SC = 2;         % scaling parameter for robust estimation
nzLIMIT = 16;   % if the number of inplanes is less that this, then
                % it puts a border by replicating the first and last
                % inplanes two times.
permuteM = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];  % For swapping x and y

% If the number of slices is too small, repeat the first and last slice
% to avoid running out of data (the derivative computation discards the
% borders in z, 2 slices at the begining and 2 more at the end)
wbh = mrMsgBox('Please wait: interpolating slices, correcting intensity & contrast...');
if inplaneSize(3) < nzLIMIT
	mrWarnDlg(['Only ',num2str(inplaneSize(3)),' slices. Padding by replicating first and last slice.']);
	inp = cat(3,inp(:,:,1),inp(:,:,1),inp,...
		inp(:,:,end),inp(:,:,end));
    if ~isempty(crop)
        crop(:,3) = crop(:,3)+2;
    end
end
mrCloseDlg(wbh);

% correct intensity & contrast of inplanes 
inpIC = intensityContrastCorrection(inp, crop);
if contrastReverse
    limit = 4; % threshold used by intensityContrastCorrection
    inpIC = inpIC * (pi/limit);
    inpIC = -inpIC;
end

% Loop until the approximate maximum displacement is less than MINDISP, or
% the maximum number of iterations is reached. The maximum displacement is
% calculated aproximately from the sum of terms: the displacement
% corresponding to the rotation of the farthest point in the inplanes, plus
% the norm of the translation.

userResponse = 'Yes';
maxdisp = mindisp + 1;
while strcmp(userResponse, 'Yes')
    niter = 0;
    wbh = mrWaitBar(0,'Computing alignment...');
    while ( ((maxdisp > mindisp) & (niter < maxIter)) | (niter==0))

        % updating number of iterations and message in progress bar
        niter = niter+1;
        mrWaitBar(niter/maxIter,wbh);

        % resample volume according to initial xform
        volInterp = interpVolume(vol, xform, inplaneSize);
        % Pad additional slices, if needed
        if inplaneSize(3) < nzLIMIT
            volInterp = cat(3,volInterp(:,:,1),volInterp(:,:,1),volInterp,...
                volInterp(:,:,end),volInterp(:,:,end));
        end
        % correct intensity & contrast of interpoted volume
        volIC = intensityContrastCorrection(volInterp, crop);
        if contrastReverse
            volIC = volIC * (pi/limit);
        end
        % *** For debugging ***
        % xform
        % figure(1); imagesc(inpIC(:,:,9)); colormap(gray); axis image;
        % figure(2); imagesc(volIC(:,:,9)); colormap(gray); axis image;
        
        % motion estimation (no multiresolution, no iterative)
        M = estMotion3(inpIC, volIC,...
            1,...	      % rotFlag
            1,...	      % robustFlag
            contrastReverse,... % treat image intensities as phase valued
            crop,...      % crop
            CB,...        % cutoff parameter for robust estimation
            SC);          % scale parameter for robust estimation
        % estMotion3 uses [y,x,z] convention so we need to swap x and y
        % before passing coordinates through M and then swap them back
        % before passing them through xform. Yuch. But not worth rewriting
        % estMotion3.
        M = permuteM * M * permuteM;
        xform = xform * M;
        maxdisp = norm(inplaneSize' - M(1:3,1:3)*inplaneSize') + norm(M(1:3,4));
    end
    mrCloseDlg(wbh);
    if (niter == maxIter)
        % ask if we should continue iterating
        quest = strvcat(['Maximum number of iterations (',num2str(maxIter),') reached.'], 'CONTINUE ITERATING?');
        userResponse = questdlg(quest, 'WARNING', 'Yes', 'No', 'No');
    else
        userResponse = 'No';
    end
    
end

