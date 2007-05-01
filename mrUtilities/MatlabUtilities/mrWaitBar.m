function h = mrWaitBar(x,t)
%
% h = mrWaitBar([x],[title or handle])
%
% Calls Matlab's waitbar or disp depending on verbose preference.
% x: percentage of waiting done
% tile or handle: when first opening a waitbar this should be a string
%    which is used as the title of the window. Thereafter, it should be the
%    handle to the waitbar.
%
% To set the 'verbose' preference:
%    setpref('mrLoadRet','verbose',0);
%    setpref('mrLoadRet','verbose',1);
%
% djh, 5/2005

if ~exist('x','var')
    x=0;
end
if ~exist('t','var')
    t = 'Please wait';
end


% update of wait bar
if ishandle(t)
    waitbar(x,t);
    drawnow;
elseif isfield(t,'disppercent')
    disppercent(x);
    % initial call
elseif ischar(t)
    % check the verbose preference
    if ispref('mrLoadRet','verbose')
        verbose = getpref('mrLoadRet','verbose');
    else
        verbose = 0;
    end
    % if verbose, make a window wait bar
    if verbose
        h = waitbar(x,t);
        drawnow;
        % otherwise use disppercent
    else
        disppercent(-inf,t);
        h.disppercent = 1;
    end
end
return


% Test/debug
setpref('mrLoadRet','verbose',1);
setpref('mrLoadRet','verbose',0);

startTime = mglGetSecs;
h = mrWaitBar(0,'Test. Please wait...');
for i=1:100
    mrWaitBar(i/100,h);
end
mrCloseDlg(h);
mglGetSecs(startTime)
