function h = mrMsgBox(msgstr,modal)
%
% h = mrMsgBox(msgstr)
%
% Calls Matlab's msgbox or disp depending on verbose preference.
% If modal is non-zero then execution is blocked until the user responds.
%
% To set the 'verbose' preference:
%    setpref('mrLoadRet','verbose',0);
%    setpref('mrLoadRet','verbose',1);
%
% djh, 5/2005

if ieNotDefined('modal')
	modal = 0;
end

if ispref('mrLoadRet','verbose')
    verbose = getpref('mrLoadRet','verbose');
else
    verbose = 0;
end

if verbose & modal
	h = msgbox(msgstr,'','modal');
	uiwait(h);
	h =[];
elseif verbose
    h = msgbox(msgstr);
	drawnow;
else
    disp(msgstr);
	h = [];
end
