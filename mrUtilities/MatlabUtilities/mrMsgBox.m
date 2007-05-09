function h = mrMsgBox(msgstr,modal)
%
% h = mrMsgBox(msgstr)
%
% Calls Matlab's msgbox or disp depending on verbose preference.
% If modal is non-zero then execution is blocked until the user responds.
%
% To set the 'verbose' preference:
%    mrSetPref('verbose',0);
%    mrSetPref('verbose',1);
%
% djh, 5/2005
%
% djh, 5/2007, modified to use mrGetPref instead of Matlab's getpref

if ieNotDefined('modal')
	modal = 0;
end

verbose = mrGetPref('verbose');

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
