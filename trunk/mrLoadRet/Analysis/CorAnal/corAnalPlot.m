function corAnalPlot(view,overlayNum,scan,x,y,z)
%
% corAnal(view,view,overlayNum,scan,x,y,z)
% 
% djh 9/2005

% If corAnal is loaded, then use it. Otherwise, error.
co = viewGet(view,'co');
amp = viewGet(view,'amp');
ph = viewGet(view,'ph');
if isempty(co) | isempty(amp) | isempty(ph)
	mrErrorDlg('corAnal must be loaded.');
end
coval = co.data{scan}(x,y,z);
ampval = amp.data{scan}(x,y,z);
phval = ph.data{scan}(x,y,z);

% Analysis parameters
ncycles = viewGet(view,'ncycles',scan);
detrend = viewGet(view,'detrend',scan);
spatialnorm = viewGet(view,'spatialnorm',scan);
junkframes = viewGet(view,'junkFrames',scan);
nframes = viewGet(view,'nFrames',scan);
framePeriod = viewGet(view,'framePeriod',scan);
highpassPeriod = round(nframes/ncycles);

% Load tseries from file. Error if file doesn't exist.
pathStr = viewGet(view,'tseriesPathStr',scan);
if ~exist(pathStr,'file')
	mrErrorDlg(['File ',pathStr,' not found']);
end
[tseries,hdr] = cbiReadNifti(pathStr,{x,y,z,[]});
tseries = squeeze(tseries);
tseries = tseries(junkframes+1:junkframes+nframes);

% Remove dc, convert to percent, detrend, and spatial normalization
ptseries = percentTSeries(tseries,...
    'detrend', detrend,...
    'highpassPeriod', highpassPeriod,...
    'spatialNormalization', spatialnorm,...
    'subtractMean', 'Yes',...
	'temporalNormalization', 'No');

% Create model fit
model = ampval * sin(2*pi*ncycles/nframes * [0:nframes-1]' - phval);


% Change xlabel to be in seconds
% Change ylabel contingent upon detrend to include units

% Select window
selectGraphWin;

% Plot it
time = linspace(framePeriod,nframes*framePeriod,nframes)';
plot(time,[ptseries,model],'LineWidth',2);

% Ticks and labels
fontSize = 14;
set(gca,'FontSize',fontSize);
headerStr = ['Times series from voxel [',num2str([x,y,z]),']   co = ',num2str(coval),'  amp = ',num2str(ampval),'  ph = ',num2str(phval)];
title(headerStr);
xtickStep = nframes*framePeriod/ncycles;
xtick = ceil([0:xtickStep:nframes*framePeriod]);
set(gca,'xtick',xtick);
set(gca,'XLim',ceil([0,nframes*framePeriod]));
set(gca,'xgrid','on')
xlabel('Time (sec)');
switch spatialnorm
	case {'Divide by mean'}
		ylabel('fMRI response (% change in image intensity)');
	otherwise
		ylabel('fMRI response(arbitrary units)');
end

return


