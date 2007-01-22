function mlrGuiSet(view,field,value);
%
% mlrGuiSet(view,field,value);
%
% view can be either a view structure or a viewNum. Either way, sets
% handles in the global variable: guidata(MLR.views{viewNum}.figure).
% 
% djh 6/2005

mrGlobals;

if ieNotDefined('view'), mrErrorDlg('No view specified.'); end
if ieNotDefined('field'), mrErrorDlg('No parameter specified'); end
if ieNotDefined('value'), val = []; end

% viewNum can be either a view structure or a view number.
if isnumeric(view)
    viewNum = view;
	% *** Delete the following line after eliminating references to view ***
	view = MLR.views{viewNum};
elseif isview(view)
    viewNum = view.viewNum;
else
    return
end
if (viewNum < 1) | (viewNum > length(MLR.views))
	mrErrorDlg('Invalid viewNum,');
end

% Return if no gui
if isempty(MLR.views{viewNum}.figure)
	return
else
	% Get the gui handles to
	handles = guidata(MLR.views{viewNum}.figure);
end

switch lower(field)

	case {'grouppopup','groupstring'}
		% Set the groupPopup string array
		set(handles.groupPopup,'String',value);

	case {'group'}
		% Choose the group
		set(handles.groupPopup,'Value',value);

	case {'roipopup','roistring'}
		% Set the roiPopup string array
		set(handles.roiPopup,'String',value);

	case {'roi'}
		% Choose the roi
		set(handles.roiPopup,'Value',value);

	case {'basepopup'}
		% Set the basePopup string array
		set(handles.basePopup,'String',value);

	case {'basevolume'}
		% Choose the baseVolume
		set(handles.basePopup,'Value',value);
		
	case {'basedims'}
		% mlrGuiSet(view,'baseDims',[ydim xdim zdim]);
		newDims = value;
		newCoords = min(handles.coords,newDims);
		handles.coords = min(handles.coords,newCoords);
		guidata(MLR.views{viewNum}.figure,handles);

	case {'basemin'}
		% mlrGuiSet(view,'baseMin',value);
		value = clipToSlider(handles.baseMinSlider,value);
		set(handles.baseMinSlider,'Value',value);
		set(handles.baseMinText,'String',num2str(value));
		 
	case {'baseminrange'}
		% mlrGuiSet(view,'baseMinRange',[min,max]);
        if value(2) > value(1)
            set(handles.baseMinSlider,'Min',value(1));
            set(handles.baseMinSlider,'Max',value(2));
        end

	case {'basemax'}
		% mlrGuiSet(view,'baseMax',value);
		value = clipToSlider(handles.baseMaxSlider,value);
		set(handles.baseMaxSlider,'Value',value);
		set(handles.baseMaxText,'String',num2str(value));

	case {'basemaxrange'}
		% mlrGuiSet(view,'baseMinRange',[min,max]);
        if value(2) > value(1)
            set(handles.baseMaxSlider,'Min',value(1));
            set(handles.baseMaxSlider,'Max',value(2));
        end

	case {'analysispopup'}
		% mlrGuiSet(view,'analysisPopup',strings);
		set(handles.analysisPopup,'String',value);

	case {'analysis'}
		% mlrGuiSet(view,'analysis',overlayNum);
		set(handles.analysisPopup,'Value',value);
	
	case {'overlaypopup'}
		% mlrGuiSet(view,'overlayPopup',strings);
		set(handles.overlayPopup,'String',value);

	case {'overlay'}
		% mlrGuiSet(view,'overlay',overlayNum);
		set(handles.overlayPopup,'Value',value);
	
	case {'overlaymin'}
		% mlrGuiSet(view,'overlayMin',value);
		value = clipToSlider(handles.overlayMinSlider,value);
		set(handles.overlayMinSlider,'Value',value);
		set(handles.overlayMinText,'String',num2str(value));

	case {'overlayminrange'}
		% mlrGuiSet(view,'overlayMinRange',[min,max]);
        if (value(2) > value(1))
            set(handles.overlayMinSlider,'Min',value(1));
            set(handles.overlayMinSlider,'Max',value(2));
        end

	case {'overlaymax'}
        % mlrGuiSet(view,'overlayMax',value);
		value = clipToSlider(handles.overlayMaxSlider,value);
		set(handles.overlayMaxSlider,'Value',value);
		set(handles.overlayMaxText,'String',num2str(value));

	case {'overlaymaxrange'}
		% mlrGuiSet(view,'overlayMinRange',[min,max]);
        if (value(2) > value(1))
            set(handles.overlayMaxSlider,'Min',value(1));
            set(handles.overlayMaxSlider,'Max',value(2));
        end

	case {'alpha','overlayalpha'}
		% mlrGuiSet(view,'alpha',value);
		value = clipToSlider(handles.alphaSlider,value);
		set(handles.alphaSlider,'Value',value);
		set(handles.alphaText,'String',num2str(value));

	case {'nscans'}
        % mlrGuiSet(view,'nscans',value);
		nScans = round(value);
		curScan = get(handles.scanSlider,'Value');
		if (nScans > 1)
			set(handles.scanSlider,'Min',1);
			set(handles.scanSlider,'Max',nScans);
			set(handles.scanSlider,'sliderStep',[1/(nScans-1) 2/(nScans-1)]);          
            set(handles.scanSlider,'Visible','on');
            curScan = min(curScan,nScans);
        else
            set(handles.scanSlider,'Min',0.9);
			set(handles.scanSlider,'Max',1.1);
            set(handles.scanSlider,'Visible','off');
            curScan = 1;
        end
        mlrGuiSet(view,'scan',curScan);

	case {'scan'}
        % mlrGuiSet(view,'scan',value);
		value = clipToSlider(handles.scanSlider,value,1);
		set(handles.scanSlider,'Value',value);
		set(handles.scanText,'String',num2str(value));
		% description
		description = viewGet(view,'description',value);
		set(handles.descriptionText,'String',description);

    case {'scantext'}
        % mlrGuiSet(view,'scanText',value);
		set(handles.scanText,'String',num2str(value));
		% description
		description = viewGet(view,'description',value);
		set(handles.descriptionText,'String',description);

	case {'nslices'}
        % mlrGuiSet(view,'nslices',value);
		value = round(value);
        if (value > 1)
            % reset range
            set(handles.sliceSlider,'Min',1);
            set(handles.sliceSlider,'Max',value);
            % reset step
            set(handles.sliceSlider,'sliderStep',[1/(value-1),2/(value-1)]);
            set(handles.sliceSlider,'Visible','on');
        else
            set(handles.sliceSlider,'Visible','off');
            set(handles.sliceText,'String',0);
        end

	case {'slice'}
        % mlrGuiSet(view,'slice',value);
		value = clipToSlider(handles.sliceSlider,value,1);
		handles.coords(handles.sliceOrientation) = value;
		set(handles.sliceSlider,'Value',value);
		set(handles.sliceText,'String',num2str(value));
		guidata(MLR.views{viewNum}.figure,handles);

    case {'slicetext'}
        % mlrGuiSet(view,'sliceText',value);
		handles.coords(handles.sliceOrientation) = value;
		set(handles.sliceText,'String',num2str(value));
		guidata(MLR.views{viewNum}.figure,handles);

	case {'sliceorientation'}
        % mlrGuiSet(view,'sliceorientation',value);
		sliceOrientation = value;
		handles.sliceOrientation = sliceOrientation;
		guidata(MLR.views{viewNum}.figure,handles);
		% Set the correct radio button and unset the other two
		switch sliceOrientation
			% axial
			case 1
				% axial
				set(handles.sagittalRadioButton,'Value',0);
				set(handles.coronalRadioButton,'Value',0);
				set(handles.axialRadioButton,'Value',1);
			case 2
				% coronal
				set(handles.sagittalRadioButton,'Value',0);
				set(handles.coronalRadioButton,'Value',1);
				set(handles.axialRadioButton,'Value',0);
			case 3
				% sagittal
				set(handles.sagittalRadioButton,'Value',1);
				set(handles.coronalRadioButton,'Value',0);
				set(handles.axialRadioButton,'Value',0);
		end

	case {'rotate'}
        % mlrGuiSet(view,'rotate',value);
        value = clipToSlider(handles.rotateSlider,value);
		set(handles.rotateSlider,'Value',value);
		set(handles.rotateText,'String',num2str(value));
		
	case {'viewnum'}
        % mlrGuiSet(view,'viewnum',value);
		handles.viewNum = value;
		guidata(MLR.views{viewNum}.figure,handles);

	otherwise
		error(['Invalid field: ',field]);
end

function value = clipToSlider(slider,value,integerFlag)
% Clips value so that it doesn't exceed slider limits.
% slider is a slider handle
% value must be a number (otherwise, use current slider value)
% integerFlag forces value to be an integer
if ieNotDefined('integerFlag')
	integerFlag = 0;
end
if ~isnumeric(value)
	value = get(slider,'Value');
end
if integerFlag
    if (value < get(slider,'Min'))
        value = ceil(get(slider,'Min'));
    end
    if (value > get(slider,'Max'))
        value = floor(get(slider,'Max'));
    end
else
    if (value < get(slider,'Min'))
        value = get(slider,'Min');
    end
    if (value > get(slider,'Max'))
        value = get(slider,'Max');
    end
end
