% selectInList.m
%       
%      usage: [numberList,names] = selectInList(thisView,type,title,preselected,groupNum)
%         by: julien ebsle
%       date: 20/12/2010
%        $Id$
%    purpose: asks users to select a list of a given type of mrLoadRet variables from the thisView
%                 type can be: 'overlays','scans','analyses','bases','rois'



function [numberList,names] = selectInList(thisView,type,title,preselected,groupNum)

switch(lower(type))
  case {'scans','scan'}
    type='Scans';
    if ieNotDefined('groupNum')
       groupNum = viewGet(thisView,'currentGroup');
    end
    numberInView = viewGet(thisView,'nScans',groupNum);
    for i = 1:numberInView
      names{i} = sprintf('%i:%s (%s)',i,viewGet(thisView,'description',i,groupNum),viewGet(thisView,'tSeriesFile',i,groupNum));
    end
    
  case {'overlays','overlay'}
    type='Overlays';
    names = viewGet(thisView,'overlayNames');
    
  case {'analyses','analysis'}
    type='Analysis';
    names = viewGet(thisView,'analysisNames');
    
  case {'bases','base'}
    type='Bases';
    names = viewGet(thisView,'baseNames');
    
  case {'rois','roi'}
    type='ROIs';
    names = viewGet(thisView,'roiNames');

end
numberInView = length(names);

%Check for zero:
if numberInView == 0
  mrErrorDlg(['No ' type ' found!']);
  return
end


if ieNotDefined('title')
  title = ['Choose ' type];
end
if ieNotDefined('preselected')
  preselected = [];
end

%add a line for all
if numberInView>1
  names = [names {'All'}];
  numberInView = numberInView+1;
end

preselection = zeros(1,numberInView);
preselection(preselected) = 1;

iSel = buttondlg(title, names,preselection);
if isempty(iSel)
  numberList = iSel; %if cancel has been pressed, this will be a 0*0 matrix, 
  %but if the top close button has been pressed, it will be a 0*1 matrix
else
  if numberInView>1 && iSel(end) %if 'All' has been selected
    numberList = 1:numberInView-1;
  else
    numberList = find(iSel); %if OK is pressed but nothing has been selected, this will output a 1*0 array
  end
end

return;
