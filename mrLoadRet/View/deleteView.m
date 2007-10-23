function deleteView(view)
%
% deleteView(view)
%
% Removes the view from MLR.views
% view can be a view or a viewNum
%
% djh, 6/2004

mrGlobals

if isnumeric(view)
	viewNum = view;
elseif isview(view)
	viewNum = view.viewNum;
else
	return;
end

% Delete it
MLR.views{viewNum} = [];

% check to see if there are any valid views, if not, clear MLR
for i = 1:length(MLR.views)
  if ~isempty(MLR.views{i}),return,end
end
clear global MLR;

% Remove it and update viewNum for all of the remaining views
% (busted because doens't update local variables) 
%
% numviews = length(MLR.views);
% MLR.views = {MLR.views{find(viewNum ~= [1:numviews])}};
% for v = 1:length(MLR.views)
%   MLR.views{v}.viewNum = v;
% 	mlrGuiSet(v,'viewNum',v);
% end

