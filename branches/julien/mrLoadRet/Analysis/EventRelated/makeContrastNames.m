function contrastNames = makeContrastNames(contrasts,EVnames)
% getContrastNames.m
%
%        $Id$
%      usage: contrastNames = getContrastNames(contrasts,EVnames)
%         by: julien besle
%       date: 26/11/2010
%    purpose: makes Contrast Names from contrasts and EV names

contrastNames = cell(1,size(contrasts,1));

for iContrast = 1:size(contrasts,1)
  if nnz(contrasts(iContrast,:))==1 && sum(contrasts(iContrast,:))==1 %if the contrast is one EV
    contrastNames{iContrast} = EVnames{logical(contrasts(iContrast,:))};
  elseif nnz(contrasts(iContrast,:))==2 && sum(contrasts(iContrast,:))==0 %if the contrast is a comparison of 2 EVs
    EV1 = find(contrasts(iContrast,:),1,'first');
    EV2 = find(contrasts(iContrast,:),1,'last');
    if params.computeTtests && strcmp(testParams.tTestSide,'Both')
      connector = ' VS ';
    elseif contrasts(iContrast,EV1)>contrasts(iContrast,EV2)
      connector = ' > ';
    else    
      connector = ' > ';
    end
    contrastNames{iContrast} = [EVnames{EV1} connector EVnames{EV2}];
  elseif all(~diff(contrasts(iContrast,logical(contrasts(iContrast,:))))) %if the contrast is a mean of several EVs
    contrastNames{iContrast} = 'Average(';
    for iEV = find(contrasts(iContrast,:))
      contrastNames{iContrast} = [contrastNames{iContrast} EVnames{iEV} ', '];
    end
    contrastNames{iContrast}(end-1) = ')';
    contrastNames{iContrast}(end) = [];
  else
    contrastNames{iContrast} = mat2str(contrasts(iContrast,:));
  end
end
