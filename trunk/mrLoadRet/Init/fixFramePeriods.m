% mrInit.m
%
%      usage: scanParams = fixFramePeriods(scanParams,weirdFramePeriods,minFramePeriod,maxFramePeriod,tseriesDir)

%         by: julien besle, taken out of mrInit
%       date: 22/07/13
%    purpose: prompts user to fix apparently abnormal frame periods, called by mrInit and saveNewTSeries


function scanParams = fixFramePeriods(scanParams,weirdFramePeriods,minFramePeriod,maxFramePeriod,tseriesDir)

question{1} = 'Abnormal frame periods have been detected';
question{2} = 'Do you want to change the following values ?';
newFramePeriod = weirdFramePeriods;

for iScan = 1:length(weirdFramePeriods)
  if weirdFramePeriods(iScan)
    if weirdFramePeriods(iScan)<minFramePeriod
      newFramePeriod(iScan) = weirdFramePeriods(iScan)*1000;
    elseif weirdFramePeriods(iScan)>maxFramePeriod
      newFramePeriod(iScan) = weirdFramePeriods(iScan)/1000;
    end
    question{end+1} = sprintf('Change %f s into %f s in file %s',weirdFramePeriods(iScan),newFramePeriod(iScan),scanParams(iScan).fileName);
  end
end
yes = askuser(question);
if yes
  for iScan = 1:length(weirdFramePeriods)
    if newFramePeriod(iScan)
      filename = fullfile(tseriesDir, scanParams(iScan).fileName);
      hdr = cbiReadNiftiHeader(filename);
      %set time units to seconds (4th bit = 8)
      niftiSpaceUnit = rem(hdr.xyzt_units, 8); 
      niftiTimeUnit = rem(hdr.xyzt_units-niftiSpaceUnit, 64);
      switch(niftiTimeUnit)
        case 8
          timeUnit = 'sec';
        case 16
          timeUnit = 'msec';
        case 32
          timeUnit = 'microsec';
      end
      hdr.xyzt_units = floor(hdr.xyzt_units/64)+8+niftiSpaceUnit;
      fprintf('Changing frame period from %f %s to %f sec in file %s\n',weirdFramePeriods(iScan),timeUnit,newFramePeriod(iScan),scanParams(iScan).fileName);
      % set the frameperiod
      scanParams(iScan).framePeriod = newFramePeriod(iScan);
      scanParams(iScan).niftiHdr.pixdim(5) = newFramePeriod(iScan); %don't forget to modify the hdr kept in scanParams, otherwise, problems with motionComp (and maybe others)
      %which makes me think that it's not the cleverest thing to keep this niftiHdr field... because someone will always be tempted to change headers outside mrLoadRet anyway
      hdr.pixdim(5) = newFramePeriod(iScan);
      % and write it back
      cbiWriteNiftiHeader(hdr,filename);
    end
  end
end
