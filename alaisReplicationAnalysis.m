% alaisBurrAnalysis.m
%
%      usage: alaisBurrAnalysis()
%         by: justin gardner
%       date: 10/11/19
%    purpose: analyze data from Alais & Burr replication (experiment alaisburr.m)
%
function e = alaisBurrAnalysis(stimfileNames,varargin)
 
% default return argument
fit = [];
 
% default to working on current directory
if nargin < 1, stimfileNames = [];end
 
% parse arguments
getArgs(varargin,{'dispFit=1','combineData=1','numBootDraws=300','numBootSamples=20'});
 
% get filenames and path
[e.path stimfileNames] = getStimfileNames(stimfileNames);
if isempty(e.path),return,end
 
% check for valid files as we go through
% nFiles will contain how many vaild files we have
e.nFiles = 0;
e.visualStaircase = {};
e.auditoryStaircase = {};
 
% cycle through all files
for iFile = 1:length(stimfileNames)
  % display what is happening
  dispHeader(sprintf('(alaisBurrAnalysis) Analyzing file (%i/%i): %s        ',iFile,length(stimfileNames),stimfileNames{iFile}));
 
  % load and parse the stimfile
  d = loadStimfile(fullfile(e.path,stimfileNames{iFile}));
 
  % valid file, so keep its information
  if ~isempty(d)
    if combineData
      % if combining, then a stimfile that contains the same
      % conditions as a previous one will be combined together
      [isNewFile e] = combineStimfiles(e,d);
    else
      % if not combining then treat every file as new
      isNewFile = 1;
    end
    if isNewFile
      % update count
      e.nFiles = e.nFiles + 1;
      % and list of filenames
      e.filenames{e.nFiles} = stimfileNames{iFile};
      % and the data
      e.d{e.nFiles} = d;
      % see if this is a staircase or psychometric function
      if ~e.d{e.nFiles}.isStaircase
    % keep which ones are psychometric functions
    e.isPsycho(e.nFiles) = 1;
      else
    % collect staircases
    if e.d{e.nFiles}.stimulus.visual
      % add to the visual staircases
      e.visualStaircase{end+1} = e.d{e.nFiles}.stimulus.stair;
    else
      % add to the auditory staircases
      e.auditoryStaircase{end+1} = e.d{e.nFiles}.stimulus.stair;
    end
    % not a psychometric function
    e.isPsycho(e.nFiles) = 0;
      end
    end
  end
end
 
% now cycle through and fit functions
for iFile = 1:e.nFiles
  if ~e.d{iFile}.isStaircase
    % fit psychometric function
    e.d{iFile} = fitPsychometricFunction(e.d{iFile});
    e.d{iFile} = bootstrap(e.d{iFile},numBootDraws,numBootSamples);
  end
end
 
% convert to struct
e.visualStaircase = cell2mat(e.visualStaircase);
e.auditoryStaircase = cell2mat(e.auditoryStaircase);
 
% if no valid files found return
if e.nFiles == 0
  disp(sprintf('(alaisBurrAnalysis) No files found'));
  return
end
 
% calculate/formalize stuff we use for analysis later
[aPSE,aSTD,vPSE,vSTD,bPSE,bSTD,delta,widthArray,deltaArray,visID,audID,bID,aErr,vErr,bErr] = barData(e); %easy-access parameters, summary statistics+their error, and d file specification
[audR2Percent,visR2Percent,bR2Percent] = goodnessOfFit(e,audID,visID,bID) %curve confidence measures of fit
%add curve confidence as d file quality
e.d{audID}.fit.percent = audR2Percent
for i = 1:length(e.d{visID}.visualWidth)
e.d{visID}.fit(i).percent = visR2Percent(i)
end
for i = 1:length(e.d{bID}.condWidth)
e.d{bID}.fit(i).percent = bR2Percent(i)
end
 
% display the fits
if dispFit
  % plot each psychometric function
  for iFile = find(e.isPsycho)
    %diplay fit
    dispFits(e.d{iFile},iFile);
  end
  % display visual staircase
  if ~isempty(e.visualStaircase)
    mlrSmartfig('alaisBurrAnalysis_visualStaircase','reuse');clf;
    doStaircase('threshold',e.visualStaircase,'type=weibull','dispFig=1','titleStr=Visual staircase','useCurrentFig=1');
  end
  % display auditory staircase
  if ~isempty(e.auditoryStaircase)
    mlrSmartfig('alaisBurrAnalysis_auditoryStaircase','reuse');clf;
    doStaircase('threshold',e.auditoryStaircase,'type=weibull','dispFig=1','titleStr=Auditory staircase','useCurrentFig=1');
  end
 
%bar graph of model comparison
for whichWidth = 1:length(vPSE)
    [mle, suboptimal, SCS, SCSdiscrepant, scsSubVisW, scsSubAudW] = calcModelThresholds(whichWidth,vSTD,vPSE,aSTD,aPSE,bPSE,delta)
    % bar graph of different model predictions
    graphModelThresholds(aSTD,vSTD,mle,SCS,SCSdiscrepant,suboptimal,bSTD,widthArray,whichWidth,visID,audID,bID,e,aErr,vErr,bErr)
end
end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    combineStimfiles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isNewFile e] = combineStimfiles(e,d)
 
% default to adding to list
isNewFile = true;
 
% first file
if ~isfield(e,'d'),return,end
 
% only combine bimodal conditions
if d.stimulus.bimodal
  % look for matching stimfiles
  for iStimfile = 1:length(e.d)
    % found a match
    if isequal(e.d{iStimfile}.experimentName,d.experimentName)
      e.d{iStimfile} = concatStimfile(e.d{iStimfile},d);
      isNewFile = false;
      return
    end
  end
  % no matches
end
 
%%%%%%%%%%%%%%%%%%%%%%%%
%    concatStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%
function d = concatStimfile(d1,d2)
 
% set the number of trials
d.nTrials = d1.nTrials + d2.nTrials;
 
% concat fields
d.reactionTime = [d1.reactionTime d2.reactionTime];
d.response = [d1.response d2.response];
 
% copy these fileds
copyFields = {'parameter','randVars'};
for iField = 1:length(copyFields)
  fieldsToConcat = fieldnames(d1.(copyFields{iField}));
  for iConcatField = 1:length(fieldsToConcat)
    d.(copyFields{iField}).(fieldsToConcat{iConcatField}) = [d1.(copyFields{iField}).(fieldsToConcat{iConcatField}) d2.(copyFields{iField}).(fieldsToConcat{iConcatField})];
  end
end
 
% grab fields
grabFields = {'stimulusType','visualWidth','displacement','experimentName','nCond','condNames','isStaircase','condWidth','condDisplacement'};
for iField = 1:length(grabFields)
  d.(grabFields{iField}) = d1.(grabFields{iField});
end
 
% concat the condTrialNums being careful to add the correct number of trials
for iCond = 1:d.nCond
  d.condTrialNums{iCond} = [d1.condTrialNums{iCond} (d2.condTrialNums{iCond}+d1.nTrials)];
end
 
%%%%%%%%%%%%%%%%
% display fits %
%%%%%%%%%%%%%%%%
function dispFits(d,iFit)
 
% open figure
mlrSmartfig(sprintf('%i_%s',iFit,fixBadChars(d.experimentName)),'reuse');clf;
 
nPlots = length(d.visualWidth);
if nPlots > 1
  % plot all together
  nPlots = nPlots+1;
  subplot(1,nPlots,1);
  % display the psychometric functions all together
  dispPsychometricFunction(d,1:d.nCond);
  % now do each width separately
  for iWidth = 1:length(d.visualWidth)
    subplot(1,nPlots,iWidth+1);
    dispPsychometricFunction(d,find(d.condWidth == d.visualWidth(iWidth)));
  end
else
  % display the psychometric functions all together
  dispPsychometricFunction(d,1:d.nCond);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    display a psychometric function    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometricFunction(d,whichConds)
 
% set colors to display in
if (length(whichConds) == 1)
  dataColors = {'k'};
  fitColors = {'r'};
else
  for iCond = 1:length(whichConds)
    dataColors{iCond} = getSmoothColor(iCond,length(whichConds),'cool');
    fitColors{iCond} = getSmoothColor(iCond,length(whichConds),'cool');
  end
end
 
% title string
titleStr = d.stimulusType;
 
for iCond = 1:length(whichConds)
  % plot fit
  plot(d.fit(whichConds(iCond)).fitX,d.fit(whichConds(iCond)).fitY*100,'-','Color',fitColors{iCond});hold on
 
  % plot in percentile
  myerrorbar(d.cond(whichConds(iCond)).uniquePosDiff,100*d.cond(whichConds(iCond)).correctBinned,'yError',100*d.cond(whichConds(iCond)).correctBinnedError,'Symbol','o','MarkerFaceColor',dataColors{iCond});
  xlabel('Probe Offset (visual degrees)');
  yaxis(0,100);
  ylabel('Frequency Identified Rightwards (%)');
  % append fit parameters to title
  if d.stimulusType(1) == 'V'
      titleStr = sprintf('%s\nWidth: %d; Mean: %0.2f; Std: %0.2f; P: %g',titleStr,d.visualWidth(whichConds(iCond)),d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,1-d.fit(whichConds(iCond)).percent);
  end
  if d.stimulusType(1) == 'A'
      titleStr = sprintf('%s\nMean: %0.2f; Std: %0.2f; P: %g',titleStr,d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,1-d.fit(whichConds(iCond)).percent);
  end
  if d.stimulusType(1) == 'B'
      titleStr = sprintf('%s\nWidth: %d; Displacement: %d; Mean: %0.2f; Std: %0.2f; P: %g',titleStr,d.visualWidth(whichConds(ceil(iCond/length(d.displacement)))),d.condDisplacement(whichConds(iCond)),d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,1-d.fit(whichConds(iCond)).percent);
  end
end
 
% display title
title(titleStr);
 
% display a legend for more than one
if length(whichConds) > 1
  % display legend
  hLegend = mylegend({d.condNames{whichConds}},dataColors);
  set(hLegend,'Location','northwest');
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitPsychometricFunction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = fitPsychometricFunction(d)
 
% fit for each set of data in the d structure
for iCond = 1:d.nCond
  % get these trial nums
  trialNums = d.condTrialNums{iCond};
 
  % get the differences
  d.cond(iCond).posDiff = d.parameter.posDiff(trialNums);
  d.cond(iCond).uniquePosDiff = unique(d.cond(iCond).posDiff);
 
  % compute whether answer are correct
  correct = d.parameter.centerWhich(trialNums) == d.response(trialNums);
 
  % bin and average
  for iVal = 1:length(d.cond(iCond).uniquePosDiff)
    % get the trials with the setting of posDiff
    whichTrials = find(d.cond(iCond).posDiff == d.cond(iCond).uniquePosDiff(iVal));
    % compute how many trials
    nTrials = length(whichTrials);
    % compute correct
    d.cond(iCond).correctBinned(iVal) = sum(correct(whichTrials))/nTrials;
    d.cond(iCond).numTrials(iVal) = nTrials;
    % compute ste
    d.cond(iCond).correctBinnedError(iVal) = d.cond(iCond).correctBinned(iVal)*(1-d.cond(iCond).correctBinned(iVal))/sqrt(nTrials);
    % remember nTrials
    d.cond(iCond).nTrials(iVal) = nTrials;
  end
 
  % fit a cumulative gaussian to data
  d.fit(iCond) = fitCumulativeGaussian(d.cond(iCond).uniquePosDiff,d.cond(iCond).correctBinned);
end
 
%%%%%%%%%%%%%%%%%%%%
%    bootstrap     %
%%%%%%%%%%%%%%%%%%%%
function d = bootstrap(d,numBootDraws,numBootSamples)
for iCond = 1:d.nCond
    dev = []; avg = [];
for repeat = 1:numBootSamples
    correct = []; posDiff = [];
    trialNums = d.condTrialNums{iCond};
    d.cond(iCond).uniquePosDiff = unique(d.cond(iCond).posDiff);
    for sample = 1:numBootDraws
        index = round(rand(1)*(length(trialNums)-2)+1);
        correct = [correct d.parameter.centerWhich(d.condTrialNums{iCond}(index)) == d.response(d.condTrialNums{iCond}(index))];
        posDiff = [posDiff d.cond(iCond).posDiff(index)];
    end
   
    for iVal = 1:length(d.cond(iCond).uniquePosDiff);
    % get the trials with the setting of posDiff
    whichTrials = find(posDiff == d.cond(iCond).uniquePosDiff(iVal));
    % compute how many trials
    nTrials = length(whichTrials);
    % compute correct
    d.cond(iCond).correctBinned(iVal) = sum(correct(whichTrials))/nTrials;
    % compute ste
    d.cond(iCond).correctBinnedError(iVal) = d.cond(iCond).correctBinned(iVal)*(1-d.cond(iCond).correctBinned(iVal))/sqrt(nTrials);
    % remember nTrials
    d.cond(iCond).nTrials(iVal) = nTrials;
  end
   
    bootstrap = fitCumulativeGaussian(d.cond(iCond).uniquePosDiff,d.cond(iCond).correctBinned);
    dev = [dev bootstrap.std];
    avg = [avg bootstrap.mean];
end
avg = sort(avg);dev = sort(dev);
d.bootstrap(iCond).meanMax = avg(round(.95*length(avg))); d.bootstrap(iCond).meanMin = avg(round(.05*length(avg)));
d.bootstrap(iCond).stdMax = dev(round(.95*length(dev))); d.bootstrap(iCond).stdMin = dev(round(.05*length(dev)));
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function d = loadStimfile(stimfileName)
 
% default to empty
d = [];
 
% add mat extension
stimfileName = setext(stimfileName,'mat');
 
% see if it exists
if ~isfile(stimfileName)
  disp(sprintf('(alaisBurrAnalysis:loadStimfile) Could not find stimfile: %s',stimfileName));
  return
end
 
% load the file
s = load(stimfileName);
if ~isfield(s,'myscreen') || ~isfield(s,'task')
  disp(sprintf('(alaisBurrAnalysis:loadStimfile) No myscreen or task in stimfile: %s',stimfileName));
  return
end
 
% now check to see if this has the alaisBurr experiment in it
if ~iscell(s.task) || ~iscell(s.task{1})
  disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Task variable has incorrect phases in stimfile: %s',stimfileName));
  return
end
 
% check task filename
taskFilename = s.task{1}{1}.taskFilename;
if isempty(strfind(lower(taskFilename),'alaisburr')) & isempty(strfind(lower(taskFilename),'estimation'))
  disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Incorrect task in stimfile: %s',taskFilename));
  return
end
 
% parse into parameters
d = getTaskParameters(s.myscreen,s.task);
d = d{1};
 
% print what we found
disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,d.nTrials,s.myscreen.SID,s.myscreen.starttime));
 
% get the variables
d.myscreen = s.myscreen;
d.task = s.task;
d.stimulus = s.stimulus;
 
% get experiment type
d.stimulusType = '';
if d.stimulus.visual
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Visual');
end
if d.stimulus.auditory
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Auditory');
end
if d.stimulus.bimodal
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Bimodal');
end
d.stimulusType = strtrim(d.stimulusType);
 
% get width
d.visualWidth = d.stimulus.width;
if isfield(s.task{1}{1}.parameter,'displacement')
  d.displacement = s.task{1}{1}.parameter.displacement;
else
  d.displacement = [];
end
 
% set experiment name and conditions
if d.stimulus.bimodal
  % set experiment name
  d.experimentName = sprintf('%s: [%s] Width: %s',d.stimulusType,mlrnum2str(d.displacement),mlrnum2str(d.visualWidth));
  % set number of conditions
  d.nCond = length(d.displacement) * length(d.visualWidth);
  % get trials for each condition
  iCond = 1;
  for iDisplacementCond = 1:length(d.displacement)
    for iWidthCond = 1:length(d.visualWidth)
      % set condition names
      d.condNames{iCond} = sprintf('Displacement: %s Width: %s',mlrnum2str(d.displacement(iDisplacementCond)),mlrnum2str(d.visualWidth(iWidthCond)));
      % get trials for this condition
      d.condTrialNums{iCond} = find((d.parameter.displacement == d.displacement(iDisplacementCond)) & (d.parameter.width == d.visualWidth(iWidthCond)));
      % remember the parameters
      d.condWidth(iCond) = d.visualWidth(iWidthCond);
      d.condDisplacement(iCond) = d.displacement(iDisplacementCond);
      % update iCond
      iCond = iCond + 1;
    end
  end
elseif d.stimulus.visual
  % set experiment name
  d.experimentName = sprintf('%s: width: %s',d.stimulusType,mlrnum2str(d.visualWidth));
  % set number of conditions
  d.nCond = length(d.visualWidth);
  for iCond = 1:d.nCond
    d.condNames{iCond} = sprintf('Width: %s',mlrnum2str(d.visualWidth(iCond)));
    % set trial nums to all
    d.condTrialNums{iCond} = find(d.parameter.width == d.visualWidth(iCond));
    % remember the parameters
    d.condWidth(iCond) = d.visualWidth(iCond);
  end
else
  % set experiment name
  d.experimentName = sprintf('%s',d.stimulusType);
  % set number of conditions
  d.nCond = 1;
  % set condition name
  d.condNames{1} = sprintf('Auditory');
  % set trial nums to all
  d.condTrialNums{1} = 1:d.nTrials;
end  
 
% check if this is a staircase
if isfield(d.stimulus,'useStaircase') && d.stimulus.useStaircase
  d.isStaircase = 1;
else
  d.isStaircase = 0;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getStimfileNames    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimfilePath stimfileNames] = getStimfileNames(stimfileNames)
 
%  check if we are examining a single mat file
if isfile(setext(stimfileNames,'mat'))
  % make sure the extension is .mat
  stimfileNames = setext(stimfileNames,'mat');
  % first check to see if it has a path
  if ~isempty(fileparts(stimfileNames))
    stimfilePath = fileparts(stimfileNames);
    stimfileNames = getLastDir(stimfileNames);
  else
    stimfilePath = pwd;
  end
else
  % not a single file, so go look for the path
  if isempty(stimfileNames)
    % get current directory
    stimfilePath = pwd;
    % check if it is a directory
  elseif isdir(stimfileNames)
    % then use that as the path
    stimfilePath = stimfileNames;
    % see if it is a subject ID
  elseif ((length(stimfileNames)>1) && (isequal(lower(stimfileNames(1)),'s'))) || isnumeric(stimfileNames)
    % turn a numeric into a string SID
    if isnumeric(stimfileNames)
      stimfileNames = sprintf('s%03i',stimfileNames);
    end
    % get the path for this subject
    stimfilePath = fullfile('~/data/alaisburr',stimfileNames);
  else
    disp(sprintf('(alaisBurrAnalysis) Could not find %s',stimfileNames));
    stimfilePath = '';
    stimfileNames = '';
    return
  end
  % get everything in the directory that is a mat file
  matfiles = dir(fullfile(stimfilePath,'*.mat'));
  stimfileNames = {matfiles(:).name};
end
 
% make sure we are returning a cell array
stimfileNames = cellArray(stimfileNames);
 
 
%%%%%%%%%%%%%%%%%%%%
%      barData     %
%%%%%%%%%%%%%%%%%%%%
function [aPSE,aSTD,vPSE,vSTD,bPSE,bSTD,delta,widthArray,deltaArray,visID,audID,bID,aErr,vErr,bErr] = barData(e)
 
% extract data from every d (type) file
for iFile = 1:e.nFiles
   
    % grab auditory PSE/STD
    if strcmp(e.d{iFile}.stimulusType,'Auditory')
        aPSE = e.d{iFile}.fit.mean;
        aSTD = e.d{iFile}.fit.std;
        aErr = e.d{iFile}.fit.covar;
        aErr.stdLow = e.d{iFile}.bootstrap.stdMin
        aErr.stdHigh = e.d{iFile}.bootstrap.stdMax
        aErr.meanLow = e.d{iFile}.bootstrap.meanMin
        aErr.meanHigh = e.d{iFile}.bootstrap.meanMax
        % save id of auditory file for later
        audID = iFile
    end
   
    % grab visual PSEs/STDs
    if strcmp(e.d{iFile}.stimulusType,'Visual')
        % returns array of values ordered by ascending visual width
        vPSE = [e.d{iFile}.fit.mean];
        vSTD = [e.d{iFile}.fit.std];
        vErr = [e.d{iFile}.fit.covar];
        vErr.stdLow = [e.d{iFile}.bootstrap.stdMin]
        vErr.stdHigh = [e.d{iFile}.bootstrap.stdMax]
        vErr.meanLow = [e.d{iFile}.bootstrap.meanMin]
        vErr.meanHigh = [e.d{iFile}.bootstrap.meanMax]
        % save array of widths for later usage
        widthArray = e.d{iFile}.visualWidth
        % save id of visual file for later
        visID = iFile
    end
   
    % grab bimodal PSEs/STDs
    if strcmp(e.d{iFile}.stimulusType,'Bimodal')
        % returns array ordered by ascending visual width within ascending ordering by delta
        bPSE = [e.d{iFile}.fit.mean];
        bSTD = [e.d{iFile}.fit.std];
        bErr = [e.d{iFile}.fit.covar];
        bErr.stdLow = [e.d{iFile}.bootstrap.stdMin]
        bErr.stdHigh = [e.d{iFile}.bootstrap.stdMax]
        bErr.meanLow = [e.d{iFile}.bootstrap.meanMin]
        bErr.meanHigh = [e.d{iFile}.bootstrap.meanMax]
        % grab delta value and array (assumes only 3 discrepansies)
        delta = max(e.d{iFile}.displacement)
        deltaArray = [e.d{iFile}.displacement]
        % save id of bimodal file for later
        bID = iFile
    end
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     calcModelThresholds    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mle,suboptimal,SCS,SCSdiscrepant,scsSubVisW,scsSubAudW] = calcModelThresholds(whichWidth,vSTD,vPSE,aSTD,aPSE,bPSE,delta)
 
% calculate the MLE threshold
mle = sqrt((vSTD(whichWidth)*vSTD(whichWidth)*aSTD*aSTD)/(vSTD(whichWidth)*vSTD(whichWidth)+aSTD*aSTD))
 
% calculate weights for SCS/suboptimal models using conflict condition
scsSubVisW = (abs(((aPSE-5)-bPSE(whichWidth))/10)+abs(((aPSE+5)-bPSE(whichWidth+2*length(vPSE)))/10))/2
scsSubAudW = 1 - scsSubVisW
 
% suboptimal threshold
suboptimal = sqrt(scsSubVisW*vSTD(whichWidth)^2+(scsSubAudW*aSTD)^2)
 
% SCS threshold
SCS = sqrt((scsSubVisW*(vSTD(whichWidth)*vSTD(whichWidth)))+(scsSubAudW*(aSTD*aSTD))+scsSubVisW*scsSubAudW*(vPSE(whichWidth)-aPSE)^2)
SCSdiscrepant = sqrt((scsSubVisW*(vSTD(whichWidth)*vSTD(whichWidth)))+(scsSubAudW*(aSTD*aSTD))+scsSubVisW*scsSubAudW*(vPSE(whichWidth)-aPSE)^2+scsSubVisW*scsSubAudW*((10)^2))
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   graphModelThresholds   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph thresholds for unimodals, each bimodal conflict condition, and model-predicted thresholds
function graphModelThresholds(aSTD,vSTD,mle,SCS,SCSdiscrepant,suboptimal,bSTD,widthArray,whichWidth,visID,bID,audID,e,aErr,vErr,bErr)
 
% formalize bar graph inputs -- this is all hard coded for 3 deltas!!!!
avgbSTD = bSTD(whichWidth+length(vSTD))
bSTDdiscrep = (bSTD(whichWidth)+bSTD(whichWidth+length(vSTD)+length(vSTD)))/2
graphStats = [vSTD(whichWidth) aSTD avgbSTD bSTDdiscrep mle suboptimal SCS SCSdiscrepant];
errhigh = [vErr.stdHigh(whichWidth) aErr.stdHigh bErr.stdHigh(whichWidth+length(vSTD)) (bErr.stdHigh(whichWidth)+bErr.stdHigh(whichWidth+length(vSTD)+length(vSTD)))/2 graphStats(5) graphStats(6) graphStats(7) graphStats(8)]-graphStats
errlow = [vErr.stdLow(whichWidth) aErr.stdLow bErr.stdLow(whichWidth+length(vSTD)) (bErr.stdLow(whichWidth)+bErr.stdLow(whichWidth+length(vSTD)+length(vSTD)))/2 graphStats(5) graphStats(6) graphStats(7) graphStats(8)]-graphStats
graphWidth = widthArray(whichWidth);
 
% initiate figure
figure(whichWidth+3);
nplots = 4
conditionString = 'Visual Width Condition: %g'
conditionTitle = sprintf(conditionString,graphWidth)
title(conditionTitle)
subplot(1,4,1)
 
% call bar graph
bar([1 2 3 4 5 6 7 8],graphStats); hold on;
er = errorbar([1 2 3 4 5 6 7 8],graphStats,errhigh,errlow);er.LineStyle='none';
hold on
 
% show error bars, taken from covar matrix (bimodal averaged)
%avgBerr = sqrt(((bErr(2,whichWidth*2))^2 + (bErr(2,(2*length(vSTD)+whichWidth*2)))^2 + (bErr(2,(3*length(vSTD)+whichWidth*2)))^2)/3)
%errorbar(graphStats,[vErr(2,(whichWidth*2)) aErr(2,2) avgBerr 0 0 0 0],'.')
 
% calculate p values for model predictions
%[mleH mleP] = ztest(mle,avgbSTD,avgBerr,'tail','both')
%[subH subP] = ztest(suboptimal,avgbSTD,avgBerr,'tail','both')
%[scsH scsP] = ztest(SCS,avgbSTD,avgBerr,'tail','both')
%format shortE
%text(3.55,mle+.2,[num2str(mleP)],'fontsize',7)
%text(4.55,suboptimal+.2,[num2str(subP)],'fontsize',7)
%text(5.55,SCS+.2,[num2str(scsP)],'fontsize',7)
 
% label
set(gca,'XTickLabel',{'Visual','Auditory','Bimodal: Non-discrepant','Bimodal: Discrepant','Optimal Integration','Suboptimal Integration','Stochastic Switching: Non-discrepant','Stochastic Switching: Discrepant'},'XTickLabelRotation',30)

ylabel('Threshold (Visual Degrees)')
xlabel('Observed Data/Model Predictions')
 
intString = 'Threshold Comparison, Width: %g'
graphTitle = sprintf(intString,graphWidth)
title(graphTitle)
hold off
 
% display visual curve
subplot(1,4,2)
dispPsychometricFunction(e.d{visID},[whichWidth])
 
% display auditory curve (indexed as bimodal for some reason)
subplot(1,4,3)
dispPsychometricFunction(e.d{bID},[1])
 
% display bimodal curves (indexed as ??? for some reason)
subplot(1,4,4)
dispPsychometricFunction(e.d{audID},[(whichWidth) (whichWidth+length(vSTD)) (whichWidth+2*length(vSTD))])