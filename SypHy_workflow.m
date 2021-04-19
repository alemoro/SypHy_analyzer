%% SypHy workflow FGA 2.0
% This workflow is intented for analysis of SypHy data after ROI placement.

%% First enter the required parameters
% baselineRange = start and end of the desired baseline in frames
% stimulationRange1 = start and end of the stimulation in frames
% stimulationRange2 = start and end of the stimulation in frames
% ammoniaWashRange = start and end of the NH4Cl superfusion in frames
% acidWashRange = start and end of the pH5.5 superfusion in frames
syphyOption = struct('imagingFrequency', 2,...
    'baselineRange', 1:30,...
    'stimulationRange1', 41:46,...
    'stimulationRange2', 161:166,...
    'ammoniaWashRange', 280:310,...
    'acidWashRange', 320:340,...
    'thresholdSigma', 3.5);

%% Import the data
syphyT = importSypHyData(syphyOption);
%% Filter the active synapses from the inactive onen
syphyT = filterSynapses(syphyT, syphyOption, 1);
%% Calculate the stimulation intensities
syphyT = getStimulationIntensity(syphyT, syphyOption, 'ph5Level', true);
%% Calculate the decay time after stimulation
syphyT = getDecayTime(syphyT, syphyOption);
%% Plot traces
% Plot 'Raw' for the raw data, 'DFF0' for the data normalized to the
% baseline, 'toNH4' for the data normalized to the NH4 response
plotTraces(0:0.5:190, syphyT, syphyOption, 'toNH4');
%% List of funtions
% data importing into table and calculate DF/F0 and normalization to NH4
function syphyT = importSypHyData(syphyOption)
% get the Excel file from the user
[dataFile, dataPath, ~] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)';
    '*.mat', 'Session files (*.mat)'},...
    'Select a SypHy dataset');
globalFile = fullfile(dataPath, dataFile);
% get the all the Excel sheet names
[~, xlFile, ~]        = fileparts(globalFile);
[~, xlAllSheets] = xlsfinfo(globalFile);
% Filter the Excel sheet which contain data, filter based on the sheet name (as 210219_myCondition_cs01_001)
datafilter   = ~cellfun(@isempty,regexpi(xlAllSheets, '^\d{6}_[a-zA-Z_0-9\-\@]*_\w*(\d\[gr])?')); % recognize the feature date_condition_coverslip_cell as in the help
xlDataSheets = xlAllSheets(datafilter);
% calculate the number of sheets containing data and go through the data sheets
numSheets = numel(xlDataSheets);
wait = waitbar(0,'','Name','Loading Data');
syphyT = cell(numSheets+1, 8);
syphyT(1,:) = {'Filename', 'CellID', 'Condition', 'week', 'Time', 'RawData', 'DFF0Data', 'NormalizedToNH4'};
for sheet = 1:numSheets
    waitbar(sheet/numSheets, wait, sprintf('Loading sheet %d/%d',sheet,numSheets))
    sheetname = xlDataSheets{sheet};
    % get the names (date, condition, coverslip, cell)
    [~, ~, data] = xlsread(globalFile, sheetname);
    tempData = cell2mat(data(:,2:end));
    nameparts = regexpi(sheetname, '_', 'split');
    syphyT{sheet+1,1} = xlFile;
    syphyT{sheet+1,2} = sheetname;
    syphyT{sheet+1,3} = nameparts{2};
    syphyT{sheet+1,4} = weeknum(datenum(nameparts{1}, 'yymmdd'));
    syphyT{sheet+1,5} = {(0:size(data,1)) / syphyOption.imagingFrequency};
    syphyT{sheet+1,6} = {tempData};
    % calculate a DF/F0 and a data normalized 0-1 from baseline to NH4 response
    nFrames = size(tempData,1);
    baseline = repmat(mean(tempData(syphyOption.baselineRange,:)),nFrames,1);
    maxNH4 = repmat(max(tempData(syphyOption.ammoniaWashRange,:)),nFrames,1);
    syphyT{sheet+1,7} = {(tempData - baseline) ./ baseline};
    syphyT{sheet+1,8} = {(tempData - baseline) ./ (maxNH4 - baseline)};
end
% show that we're done
close(wait);
syphyT = cell2table(syphyT(2:end, :), 'VariableNames', syphyT(1,:));
syphyT.Condition = categorical(syphyT.Condition);
end

% Select active synapses
function syphyT = filterSynapses(syphyT, syphyOption, startFrom)
answer = questdlg('Manually inspect the data?', 'Filter synapses');
switch answer
    case 'Yes'
        bManual = true;
    case 'No'
        bManual = false;
    otherwise
        return
end
nCell = size(syphyT,1);
allKeep = cell(nCell,1);
allSynPC = zeros(nCell,1);
for cC = startFrom:nCell
    tempData = cell2mat(syphyT{cC, 'DFF0Data'});
    timeData = cell2mat(syphyT{cC, 'Time'});
    timeData = timeData(1:end-1);
    nRoi = size(tempData,2);
    trace = 1;
    bKeep = ones(nRoi,1); % Used to store the real data
    manualKeep = nan(nRoi,1); % Used to store the temporany data
    cmap = [150, 30, 30; 36, 163, 76] / 255;
    while trace <= nRoi
        baseline = mean(tempData(syphyOption.baselineRange, trace));
        baseStDev = std(tempData(syphyOption.baselineRange, trace));
        traceThreshold = baseline + syphyOption.thresholdSigma * baseStDev;
        % Set 1 filter based on the first evoke threshold and a second
        % based on the relative intensity between the stimulation and the
        % NH4 response
        fltr1 = any(tempData(syphyOption.stimulationRange1,trace) > traceThreshold);
        fltr2 = any(tempData(1:syphyOption.ammoniaWashRange(1),trace) >= max(tempData(syphyOption.ammoniaWashRange,trace)));
        if ~fltr1 || fltr2
            bKeep(trace) = 0;
        end
        if ~bManual
            % only used the filters above and go to the next trace
            trace = trace+1;
        else
            % create a figure to display the trace with key pressed functions
            if isnan(manualKeep(trace))
                manualKeep(trace) = bKeep(trace);
            end
            fig = figure;
            plot(timeData, tempData(:,trace), 'Color', cmap(manualKeep(trace)+1,:), 'LineWidth', 2)
            hold on;
            plot(timeData, ones(size(timeData))*traceThreshold, '--k')
            title(['Roi: ', num2str(trace),'/', num2str(nRoi),'. k: Keep - d: Discard - b: Back -n: Next']);
            xlabel('Time (s)'); ylabel('Normalize trace')
            w = waitforbuttonpress;
            if w == 1
                keyPress = fig.CurrentCharacter;
                switch keyPress
                    case 'k'
                        manualKeep(trace) = 1;
                        trace = trace+1;
                    case 'd'
                        manualKeep(trace) = 0;
                        trace = trace+1;
                    case 'b'
                        trace = trace-1;
                    case 'n'
                        trace = trace+1;
                end
            end
            close;
        end
    end
    if ~isnan(manualKeep(1))
        bKeep = manualKeep;
    end
    allKeep{cC} = logical(bKeep');
    allSynPC(cC) = sum(bKeep) / nRoi * 100;
end
syphyT.ActiveSynapses = allKeep;
syphyT.ActiveSynapsesPC = allSynPC;
end

% Function to calculate the intensity during stimulation
function syphyT = getStimulationIntensity(syphyT, syphyOption, varargin)
nStimulation = 1;
if ~isempty(syphyOption.stimulationRange2)
    nStimulation = 2;
end
% nCell = size(syphyT,1);
% allFF0 = zeros(nCell,nStimulation);
% allRelfrac = zeros(nCell,nStimulation);
% allNH4 = zeros(nCell,1);
% for cC = 1:nCell
%     tempFF0 = cell2mat(syphyT{cC, 'DFF0Data'});
%     tempNH4 = cell2mat(syphyT{cC, 'NormalizedToNH4'});
%     activeFltr = cell2mat(syphyT{cC, 'ActiveSynapses'});
%     allNH4(cC,1) = mean(max(tempFF0(syphyOption.ammoniaWashRange, activeFltr)));
%     allFF0(cC,1) = mean(max(tempFF0(syphyOption.stimulationRange1, activeFltr)));
%     allRelfrac(cC,1) = mean(max(tempNH4(syphyOption.stimulationRange1, activeFltr)));
%     if nStimulation == 2
%         allFF0(cC,2) = mean(max(tempFF0(syphyOption.stimulationRange2, activeFltr)));
%         allRelfrac(cC,2) = mean(max(tempNH4(syphyOption.stimulationRange2, activeFltr)));
%     end
% end
syphyT.NH4Intensity = cellfun(@(x,y) mean(mean(x(syphyOption.ammoniaWashRange,y))), syphyT.DFF0Data, syphyT.ActiveSynapses);
syphyT.FF0Stim = cellfun(@(x,y) mean(max(x(syphyOption.stimulationRange1,y))), syphyT.DFF0Data, syphyT.ActiveSynapses);
syphyT.RelFrac = cellfun(@(x,y) mean(max(x(syphyOption.stimulationRange1,y))), syphyT.NormalizedToNH4, syphyT.ActiveSynapses);
if nStimulation == 2
    syphyT.FF0Stim2 = cellfun(@(x,y) mean(max(x(syphyOption.stimulationRange2,y))), syphyT.DFF0Data, syphyT.ActiveSynapses);
    syphyT.RelFrac2 = cellfun(@(x,y) mean(max(x(syphyOption.stimulationRange2,y))), syphyT.NormalizedToNH4, syphyT.ActiveSynapses);
end
if varargin{find(strcmpi(varargin, 'ph5Level'))+1}
    syphyT.pH5 = cellfun(@(x,y) mean(mean(x(syphyOption.acidWashRange,y))), syphyT.DFF0Data, syphyT.ActiveSynapses);
end
end

% Function to calculate the decay constant tau from a single exponential fit of the average trace after the stimulation
function syphyT = getDecayTime(syphyT, syphyOption)
nStimulation = 1;
stimRanges = syphyOption.stimulationRange1;
decayTimes = stimRanges(end)+1:syphyOption.ammoniaWashRange(1)-1;
if ~isempty(syphyOption.stimulationRange2)
    nStimulation = 2;
    stimRanges = [syphyOption.stimulationRange1; syphyOption.stimulationRange2];
    decayTimes = [{stimRanges(1,end)+1:stimRanges(2,1)-1}; {stimRanges(2,end)+1:syphyOption.ammoniaWashRange(1)-1}];
end
nCell = size(syphyT,1);
decayTau = zeros(nCell,nStimulation);
decayData = cell(nCell,nStimulation);
for c = 1:nCell
    tempData = cell2mat(syphyT{c, 'RawData'});
    tempFltr = cell2mat(syphyT{c, 'ActiveSynapses'});
    tempData = tempData(:,tempFltr);
    if ~isempty(tempData)
        % get the raw average for the active synapses
        rawAverage = nanmean(tempData,2);
        % normilize between baseline and end of the stimulation
        nFrames = size(rawAverage,1);
        baseline = repmat(mean(rawAverage(syphyOption.baselineRange,:)),nFrames,1);
        for s = 1:nStimulation
            maxStim = repmat(max(rawAverage(stimRanges(s,:),:)),nFrames,1);
            meanActive = (rawAverage - baseline) ./ (maxStim - baseline);
            % now calculate the fit
            ft = fittype( 'exp1' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Normalize = 'on';
            decayTime = decayTimes{s};
            decayTrace = meanActive(decayTime,end);
            [fitresult, gof] = fit(decayTime', decayTrace, ft, opts);
            decayTau(c,s) = -1*(fitresult.b)^-1;
            decayData{c,s} = meanActive';
        end
    end
end
syphyT.decayData = decayData(:,1);
syphyT.decayTau = decayTau(:,1);
if nStimulation == 2
    syphyT.decayData2 = decayData(:,2);
    syphyT.decayTau2 = decayTau(:,2);
end
end

% Function for plotting the stimulation and the average traces
function plotTraces(time, tempT, tempO, toPlot)
% first adjust the data to plot
switch toPlot
    case 'Raw'
        tempT.plotData = cellfun(@(x,y) mean(x(:,y),2)', tempT.RawData, tempT.ActiveSynapses, 'UniformOutput', false);
        yText = 'SypHy raw intensity (a.u.)';
        hMax = 1e10;
    case 'DFF0'
        tempT.plotData = cellfun(@(x,y) mean(x(:,y),2)', tempT.DFF0Data, tempT.ActiveSynapses, 'UniformOutput', false);
        yText = 'SypHy \DeltaF/F_0 (a.u.)';
        hMax = 10;
    case 'toNH4'
        tempT.plotData = cellfun(@(x,y) mean(x(:,y),2)', tempT.NormalizedToNH4, tempT.ActiveSynapses, 'UniformOutput', false);
        yText = 'SypHy normalized to NH4 (a.u.)';
        hMax = 1.5;
    otherwise
        errordlg("Please use 'Raw', 'DFF0', or 'toNH4'.", 'Data not found');
end
figure; hold on;
if isempty(tempO.stimulationRange2)
    nStimulation = 1;
    stimRanges = tempO.stimulationRange1 / 2;
else
    nStimulation = 2;
    stimRanges = [tempO.stimulationRange1; tempO.stimulationRange2] / 2;
end
for s = 1:nStimulation
    patch([stimRanges(s,1) stimRanges(s,end) stimRanges(s,end) stimRanges(s,1)], [0 0 hMax hMax], [.7 .9 1], 'EdgeColor', 'none', 'FaceAlpha',.3)
end
% color map
cmap = [0.000 0.000 0.000;
    0.058 0.506 0.251;
    0.196 0.600 0.600;
    0.942 0.401 0.250;
    0.700 0.900 1.000];
conds = unique(tempT.Condition);
nCond = numel(conds);
nFrames = size(tempT.plotData{1},2);
time = time(1:nFrames);
fillX = [time fliplr(time)];
for c=1:nCond
    tempFltr = tempT.Condition == conds(c);
    tempData = cell2mat(tempT.plotData(tempFltr,:));
    tempMean = nanmean(tempData);
    tempSEM = nanstd(tempData) ./ sqrt(sum(~isnan(tempData(:,1))));
    fillY = [(tempMean-tempSEM) fliplr(tempMean+tempSEM)];
    fill(fillX,fillY,cmap(c,:),'EdgeColor','none','FaceAlpha',.1)
    hL(c) = plot(time, tempMean, 'Color', cmap(c,:), 'LineWidth', 2);
    legName{c} = sprintf('%s [%d]', char(conds(c)), sum(~isnan(tempData(:,1))));
end
legend(hL, legName); legend boxoff;
set(gca, 'TickDir', 'out');
xlabel('Time (s)');
ylabel(yText);
end