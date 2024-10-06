clc, close all

% --------------- Member initialisations Start here -----------------------

%thresholdCorr = -1;            % correlation Threshold value
%refroi = 1;                % Cell reference number
sr = 5;                     % sampling rate (Hz)
dt = 1/sr;                  % sample time seconds
traces = Traces;            % file name population of calcium signal   from ROI(numeric matrix)
coord = Coords;             % file name of coordinates center from ROI (numeric matrix)
im = MAX_NewStack;          % file name image from the field of tissue
thresholdCalPerc = 0.2;     % considering 20% of the max normalised values for each cell.
numberOfDummyMatrices = 10;
useRomanoMethod = 1;        % 0 for no 1 for yes
tauDecay = 1.8;             % in seconds

% --------------- Member initialisations end here -----------------------

% Primary step: Get the traces normalised. Normalised traces are requried
% for all further workings.
tracesCombo = getNormalisedTraces(traces, sr, dt, tauDecay, useRomanoMethod);
traces_normalized = tracesCombo{1};
transient = tracesCombo{2};
time = tracesCombo{3};

if useRomanoMethod
    % Alternate way of Romano method.
    [deltaFoF, mu, sigma] = estimateBaselineNoise(traces_normalized);
    movements = zeros(size(traces_normalized, 1));
    [densityData, densityNoise, xev, yev] = NoiseModel(deltaFoF, sigma, movements);
    [mapOfOdds] = SignificantOdds(deltaFoF, sigma, movements, densityData, densityNoise, xev, 0);
    deltaFoF = max(deltaFoF, 0);
    [raster, mapOfOddsJoint] = Rasterize(deltaFoF, sigma, movements, mapOfOdds, xev, yev, 5, tauDecay);
    traces_normalized = deltaFoF;
    % Alternate way of Romano method ends here.
end

figure(1)
plotTraces(traces_normalized, transient, time);

% Draw stacked plot
figure(2)
stackedplot(traces_normalized);

correlOriginalTraces = corrcoef(traces_normalized);
figure(3);
heatmap(correlOriginalTraces);

% Binarise the traces data
%correl_exp = binariseAndReturnCorrCoeff(thresholdCalPerc, traces_normalized, false);
%figure(4)
%heatmap(correl);
 
% Fetch mix-max values for each cell from original sketch matrix
%matrixBoundsForSimul = getCellWiseMixMaxBounds(traces_normalized);

% Create dummy cell traces.
%generateSimulatedTracesBasedOnCellBounds(numberOfDummyMatrices, traces_normalized, matrixBoundsForSimul);

% Get only binarised Matrice from original matrice
if useRomanoMethod
    binarisedTraces = raster;
else
    binarisedTraces = binariseTraces(thresholdCalPerc, traces_normalized, false);
end
rasterplot(find(binarisedTraces), ...
            size(binarisedTraces, 2), size(binarisedTraces, 1));

% Fetch the events in terms on timeline and maintain in a row having arrays for
% each cell representation.
eventTimeRowRecord = getCellTimelineEventsFromBinarisedTrace(binarisedTraces);

% Generate simulated trace matrices based 
%simulatedTraceMatrices = generateSimulatedTracesFromEventRecord(numberOfDummyMatrices, ...
%    binarisedTraces, eventTimeRowRecord);

simulatedData = generateSimulatedTracesFromEventRecord(numberOfDummyMatrices, ...
    binarisedTraces, eventTimeRowRecord);

simulatedMatrices = simulatedData{1};

% Plot simulated matrices
% for i = 1:numel(simulatedMatrices)
%     matrice = simulatedMatrices{i};
%     figure(i);
%     rasterplot(find(matrice), ...
%             size(matrice, 2), size(matrice, 1));
% end

% Generate percentage aggreement matrice on original traces binarised data.
originalPerAggrMatrice = generatePercentagAggrMatrice(binarisedTraces);

% Generate percentage aggreement matrice on simulated traces binarised data.
perAggrMatriceList = cell(numel(simulatedMatrices), 1);
for i = 1:numel(simulatedMatrices)
    simulatedMatrice = simulatedMatrices{i};
    perAggrMatrice = generatePercentagAggrMatrice(simulatedMatrice);
    perAggrMatriceList{i} = perAggrMatrice;
end

% testSignificanceList = cell(numel(perAggrMatriceList), 1);
% combinedTestSigMatrice = zeros(size(binarisedTraces, 2));
% for i = 1:numel(perAggrMatriceList)
%     simulatedPerAggrSampleMatrice = perAggrMatriceList{i};
%     meanTmp = mean(mean(simulatedPerAggrSampleMatrice));
%     stdTmp = mean(std(simulatedPerAggrSampleMatrice));
%     pValue = meanTmp + 2 * stdTmp;
%     testSignificance = originalPerAggrMatrice > pValue;
%     testSignificanceList{i} = double(testSignificance);
%     combinedTestSigMatrice = combinedTestSigMatrice + testSignificance;
% end
% 
% figure(5);
% heatmap(combinedTestSigMatrice);

% histogram(originalPerAggrMatrice);
% grid on;
% %xlabel('Data Value', 'FontSize', 15);
% %ylabel('Count', 'FontSize', 15);
% % Compute mean and standard deviation.
% mu = mean(mean(originalPerAggrMatrice));
% sigma = mean(std(originalPerAggrMatrice));
% meanTwo = mu + 2 * sigma;
% % Indicate those on the plot.
% xline(mu, 'Color', 'g', 'LineWidth', 2);
% xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
% xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
% %xline(meanTwo - sigma, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
% xline(meanTwo, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');

% Prepare pairwise matrice population of all simulated matrices.
combinedPASMatrice = [];
for i = 1:numel(perAggrMatriceList)
    simulatedPerAggrSampleMatrice = perAggrMatriceList{i};
    combinedPASMatrice = [combinedPASMatrice, simulatedPerAggrSampleMatrice];
end

titleText = 'Pairwise % Aggreement';

% Plot pairwise % aggreement
figure(6);
combinedAggrThreshold = plotHistogramWithStdAndMean(combinedPASMatrice, titleText);

figure(7);
plotHeatmapWithSignificance(originalPerAggrMatrice, combinedAggrThreshold, titleText);

meansOfSimulatedMatrices = [];
for i = 1:numel(perAggrMatriceList)
    simulatedPerAggrSampleMatrice = perAggrMatriceList{i};
    meanTmp = mean(mean(simulatedPerAggrSampleMatrice));
    meansOfSimulatedMatrices = [meansOfSimulatedMatrices, meanTmp];
end

% Plot population % aggreement
figure(8);
originalPAMatMean = mean(mean(originalPerAggrMatrice));
aggrTshldOfMeans = plotHistogramWithStdAndMean(meansOfSimulatedMatrices, 'Population % Aggreement', originalPAMatMean);

%figure(9);
%plotHeatmapWithSignificance(originalPerAggrMatrice, aggrTshldOfMeans, titleText);

observedSignificanceMatrix = originalPerAggrMatrice > combinedAggrThreshold;
figure(10);
matTitle = sprintf('Observed Signicance Matrix, Combined Threshold = %.3f', combinedAggrThreshold);
htMap = heatmap(double(observedSignificanceMatrix), 'CellLabelColor','none');
htMap.Title = matTitle;

trueValObservedSignificanceMatrix = double(observedSignificanceMatrix);

%figure(11);
%gplot(observedSignificanceMatrix, coord, '--ob');
%networkvisualizer(observedSignificanceMatrix, )

%prepareRandomNetwork = @(n1, n2, numedges) logical(sparse(...
 %  randi([1 n1], numedges, 1), randi([1 n2], numedges, 1), 1, n1, n2));
%W = prepareRandomNetwork(50, 50, 100);
%net = networkvisualizer(observedSignificanceMatrix, size(Coords,1));

% ------------------------------ List out all local functions below this section ------------------------------|

function perAggrMatrice = generatePercentagAggrMatrice(simulatedMatrice)
    perAggrMatrice = zeros(size(simulatedMatrice, 2));
    for cellIndex = 1: size(simulatedMatrice, 2)
        for compCellIndex = 1: size(simulatedMatrice, 2)
            activeMatchCounter = 0;
            totalActiveCounter = 0;
            if (cellIndex == compCellIndex)
                perAggrMatrice(cellIndex, compCellIndex) = 0;
            else
                for cellRowIndex = 1: size(simulatedMatrice, 1) 
                    if (simulatedMatrice(cellRowIndex, cellIndex) == 1 && ...
                            simulatedMatrice(cellRowIndex, compCellIndex) == 1)
                        activeMatchCounter = activeMatchCounter + 1;
                    end
                    if (simulatedMatrice(cellRowIndex, cellIndex) == 1 || ...
                            simulatedMatrice(cellRowIndex, compCellIndex) == 1)
                        totalActiveCounter = totalActiveCounter + 1;
                    end
                end
                if (activeMatchCounter > 0 && totalActiveCounter > 0) 
                    perAggrOfComparedCells = activeMatchCounter / totalActiveCounter;
                    perAggrMatrice(cellIndex, compCellIndex) = perAggrOfComparedCells;
                else
                    perAggrMatrice(cellIndex, compCellIndex) = 0;
                end
            end
        end
    end
end

function plotTraces(traces_normalized, transient, time)
    k = 0;
    hold on
    set(gcf, 'color', 'w', 'units', 'centimeters', ...
        'position', [1 1 30 18], ...
        'defaultaxesfontname', 'arial', ...
        'defaultaxesfontsize', 16)
    for i = 1:size(traces_normalized, 2)
        a = traces_normalized(:, i);
        evt = a(transient{1, i});
        plot(time, a)
        plot(time(1, transient{1, i}), evt, 'ro', 'MarkerSize', 5, 'LineWidth', 1);
    end
    axis off
    
    plot([-10 20], [-0.1 -0.1], 'k', 'linewidth', 3)
    plot([-10 -10], [-0.1 0.4], 'k', 'linewidth', 3)
    xtext = [-50 0.01];
    ytext = [0 -0.08];
    str1 = {'50% df/f', '30 sec'};
    text(xtext, ytext , str1)
    %add mean trace
    traces_normalized_tab = traces_normalized';
    M = mean(traces_normalized_tab);
    plot(time, M, 'color', 'black', 'LineWidth', 3)
    hold off
end

% Function created a binary representation of the normalised trace
% matrix(traces_normalized) and returns a correlated coefficient matrix.
%
% Params: (thresholdCalPerc, traces_normalized, plotRaster)
%          plotRaster - Boolean value, pass 'true' to plot raster from binary matrix.
%
function correl = binariseAndReturnCorrCoeff(thresholdCalPerc, traces_normalized, plotRaster)
    %Step 1: Create a row matrix with highest intensities for each cell value
    %        considering 20% of the max normalised values for each cell.
    thresholdValuesForBinarise = thresholdCalPerc * max(traces_normalized);  
    
    %Step 2: Binarise the normalised traces based on the threshold determined
    %for each cell in matrix [thresholdValuesForBinarise].
    binarisedTraces = traces_normalized > thresholdValuesForBinarise;
    binarisedTracesMatrix = double(binarisedTraces);
    if plotRaster
        rasterplot(find(binarisedTracesMatrix), ...
            size(binarisedTracesMatrix, 2), size(binarisedTracesMatrix, 1));
    end
    
    % Do the correlation coefficient on the binarised matrix
    correl = corrcoef(binarisedTracesMatrix);
end

% Function created a binary representation of the normalised trace
% matrix(traces_normalized) and returns matrix.
%
% Params: (thresholdCalPerc, traces_normalized, plotRaster)
%          plotRaster - Boolean value, pass 'true' to plot raster from binary matrix.
%
function binarisedTracesMatrix = binariseTraces(thresholdCalPerc, traces_normalized, plotRaster)
    %Step 1: Create a row matrix with highest intensities for each cell value
    %        considering 20% of the max normalised values for each cell.
    %thresholdValuesForBinarise = thresholdCalPerc * max(traces_normalized);  
    thresholdValuesForBinarise = thresholdCalPerc * max(traces_normalized(:));  
    
    %Step 2: Binarise the normalised traces based on the threshold determined
    %for each cell in matrix [thresholdValuesForBinarise].
    binarisedTraces = traces_normalized > thresholdValuesForBinarise;
    binarisedTracesMatrix = double(binarisedTraces);
    if plotRaster
        rasterplot(find(binarisedTracesMatrix), ...
            size(binarisedTracesMatrix, 2), size(binarisedTracesMatrix, 1));
    end
end

% Function returns the minimum and maximum values of each cell from the
% traces provided as input. Returns a row vector having each value of time
% array [min, max]
%
% Params: (traces_normalized)
%
function matrixBoundsForSimul = getCellWiseMixMaxBounds(traces_normalized)
    matrixBoundsForSimul = [];
    for i = 1:size(traces_normalized, 2)
        tempCell = [];
        tempCell(end + 1) = min(traces_normalized(:, i));
        tempCell(end + 1) = max(traces_normalized(:, i));
        matrixBoundsForSimul(:, i) = tempCell;
    end
end

% Function generates matrices with random generated values based on the limits provided in the min-max
% vector for each cell. The generated matrices are stored in a Cell
% datastructure and returned.
%
% Params: (numberOfDummyMatrices, traces_normalized, matrixBoundsForSimul)
%          numberOfDummyMatrices - Number of simulated matrices to generate.
%
function dummyCellMatrixList = generateSimulatedTracesBasedOnCellBounds(numberOfDummyMatrices, traces_normalized, matrixBoundsForSimul)
    dummyCellMatrixList = cell(numberOfDummyMatrices, 1);
    for i = 1: numberOfDummyMatrices
        x = zeros(size(traces_normalized));
        for j = 1: size(matrixBoundsForSimul, 2)
            boundsColumn = matrixBoundsForSimul(:, j);
            a = boundsColumn(1, :);
            b = boundsColumn(2, :);
            x(:, j) = a + (b-a).*rand(size(traces_normalized, 1), 1);
        end
        dummyCellMatrixList{i} = x; 
    end
end

% Function extracts the timeline events in 'COUNT' terms for each cell from a 
% binarised matrix and stores in an array.
% Active events are kept as positive integers and Non-active events are
% kept as negative integers for identification.
%
% Params: (binarisedTraces)
%
function eventTimeRowRecord = getCellTimelineEventsFromBinarisedTrace(binarisedTraces)
    eventTimeRowRecord = cell(size(binarisedTraces, 2), 1);
    for i = 1:size(binarisedTraces, 2)
        isOn = 0;
        isOff = 0;
        cellTempArray = [];
        cellRep = binarisedTraces(:, i);
        for j = 1:size(cellRep, 1)
            if (cellRep(j, :) == 1)
                if (isOff > 0)
                    cellTempArray(end + 1) = -isOff;
                    isOff = 0; %reset timeline tracker for In-Active
                end
                isOn = isOn + 1;
            end
            if (cellRep(j, :) == 0)
                if (isOn > 0)
                    cellTempArray(end + 1) = isOn;
                    isOn = 0; %reset timeline tracker for Active
                end
                isOff = isOff + 1;
            end
            % Handle end sequence
            if (j == size(cellRep, 1))
                if (isOff > 0)
                    cellTempArray(end + 1) = -isOff;
                    isOff = 0; %reset timeline tracker for In-Active
                end
                if (isOn > 0)
                    cellTempArray(end + 1) = isOn;
                    isOn = 0; %reset timeline tracker for Active
                end
            end
        end
        eventTimeRowRecord{i} = cellTempArray;
    end
end

% Function generates trace matrices by shuffling the timeline events. The generated matrices are stored in a Cell
% datastructure and returned. This function requires the
% 'eventTimeRowRecord' to be captured as a pre-curser.
%
% Params: (numberOfDummyMatrices, binarisedTraces, eventTimeRowRecord)
%          numberOfDummyMatrices - Number of simulated matrices to generate.
%
function resultCombo = generateSimulatedTracesFromEventRecord(numberOfDummyMatrices, binarisedTraces, eventTimeRowRecord)
    tempCellMatrixList = cell(numberOfDummyMatrices, 1);
    eventTimeRowRecordShuffled = cell(size(binarisedTraces, 2), 1);
    resultCombo = cell(2, 1);
    for i = 1:numberOfDummyMatrices
        tempMatrice = zeros(size(binarisedTraces)); % Create an empty matrice of original trace size
        for cellIndex = 1:numel(eventTimeRowRecord)
            tempArray = eventTimeRowRecord{cellIndex}; % Get timeline sequence of each cell [Array format]
            tempArray2 = tempArray(:, randperm(size(tempArray, 2))); % Shuffle the timeline sequence of each cell
            eventTimeRowRecordShuffled{cellIndex} = tempArray2;
            rowIndex = 1;
            % Create each cell column based on shuffled timeline sequence
            % by iterating over each value of the timeline array
            for j = 1:size(tempArray2, 2)
                timelineVal = tempArray2(:, j);
                absTimeVal = abs(timelineVal); % convert to absolute for iteration purpose
                for k = 1:absTimeVal
                    if (timelineVal > 0) % Active filter
                        tempMatrice(rowIndex, cellIndex) = 1;
                    end
                    if (timelineVal < 0) % Non-active filter
                        tempMatrice(rowIndex, cellIndex) = 0;
                    end
                    rowIndex = rowIndex + 1; % increment row index to add incrementaly to the column values
                end
            end
        end
        tempCellMatrixList{i} = tempMatrice;
    end
    resultCombo{1} = tempCellMatrixList;
    resultCombo{2} = eventTimeRowRecordShuffled;
end

% 
% Function normalises the traces and returns.
%
function tracesCombo = getNormalisedTraces(traces, sr, dt, tauDecay, useRomanoMethod)
    tracesCombo = cell(2, 1);
    nFrames = size(traces,1);
    fS = sr; % Sampling frequency (HZ)   
    tWin = 5;% The first seconds for calculating the F0
    bWin = round(tWin * fS);    
    nWin = floor(nFrames / bWin);
    
    % This loop create means every bWin frames for the first nwin calculating the  F0
    for k = 1:nWin
        mu(k,:) = mean(traces((1:bWin)+(k-1)*bWin,:));
    end
    
    % Finding the lowest min (mean)
    if useRomanoMethod
        [fluoTraces, F0, smoothBaseline, deletedCells] = SanityTest(traces, fS, tauDecay);
        dFF = (fluoTraces - F0)./F0;
    else
        F0 = min(mu);
        dFF = traces./F0 - 1; % DF over F
    end
    

    % Filter to smooth the calcium trace
    nFilt = 100;
    %fS = 10;
    fCut = 2;
    cutRatio = fCut/(fS/2); % Cut Ratio is the size of the frequncy we want to remove
    
    % If cutRatio is lower than 1 (fs>4) fiter is applied to smooth the signal
    % if cut ratio is =< 1 (fs <= 4) fitler isn't applied.
    if useRomanoMethod
        dFF_filt = dFF;
    else
        if (cutRatio < 1)
            dFF_filt = filtfilt(fir1(nFilt,fCut/(fS/2)),1,dFF);
        else
            dFF_filt(:,i) = dFF(:,i);
            dFF_filt = dFF;
        end
    end
    
    save('basal1_allTraces','dFF','dFF_filt')
    
    traces_normalized = dFF_filt;

    time = 0:dt:(size(traces_normalized, 1) - 1) * dt;

    Baseline_total = zeros(size(traces_normalized, 1), size(traces_normalized, 2));

    crit = zeros(size(traces_normalized, 1), size(traces_normalized, 2));

    deconv = zeros(size(traces_normalized, 1), size(traces_normalized, 2));
    
    for i = 1:size(traces_normalized, 2)
        Baseline_total(:, i) = prctile(traces_normalized(:, i), 10) * ones(size(traces_normalized, 1), 1);
        Baseline = prctile(traces_normalized(:, i), 10);
        [peaks_elec_P, criterion] = pickpeaks(traces_normalized(:, i), 0.08, 0);
        crit(:, i) = criterion;
        threshold = 5; % prctile(criterion,100);
        [pks, dep, peaks_elec_C, didx] = peakdet(criterion, threshold, 'threshold', 300); % perform peak-dep detection
        transient{i} = peaks_elec_C;
    end
    tracesCombo{1} = traces_normalized;
    tracesCombo{2} = transient;
    tracesCombo{3} = time;
end

function aggrThreshold = plotHistogramWithStdAndMean(matrice, titleText, additionalPlot)
    histogram(matrice);
    grid on;
    % Compute mean and standard deviation.
    mu = mean(mean(matrice));
    sigma = mean(std(matrice));
    aggrThreshold = mu + 2 * sigma;
    % Indicate those on the plot.
    xline(mu, 'Color', 'g', 'LineWidth', 2);
    xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    %xline(meanTwo - sigma, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
    xline(aggrThreshold, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
    if (exist('additionalPlot', 'var'))
        xline(additionalPlot, 'Color', 'black', 'LineWidth', 3, 'LineStyle', '--');
    end
    title(titleText, 'FontSize', 15);
end

function plotHeatmapWithSignificance(originalPerAggrMatrice, threshold, titleText)
    hm = heatmap(originalPerAggrMatrice, 'ColorMethod','median');
    xh = hm.Colormap > threshold;
    [rowIdx, colIdx] = find(xh == 1);
    hm.Colormap(rowIdx, colIdx) = 0.5;
    hm.Title = titleText;
end

function [deltaFoF, mu, sigma] = estimateBaselineNoise(tracesNormalised)
    numCells=size(tracesNormalised,2);
    numFrames=size(tracesNormalised,1);
    
    % We calculate the ROI's baseline noise by fitting a
    % gaussian to the distribution (local density estimation) of negative values of deltaFoF
    
    for numcell = 1:numCells
        dataCell = tracesNormalised(:,numcell);
        [smoothDist, x] = ksdensity(dataCell);
        [valuePeak, indPeak] = max(smoothDist);
        xFit = x(1:indPeak);
        dataToFit = smoothDist(1:indPeak) / numFrames;
        [sigma(numcell), mu(numcell), A] = CustomGaussfit(xFit', dataToFit);
        
        if ~isreal(sigma(numcell))
            dev = nanstd(dataCell);
            outliers = abs(tracesNormalised) > 2 * dev;
            
            deltaF2 = dataCell;
            deltaF2(outliers) = NaN;
            sigma(numcell) = nanstd(deltaF2);
            mu(numcell) = nanmean(deltaF2);
            
        end

        distFit = A*exp(-(x-mu(numcell)).^2./(2*sigma(numcell)^2));
    end
    deltaFoF = bsxfun(@minus, tracesNormalised, mu);
end