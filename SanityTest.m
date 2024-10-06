function [fluoTraces,F0,smoothBaseline,deletedCells] = SanityTest(traces, fps, tauDecayVal)

parOn = 0; 
tauDecay = tauDecayVal;
if parOn
    numCores = feature('numcores');
    isOpen = parpool('size') > 0;
    if ~isOpen
        parpool('open','local', round(numCores * 0.8))
    end
    
end

fluoTraces = traces;
% Calculation the baseline in a running time window whose length is the maximum of 15 s or 40 time decays

twdw = max(15, 40 * tauDecay); 
wdw = round(fps * twdw);
numCells = size(fluoTraces, 2);
numFrames = size(fluoTraces, 1);
smoothBaseline = zeros(size(fluoTraces));

disp('1.1 Calculating fluorescence baseline of ROIs...')
if parOn
      if numFrames > 2*wdw
        parfor j = 1:numCells
            dataSlice = fluoTraces(:,j);
            temp = zeros(numFrames-2*wdw,1);
            for i = wdw+1:numFrames-wdw
                temp(i-wdw) = prctile(dataSlice(i-wdw:i+wdw),8);
            end
            smoothBaseline(:,j) = [temp(1)*ones(wdw,1) ; temp; temp(end)*ones(wdw,1)];
            smoothBaseline(:,j) = runline(smoothBaseline(:,j),wdw,1);
        end
        
    else
        parfor j=1:numCells
            smoothBaseline(:,j) = [ones(numFrames,1)*prctile(fluoTraces(:,j),8)];
        end
    end
else
    if numFrames > 2*wdw
        for j=1:numCells
            dataSlice=fluoTraces(:,j);
            temp=zeros(numFrames-2*wdw,1);
            for i=wdw+1:numFrames-wdw
                temp(i-wdw)=prctile(dataSlice(i-wdw:i+wdw),8);
            end
            smoothBaseline(:,j)=[temp(1)*ones(wdw,1) ; temp; temp(end)*ones(wdw,1)];
            smoothBaseline(:,j)=runline(smoothBaseline(:,j),wdw,1);
        end
        
    else
        for j=1:numCells
            smoothBaseline(:,j) = [ones(numFrames,1)*prctile(fluoTraces(:,j),8)];
        end
    end
    
end
 
indStart=max(round(fps * 1), 1);
indEnd=min(round(fps * size(fluoTraces, 1)), numFrames);
F0=repmat(mean(fluoTraces(indStart:indEnd,:)),size(fluoTraces,1),1);
% We look for cells that are not appropiate and remove them.

[deletedCells]= 0; %checkCells(fluoTraces,smoothBaseline, fluoTraces);
toKeep=setdiff(1:numCells,deletedCells);

fluoTraces=fluoTraces(:,toKeep);
smoothBaseline=smoothBaseline(:,toKeep);
F0=F0(:,toKeep);