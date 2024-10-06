function [toDelete] = checkCells(fluoTraces,baseline,data)

%v = ver;
%parOn=any(strcmp('Parallel Computing Toolbox', {v.Name}));
parOn = 0; 
cutOffIntensity = 0;
cutOffPixels = 0;
%parOn=1; %UNCOMMENT THIS FOR PARALLELIZATION
if parOn
    numCores = feature('numcores');
    
    isOpen = parpool('size') > 0;
    if ~isOpen
        parpool('open','local', round(numCores * 0.8))
    end
    
end


numCells=size(baseline,2);
numFrames=size(baseline,1);
baselineMeans=zeros(size(baseline));

if parOn
    parfor i=1:numCells
        baselineMeans(:,i)=runline(baseline(:,i),min(2000,round(numFrames/4)),1);
    end
else
    for i=1:numCells
        baselineMeans(:,i)=runline(baseline(:,i),min(2000,round(numFrames/4)),1);
    end
end
baselineVariations=zscore((baseline-baselineMeans)./baselineMeans);

% Check for ROIs with few pixels
thresh = 0;
deviations=baselineVariations<thresh;
toDeleteArtifacts=find(logical(sum(deviations,1)));
toDeleteArtifacts=unique([toDeleteArtifacts find(isnan(sum(fluoTraces,1)))]);
for i=1:numCells
    pixs(i)=length(data(:, i));
end
toDeletePixs=find(pixs < cutOffPixels);


% Check for ROIs with weak baseline
toDeleteDim1=[];
for j=1:numCells
    dataSlice=fluoTraces(:,j)-baseline(:,j);
    
    [counts,x]=hist(dataSlice,min(round(length(dataSlice)/20),100));
    counts=smooth(counts);
    if any(isnan(x))
        toDeleteDim1=[toDeleteDim1 j];
        continue
    end

    
    
    [valueCenter,indCenter]=max(counts);
    
    if counts(1)>valueCenter * cutOffIntensity/100
        toDeleteDim1=[toDeleteDim1 j];
    end
    
end

toDeleteDim2=[];
for j=1:numCells
    if any(baseline(:,j)==0)
        toDeleteDim2=[toDeleteDim2 j];
    end
end

toDeleteDim=union(toDeleteDim1,toDeleteDim2);
toDelete=union(union(toDeleteDim,toDeleteArtifacts),toDeletePixs);
okCells=setdiff(1:numCells,toDelete);