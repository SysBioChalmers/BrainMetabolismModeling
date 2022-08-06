%This is copied from GeckoLight. The difference is that it splits up the costs into the 4 categories
function [MWDivKcatsRes,rxnWcLevels] = GetMWAndKcatNeurons(uniprots, kcats, wcLevels, swissprot, standardMW, mtGenes, pureCytRxns, mitoRxns)

swissProtHashMap = containers.Map(swissprot(:,1),1:length(swissprot(:,1)));

[nRxns,~] = size(kcats);
MWDivKcatsRes = NaN(nRxns,4);
rxnWcLevels = NaN(nRxns,1);
[~,n]   = size(uniprots);
y       = 0;
N1 = 0;

indOther = 1;
indCyto = 2;
indMito = 3;
indMT = 4;%Mitochondrial genes, they are produced locally in the mitochondria and are penalized less.


for i = 1:length(MWDivKcatsRes)

    %split up uniprots
    for j = 1:n
        if (isempty(uniprots{i,j}))
            break;
        end
        %Update vector enzymes and count the number of isozymes (x):
        if kcats(i,j) > 0
            uniprots{i,j} = strsplit(uniprots{i,j},' ');
        end
    end

    if ~isempty(uniprots{i,1})
        x = 0;
        MWDivKcats = NaN(n, 5);
        wcLev = NaN(n, 1);
        for j = 1:n
            if isempty(uniprots{i,j}) %this means we have reached the end of the list
                break;
            end
            if ~isempty(uniprots{i,j}) && kcats(i,j) > 0 
                x         = x+1;
                kcat = kcats(i,j);
                wcLev(j) = wcLevels(i,j);
                protNames = uniprots{i,j};
                if length(protNames) > 1
                    breakpt = 1;
                end

                %so, we reason that if the MW is not found, use let's a standard value
                %if a reaction requires 3 proteins, and we only find MWs for 2
                %of them, a standard value is a better guess than zero
                MWs = repmat(standardMW, length(protNames),1);
                %MWs = NaN(length(protNames),1);%this is temporary to compare with the other EC model
                %find MW for the proteins

                indices = values(swissProtHashMap,protNames); %so, all must be found; I think they should be found? If not, replace with the loop below
                isMt = false(length(protNames),1);
                for k = 1:length(protNames)
                    if swissprot{indices{k},5} ~= 0
                        MWs(k) = swissprot{indices{k},5}/1000;	%g/mmol
                        isMt(k) = ismember(swissprot{indices{k},3}, mtGenes);
                    end
                end
                
                %Here, we instead of just summing up the MWs/kcat, we also split it up into different categories 
                
                %these checks are temporary, just remove later
                isMt(isnan(MWs)) = [];
                MWs(isnan(MWs)) = [];
                if (~isempty(MWs))
                    sumMW = sum(MWs);
                    mtFrac = sum(MWs(isMt))/sumMW;
                    MWDivKcats(j, indMT) = mtFrac*sumMW/kcat;
                    
                    MWDivKcats(j, indOther) = 0;
                    MWDivKcats(j, indCyto) = 0;
                    MWDivKcats(j, indMito) = 0;
                    MWDivKcats(j, 5) = sumMW/kcat; %used further down, here to simplify that code
                    
                    %if mtFrac > 0
                    %    brk = 1;
                    %end
                    
                    if pureCytRxns(i)
                        MWDivKcats(j, indCyto) = (1-mtFrac)*sumMW/kcat;
                    elseif mitoRxns(i)
                        MWDivKcats(j, indMito) = (1-mtFrac)*sumMW/kcat;
                    else
                        MWDivKcats(j, indOther) = (1-mtFrac)*sumMW/kcat;
                    end
                end
            end
        end
        if (x > 0)
            if (any(~isnan(MWDivKcats(1:x))))
                %Apply some extra logic here due to the high uncertainty in
                %wildcard searches for EC numbers
                %So, if we have isozymes with a match on EC number, we
                %ignore the isozymes with wildcard match (i.e. wc > 0) since we don't
                %really trust those values
                lowestWC = min(wcLev(~isnan(wcLev)));
                if lowestWC == 0
                    [~, minInd] = min(MWDivKcats((~isnan(MWDivKcats(:,5))) & wcLev == 0,5));
                    MWDivKcatsRes(i,:) = MWDivKcats(minInd,1:4);
                    %[MWDivKcatsRes(i), minInd] = min(MWDivKcats((~isnan(MWDivKcats)) & wcLev == 0));
                else
                    [~, minInd] = min(MWDivKcats(~isnan(MWDivKcats(:,5)),5));
                    MWDivKcatsRes(i,:) = MWDivKcats(minInd,1:4);

                    %[MWDivKcatsRes(i), minInd] = min(MWDivKcats(~isnan(MWDivKcats)));
                end
                rxnWcLevels(i) = wcLev(minInd);
            end
        end
    end
    if rem(i,100) == 0 || i == nRxns 
        fprintf('.')
    end

end

fprintf(' Done!\n')

end