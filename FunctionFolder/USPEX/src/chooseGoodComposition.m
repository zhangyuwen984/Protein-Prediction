function ind = chooseGoodComposition(tournament, POP)


count = 0;
goodComposition = 0;

maxTime = length(POP)*10;

while ~goodComposition
    
    toMutate = find (tournament>RandInt(1,1,[0,max(tournament)-1]));
    ind = toMutate(end);
    
    numBlocks = POP(ind).numBlocks;
    if ~isempty(numBlocks)  %% for fix composition condition
        goodComposition = CompositionCheck(numBlocks);
    else
        goodComposition = 1;
    end
    
    count = count +1;
    if count > maxTime  
        ind = -1;
        break;
    end
end