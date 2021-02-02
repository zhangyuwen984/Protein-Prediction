function [ToStop, Message] = StopRun(fitness, generation, numGenerations, stopCrit, stopFitness)
Message = '';
ToStop = 0;
if generation >= numGenerations
    ToStop = 1;
    Message = ['Maximum Generation number achieved: ' num2str(numGenerations)];
elseif ~isempty(stopFitness)
    if min(fitness) <= stopFitness
        Message = ['Known ground state achieved: ' num2str(stopFitness) 'eV'];
        ToStop = 1;
    end
elseif generation >= stopCrit
    same_in = SameBest();
    if same_in == stopCrit - 1
        Message = ['Halting criteria achieved: ' num2str(stopCrit)];
        ToStop = 1;
    end
end

