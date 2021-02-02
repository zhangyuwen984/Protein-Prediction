function update_USPEX_GENERATION(IND, fitness, type)

%To update all the necessary items at the end of each gen to USPEX_STRUC
%Output the BESTINDIVIDUALS info
%Lastly updated by Qiang Zhu (2014/02/18)

global USPEX_STRUC
global POP_STRUC


USPEX_STRUC.GENERATION(POP_STRUC.generation).quasiEntropy = oldQuasiEntropy();
USPEX_STRUC.GENERATION(POP_STRUC.generation).Fitness = min(fitness);

if type == 2  %bulk varcomp
   USPEX_STRUC.GENERATION(POP_STRUC.generation).convex_hull= POP_STRUC.convex_hull;
   USPEX_STRUC.GENERATION(POP_STRUC.generation).composEntropy= composEntropy();
end

ID = [];
for i=1:length(IND)
    if IND(i) > 0
       ID = [ID  POP_STRUC.POPULATION(IND(i)).Number];
    end
end
USPEX_STRUC.GENERATION(POP_STRUC.generation).BestID = ID;
USPEX_STRUC.GENERATION(POP_STRUC.generation).ID = length(USPEX_STRUC.POPULATION); %New field to count

[nothing, ranking] = sort(POP_STRUC.DoneOrder);
for i = 1 : length(fitness)
    if POP_STRUC.DoneOrder(ranking(i)) > 0
       ID = POP_STRUC.POPULATION(ranking(i)).Number;
       USPEX_STRUC.POPULATION(ID).Fitness = fitness(ranking(i));
    end
end

%safesave([resFolder '/USPEX.mat'], USPEX_STRUC); %We don't save it here, but after everything is all set
