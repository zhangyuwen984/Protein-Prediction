function WriteCompStatistic( resFolder )

global USPEX_STRUC
global ORG_STRUC


numIons    = ORG_STRUC.numIons;
statistic  = USPEX_STRUC.statistic;
convexHull = USPEX_STRUC.GENERATION(end).convex_hull;
lenOfBlocks= length(USPEX_STRUC.POPULATION(1).numBlocks);


lenOfComp = 1:size(statistic,1);

composition = statistic(:,1:lenOfBlocks);
ratioOfComp = statistic(:,lenOfBlocks+1:2*lenOfBlocks);

summaryInfo = [];

fp = fopen([resFolder '/compositionStatistic'],'w');
fprintf(fp, '       Comp/Ratio        Total     Random Heredity Mutation Seeds  COPEX  Best/Convex\n');
for whichComp = 1:length( lenOfComp )
    if lenOfComp(whichComp) == 0
        continue;
    else
        sameRatioComp = whichComp;
    end
    
    currentRatio = ratioOfComp( whichComp, : );
    sumStatistic = statistic(whichComp, :);
    
    for i = whichComp+1:length( lenOfComp )
        if abs( ratioOfComp(i,:) - currentRatio ) < 1e-3
            sameRatioComp = [sameRatioComp, i];
            sumStatistic = sumStatistic + statistic(i,:);
            sumStatistic(lenOfBlocks+1:2*lenOfBlocks) = currentRatio;
            lenOfComp(i) = 0;
        end
    end
    
    
    % Summay of the statistic
    fprintf(fp, ' [ ');
    fprintf(fp, ' %6.4f ', currentRatio);
    fprintf(fp, ' ]' );
    fprintf(fp, '%5d(%3d)', sumStatistic(2*lenOfBlocks+1)-sumStatistic(end),sumStatistic(2*lenOfBlocks+1));
    fprintf(fp, ' %7d', sumStatistic(2*lenOfBlocks+2:end));
    fprintf(fp,'\n');
    % Detail of the statistic
    for k = 1:length( sameRatioComp )
        ID = sameRatioComp(k);
        fprintf(fp, '    ');
        fprintf(fp, ' %4d ', composition(ID,1:lenOfBlocks) );
        fprintf(fp, '     ');
        fprintf(fp, '%5d(%3d)',statistic(ID, 2*lenOfBlocks+1)-statistic(ID, end), statistic(ID, 2*lenOfBlocks+1) );
        fprintf(fp, ' %7d',  statistic(ID, 2*lenOfBlocks+2:end));
        fprintf(fp, '\n');
    end
    
    % Data collect for bar plot
    summaryInfo = [summaryInfo; sumStatistic];
end

fclose(fp);


avgStructure = sum( summaryInfo(:,lenOfBlocks*2+1)- summaryInfo(:,end) )/size(summaryInfo,1);


%
% 1D energy statistic figure
%
if lenOfBlocks==1
    try
        numX = ceil(length(USPEX_STRUC.POPULATION)/USPEX_STRUC.POPULATION(end).gen);
        
        h1 = figure;
        set(gcf,'Visible','off');
        len = 0;
        fitness = zeros(1,len);
        for i = 1:length(USPEX_STRUC.POPULATION)
            if ~strcmp(USPEX_STRUC.POPULATION(i).howCome, 'keptBest')
                fitness(i) = USPEX_STRUC.POPULATION(i).Fitness;
                len = len+1;
            else
                fitness(i) = 1.0E+5;
            end
        end
        maxFitness = max(fitness(find(fitness<1.0e+4)));
        minFitness = min(fitness(find(fitness<1.0e+4)));
        dFitness = (maxFitness-minFitness)/numX ;
        fitnessX= minFitness: (maxFitness-minFitness)/numX :maxFitness+dFitness;
        fitnessY= zeros(1,numX+2);
        for i = 1:numX+1
            fitnessY(i) = length( find(fitness>=fitnessX(i) & fitness< fitnessX(i+1)) );
        end
        bar(fitnessX, fitnessY);
        xlabel('Fitness');
        ylabel('Number of structures');
        %axis([ min(ax)-0.01*max(ax) max(ax)*1.01 -maxTot*0.05 maxTot ])
        print(h1,'-dpdf' ,        [resFolder '/fitnessStatistic.pdf']);
    catch
        %====== Do nothing
    end
end

%
% 2D composition statistic figure
%
if  lenOfBlocks==2
    try
    h1 = figure;
    set(gcf,'Visible','off');
    
    ax = summaryInfo(:,lenOfBlocks*2);
    ay = summaryInfo(:,lenOfBlocks*2+2:end);
    ay(:,end)=0;
    maxTot = max(summaryInfo(:,lenOfBlocks*2+1)- summaryInfo(:,end));
    
    for i = 1:length( ax )
        if summaryInfo(i, lenOfBlocks*2+1)<0.001
            ay(i,end) = -maxTot*0.05;
        end
    end
    ay(:,end+1)=0;
    bar(ax,ay,'stacked')
    
    lineX=0;
    lineY=0;
    
    for i = 1:size( convexHull,1 )
        lineX(i) = convexHull(i,2)/sum( convexHull(i,1:2) );
        lineY(i) = convexHull(i,3)*sum(convexHull(i,1:2)*numIons)/sum( convexHull(i,1:2) );
    end
    [lineX, seq] = sort( lineX );
    lineY=lineY(seq)-lineY(1);
    
    nlineY = length(lineY);
    for i = 1:nlineY
        lineY(i) = lineY(i)-lineY(nlineY)*lineX(i)+1;
    end
    
    % maxY = maxTot
    % minY = maxTot*1/3
    lineY = ( lineY/max(lineY) )*maxTot*2/3+maxTot*1/3;
    hold on;
    plot(lineX, lineY, 'k-');
    scatter(lineX,lineY,'MarkerEdgeColor','k','MarkerFaceColor','b');
    
    plot([-0.1 1.1], [avgStructure avgStructure],'--')
    
    colormap(hsv);
    hcb = colorbar('YTick', [1:1:6],'YTickLabel',{'Random','Heredity','Mutation',' Seeds', 'COPEX', 'Never Sampled'});
    set(hcb,'YTickMode','manual')
    %legend('Random','Heredity','Mutation','Seeds', 'COPEX', 'No Structure', 'convexHull','Location','Best');
    
    xlabel([ blockSymbol(ORG_STRUC.atomType, ORG_STRUC.numIons(1,:)), '+',  blockSymbol(ORG_STRUC.atomType, ORG_STRUC.numIons(2,:)),]);
    ylabel('Numbers of Structures ');
    axis([ min(ax)-0.01*max(ax) max(ax)*1.01 -maxTot*0.05 maxTot ])
    print(h1,'-dpdf' ,        [resFolder '/compositionStatistic.pdf']);
    catch
        %====== Do nothing
    end
end

%
% 3D composition statistic figure
%
if  lenOfBlocks==3
    try
    h1 = figure;
    set(gcf,'Visible','off');
    
    a=[-sqrt(3)/2,-1/2];  % 1 0 0
    b=[0 , 1 ];           % 0 1 0
    c=[sqrt(3)/2, -1/2];  % 0 0 1
    
    matrixPLOT=[a;b;c];
    
    line([a(1),b(1)],[a(2),b(2)]);
    line([a(1),c(1)],[a(2),c(2)]);
    line([c(1),b(1)],[c(2),b(2)]);
    axis([-1 1 -0.8 1.2]);
    axis square
    %colormap(hsv);
    hold on
    
    numInd =zeros(1,size(summaryInfo,1));
    fitness=zeros(1,size(summaryInfo,1));
    
    for i = 1:size(summaryInfo,1)
        numInd(i) =sum(summaryInfo(i,2*lenOfBlocks+1))-sum(summaryInfo(i,end));
        fitness(i)=100000;
    end
    
    maxInd = max(numInd);
    
    
    for i = 1:length(USPEX_STRUC.POPULATION)
        comp = USPEX_STRUC.POPULATION(i).numBlocks;
        ratioOfComp = comp/norm(comp);
        for j = 1:size(summaryInfo,1)
            ratio = summaryInfo(j,1:lenOfBlocks)/norm( summaryInfo(j,1:lenOfBlocks) );
            if norm(ratioOfComp-ratio)<1e-4
                if fitness(j)> USPEX_STRUC.POPULATION(i).Fitness;
                    fitness(j) = USPEX_STRUC.POPULATION(i).Fitness;
                    break
                end
            end
        end
    end
    
    maxFitness = 0.15;
    [result,seq0]= find(fitness>99999);
    fitness(seq0) = 0;
    
    averageFitness= sum( fitness )/(length(fitness)-length(seq0));
    [result,seq]= find(fitness>averageFitness);
    fitness(seq)= averageFitness;
    fitness = fitness/averageFitness*maxFitness*2/3;
    fitness(seq0)= maxFitness;
    [result,seq]=sort(numInd,'descend');
    for i = seq
        pos=summaryInfo(i,lenOfBlocks+1:2*lenOfBlocks)*matrixPLOT;
        
        mkSize=numInd(i)*6;
        if mkSize > 0
        end
        fitPlot = fitness(i);
        if mkSize==0;
            plot(pos(1),pos(2),'r.','MarkerSize',5);
        else
            if fitPlot>1e-4
                scatter(pos(1),pos(2),mkSize,fitPlot,'filled','MarkerEdgeColor','k');
            else
                scatter(pos(1),pos(2),mkSize,fitPlot,'filled','MarkerEdgeColor','r');
            end
        end
    end
    
    %scatter(-0.8,-0.7,  6,                 0.075, 'filled','MarkerEdgeColor','k');
    scatter(-0.82,-0.7, 1*6,  0.0, 'filled','MarkerEdgeColor','k');
    scatter(-0.42,-0.7, round(maxInd/3)*6,  0.025, 'filled','MarkerEdgeColor','k');
    scatter(-0.02,-0.7, round(maxInd/2)*6,  0.00,  'filled','MarkerEdgeColor','r');
    scatter( 0.38,-0.7, round(maxInd*2/3)*6,0.05,  'filled','MarkerEdgeColor','k');
    scatter( 0.78,-0.7, maxInd*6,           0.075, 'filled','MarkerEdgeColor','k');
    %text(-0.8,-0.7,' =1' )
    text(-0.8,-0.7,['  =' num2str(1)] );
    text(-0.4,-0.7,['  =' num2str(round(maxInd/3))] );
    text( 0.0,-0.7,['  Stable '] );
    text( 0.4,-0.7,['  =' num2str(round(maxInd*2/3))]);
    text( 0.8,-0.7,['  =' num2str(maxInd)] );
    
    text(-sqrt(3)/2, -0.58, blockSymbol(ORG_STRUC.atomType, ORG_STRUC.numIons(1,:))  );
    text( sqrt(3)/2, -0.58, blockSymbol(ORG_STRUC.atomType, ORG_STRUC.numIons(3,:))  );
    text( 0 , 1.08,         blockSymbol(ORG_STRUC.atomType, ORG_STRUC.numIons(2,:))  );
    
    caxis([0 maxFitness]);
    hcb = colorbar('YTick',0:maxFitness/6:maxFitness,'YTickLabel',...
        {'0 (Stable)','0.25','0.5','0.75', '>1 AverageFitness','Nerver Sampled',''});
    
    axis off
    print(h1,'-dpdf' ,        [resFolder '/compositionStatistic.pdf']);
    catch
        %====== Do nothing
    end
end


%--------------------------
function str = blockSymbol(atomType, numIons)


str='';
for i = 1:length(numIons)
    if numIons(i)>0
        if numIons(i)>1
            label = ['{' num2str(numIons(i)) '}'];
            str = [str, megaDoof(atomType(i)), '_',label];
        else
            str = [str, megaDoof(atomType(i))];
        end
        
    end
end
