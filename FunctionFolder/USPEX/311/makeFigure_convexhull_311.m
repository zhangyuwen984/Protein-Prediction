function makeFigure_convexhull_311( generation, convex_hull, comp, points, points_new)

%  This function is extracted from extenededConvexHull_3XX.m
%  To make figures for 2D and 3D(in the future) convex_hull

global ORG_STRUC

if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end

atomType = ORG_STRUC.atomType;
numIons  = ORG_STRUC.numMols;

x  = points(1,:);
y  = points(2,:);

x_ = points_new(1,:);
y_ = points_new(2,:);


if size(comp,2) == 2
    A0=convex_hull(1,end-1); %ev/atom
    B0=convex_hull(2,end-1); %eV/atom
end


try
    if size(numIons,1) == 2
        % add point to convex hull array, for plotting
        inConvexHull = zeros(2, size(convex_hull,1));
        for i = 1:size(convex_hull,1)
            N_block = convex_hull(i,1:2);
            N_atom = N_block*numIons;
            E = convex_hull(i,3);
            inConvexHull(1,i) = N_block(2)/sum(N_block); % y/(x+y) by block
            inConvexHull(2,i) = (E*sum(N_atom) - A0*sum(numIons(1,:))*N_block(1) - B0*sum(numIons(2,:))*N_block(2)) / sum(N_block); % eV/block
        end
        
        [nothing, chRanking] = sort(inConvexHull(1,:)); % we want line to go from left to right :)
        h1 = figure;
        set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
        
        if OctaveMode == 0
            scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
        else
            scatter(x,y);
        end
        
        xlabel('Composition');
        ylabel('Enthalpy of formation (eV/block)');
        line(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)));
        hold on
        box on;
        
        if OctaveMode == 0
            scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)),'MarkerEdgeColor','g','MarkerFaceColor','k');
            scatter(x_(1:end), y_(1:end), 'MarkerEdgeColor','k','MarkerFaceColor','r');
        else
            scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)));
            scatter(x_(1:end), y_(1:end));
        end
        
        %   print(h1,'-dtiff','-r300','extendedConvexHull.tif');
        print(h1,'-dpdf' ,        'extendedConvexHull.pdf');
        print(h1,'-dpdf' ,       ['generation' num2str(generation)  '/extendedConvexHull.pdf']);
        hold off
    end
catch
end
