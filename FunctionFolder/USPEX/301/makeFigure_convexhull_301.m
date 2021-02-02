function makeFigure_convexhull_301(generation, convex_hull, comp, points, points_new, Energy_base)
%  This function is extracted from extenededConvexHull_3XX.m
%  To make figures for 2D and 3D(in the future) convex_hull
global ORG_STRUC
if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end
atomType = ORG_STRUC.atomType;
numIons  = ORG_STRUC.numIons;
x  = points(1:end-1,:);
y  = points(end,:);
x_ = points_new(1:end-1,:);
y_ = points_new(end,:);
if size(numIons,1) == 3
   disp('triangle diagram')
   plot_tenary_diagram(atomType, numIons, convex_hull, x', y', Energy_base);
elseif size(numIons,1) == 2
    % add point to convex hull array, for plotting
    inConvexHull = zeros(2, size(convex_hull,1));
    ConvexHull_atoms = zeros(size(convex_hull,1), size(atomType,2));
    
    for i = 1:size(convex_hull,1)
        N_atom  = convex_hull(i,1:end-2);
        E = convex_hull(i,end-1);  %eV/atom
        [inConvexHull(1,i),inConvexHull(2,i)]  = Get_XY(numIons, N_atom, Energy_base, E);
        ConvexHull_atoms(i,:) = N_atom;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        h1 = figure;
        set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
        
        if OctaveMode == 0
            scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
        else
            scatter(x,y);
        end
        
        hold on;
        box on;
        [nothing, chRanking] = sort(inConvexHull(1,:)); % we want line to go from left to right :)
        line(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)), ...
            'LineWidth', 1.0, 'Color', 'k');
        
        if OctaveMode == 0
            scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)),...
                'MarkerEdgeColor','g','MarkerFaceColor','k');
        else
            scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)));
        end
        
        dy = (max(y)-min(y))/10;
        xlim([min(x)-0.01, max(x)+0.01]);  %to Prevent the range is too big
        ylim([min(y)-dy,   max(y)+dy]);  %to Prevent the range is too big
        makeCaption(atomType, numIons);
        set(gca,'FontSize',12,'linewidth',1.0);
        print(h1,'-dpdf' ,        'extendedConvexHull.pdf');
		position = inConvexHull;
		position(2,:) = position(2,:);
		Label = makeLabel(atomType, ConvexHull_atoms, position);
		for loop1 = 1:size(ConvexHull_atoms, 1)
			text(ConvexHull(1,loop1), ConvexHull(2,loop1)-dy, label{loop1}, 'HorizontalAlignment', 'Center');
        end
        
        if OctaveMode == 0
            scatter(x_(1:end), y_(1:end), 'MarkerEdgeColor','k','MarkerFaceColor','r');
        else
            scatter(x_(1:end), y_(1:end));
        end
        
        print(h1,'-dpdf' ,       ['generation' num2str(generation)  '/extendedConvexHull.pdf']);
        hold off
    catch
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeCaption(atomType, numIons)
if size(numIons,2)==2
    name = ['Composition ratio: '];
    name = [name megaDoof(atomType(2)) '/(' megaDoof(atomType(1)) '+' megaDoof(atomType(2)) ')'];
    xlabel(name);
    ylabel('Enthalpy of formation (eV/atom)');
else
    name = '';
    for type = 1:size(numIons,2)
        num = numIons(1, type);
        if num > 1
            name = [name megaDoof(atomType(type)) num2str(num)];
        elseif num>0
            name = [name megaDoof(atomType(type))];
        end
    end
    name = [name '+'];
    for type = 1:size(numIons,2)
        num = numIons(2, type);
        if num > 1
            name = [name megaDoof(atomType(type)) num2str(num)];
        elseif num>0
            name = [name megaDoof(atomType(type))];
        end
    end
    xlabel(name);
    ylabel('Enthalpy of formation (eV/block)');
end

