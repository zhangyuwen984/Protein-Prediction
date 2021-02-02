function plot_tenary_diagram(atomType, numIons, Convex_Hull, X, Y, Energy_base)
%plot triangle phase diagram based on convex_hull 
%INPUT: 
%1,atomType [3 5 7] [Li, B, N] 
%2,Convex_hull [composition; energy] [1 0 0 0]
OUT = [X(1,:) Y(1)];
for i=2:length(Y)
   same = 0;
   for j=1:i-1
       if norm(X(i,:)-X(j,:))<0.0001
          same = 1;
          break;
       end
   end
   if same==0    
      tmp = [X(i,1), X(i,2), Y(i)];
      OUT = [OUT;tmp];
   end
end
%1-----------phase diagram
h = figure;
title(['Phase diagram']);
for i = 1:size(Convex_Hull,1)
    E = Convex_Hull(i,end-1);
    comp = Convex_Hull(i,1:3);
    [P(i,1:2), P(i,3)] = Get_XY(numIons, comp, Energy_base, E);
end
K = convhull(P);
trisurf(K,P(:,1),P(:,2),P(:,3),...
       'FaceColor','white','LineWidth',1.5);
position = [-0.9, -0.5; 0, 1.04; 0.9, -0.5];
Label = makeLabel(atomType, numIons,position);
Name = [Label{1} '-' Label{2} '-' Label{3} '-diagram'];
view(180,-90); %view position
axis off
print(h,'-dpdf', '-r300', Name);

%2---------------decomposition map
h = figure;
title([   'Formation Energy (eV/Block) ']);
hold on;
a=[-sqrt(3)/2,-1/2];  % 1 0 0
b=[0 , 1 ];           % 0 1 0
c=[sqrt(3)/2, -1/2];
line([a(1),b(1)],[a(2),b(2)])
line([a(1),c(1)],[a(2),c(2)])
line([c(1),b(1)],[c(2),b(2)])
axis([-1 1 -0.8 1.2])
axis square
axis off
scatter(-1*OUT(:,1),OUT(:,2), 25, OUT(:,3), 'o','filled');
position = [0.9, -0.5; 0, 1.04; -0.9, -0.5];
Label = makeLabel(atomType, numIons,position);
colorbar('Position',[0.850 0.233 0.035 0.631]);
Name = [Label{1} '-' Label{2} '-' Label{3} '-formationE'];
print(h,'-dpdf', '-r300', Name);
