function WriteClusters(fitness)
global POP_STRUC
global ORG_STRUC
global CLUSTERS

resFolder = POP_STRUC.resFolder;
atomType= ORG_STRUC.atomType;
numions = ORG_STRUC.numIons;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------------------------------------------------------------------------------------------------%%%-------------------------------------case of two-component (binary) clusters ------------------------------%%
%%----------------------------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(numions,2) == 2	

numcomp = CLUSTERS.number_ConvexHall;
nx = numions(2,1) - numions(1,1) + 1;
ny = numions(2,2) - numions(1,2) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% printing files best_clusters.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath  = [resFolder '/best_clusters_gen' num2str(POP_STRUC.generation) '.txt'];
fid  = fopen(fpath,  'wt');

%printing energies of the best clusters
fprintf(fid,'        |');
for i = numions(1,2) : numions(2,2)
    fprintf(fid,'    %2.0f     ',i);
end 
fprintf(fid,'\n');
for i = numions(1,2) : numions(2,2) + 1
    fprintf(fid,'----------');
end
fprintf(fid,'\n');

for i = numions(1,1) : numions(2,1)
    fprintf(fid,'  %2.0f    |',i);
    for j = numions(1,2) : numions(2,2)
        ii = i-numions(1,1)+2;
        jj = j-numions(1,2)+2;
        fprintf(fid,' %9.4f ',CLUSTERS.composition(ii,jj).bestEnthalpy);
    end
    fprintf(fid,'\n');
end

%printing magic clusters
fprintf(fid,'\n');
fprintf(fid,'    |');
for i = numions(1,2) : numions(2,2)
    fprintf(fid,' %2.0f ',i);
end 
fprintf(fid,'\n');
for i = numions(1,2) : numions(2,2) + 1
    fprintf(fid,'----');
end
fprintf(fid,'\n');

for i = numions(1,1) : numions(2,1)
  fprintf(fid,' %2.0f |',i);
  for j = numions(1,2) : numions(2,2)
    ii = i-numions(1,1)+2;
    jj = j-numions(1,2)+2;         
    d2Ex = CLUSTERS.composition(ii+1,jj).bestEnthalpy+CLUSTERS.composition(ii-1,jj).bestEnthalpy-2*CLUSTERS.composition(ii,jj).bestEnthalpy;
    d2Ey = CLUSTERS.composition(ii,jj+1).bestEnthalpy+CLUSTERS.composition(ii,jj-1).bestEnthalpy-2*CLUSTERS.composition(ii,jj).bestEnthalpy;
    if (d2Ex >= 0) && (d2Ey >= 0)
       fprintf(fid,'  m ');
    else
       fprintf(fid,'    ');
    end
  end
  fprintf(fid,'\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%% printing files best_clusters_POSCARS.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = numions(1,1) : numions(2,1)
for j = numions(1,2) : numions(2,2)
    ii = i-numions(1,1)+2;
    jj = j-numions(1,2)+2;  
    coor     = CLUSTERS.composition(ii,jj).COORDINATES;
    numIons1 = CLUSTERS.composition(ii,jj).numIons;
    lattice  = CLUSTERS.composition(ii,jj).LATTICE;
    symg     = CLUSTERS.composition(ii,jj).symg;
    order    = CLUSTERS.composition(ii,jj).order;
    count    = CLUSTERS.composition(ii,jj).Number;
    Write_POSCAR(atomType, count, symg, numIons1, lattice, coor);
    unix([' cat POSCAR       >>' resFolder '/best_clusters_gen' num2str(POP_STRUC.generation) '_POSCARS']);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PRINTING GRAPHS NEEDED FOR DEBUGGING, FOLDER "DEBUGGING_FILES" %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DebuggingFolder = resFolder; % [resFolder '/DEBUGGING_FILES'];
%if ~exist(DebuggingFolder)
%    mkdir(DebuggingFolder);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%% reading data from the file "Reference_Data.txt" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_ref = fopen([ORG_STRUC.homePath '/Reference_Data.txt'],'rt');

if (fid_ref~=-1)
N1_ref = fscanf(fid_ref,'%f',1);
N2_ref = fscanf(fid_ref,'%f',1);
ss = fgets(fid_ref);

for i=1:N1_ref
    for j=1:N2_ref
        E_LJ(i,j) = fscanf(fid_ref,'%f',1);
    end
    ss = fgets(fid_ref);
end
fclose(fid_ref);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% plotting 3d graph of energies and reference data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
h = figure;
hold;
set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
ylabel('Energy');

%%%%%%%%%%%%%%%%%% Calculation relative energy of each composition %%%%%%%%%%%%%%%%%%%
x1 = numions(1,1);
y1 = numions(1,2);
z1 = CLUSTERS.composition(2,2).bestEnthalpy;
x2 = numions(2,1);
y2 = numions(1,2);
z2 = CLUSTERS.composition(x2-x1+2,y2-y1+2).bestEnthalpy;
x3 = numions(1,1);
y3 = numions(2,2);
z3 = CLUSTERS.composition(2,y3-y1+2).bestEnthalpy;

E_rel=[];
E_rel1=[];
n=0;
for i = numions(1,1) : numions(2,1)
for j = numions(1,2) : numions(2,2)
    ii = i-numions(1,1)+2;
    jj = j-numions(1,2)+2;
    E_rel(ii,jj)=(det([x1 y1 z1;x2 y2 z2;x3 y3 z3]) - ii*det([y1 z1 1;y2 z2 1;y3 z3 1]) + jj*det([x1 z1 1;x2 z2 1;x3 z3 1]))/det([x1 y1 1;x2 y2 1;x3 y3 1]);
    n=n+1;
    E_rel1(n)=E_rel(ii,jj);
end
end

delta1=0;
for i = numions(1,1) : numions(2,1)
for j = numions(1,2) : numions(2,2)
    ii = i-numions(1,1)+2;
    jj = j-numions(1,2)+2;
    delta1 = delta1 + abs(CLUSTERS.composition(ii,jj).bestEnthalpy - E_rel(ii,jj));
end
end

deltaE = delta1/(nx*ny*(ny-1));
%deltaE=0;

%%%%%%%%%%%%%%%%%%%%%%%%% plotting grid of best structures at current generation %%%%%%%%%%%%%%%%%%%%%
for i = numions(1,1) : numions(2,1)
    xbest=[]; ybest=[];
    nn=0;
    for j= numions(1,2) : numions(2,2)
        ii = i-numions(1,1)+2;
        jj = j-numions(1,2)+2;
        nn = nn + 1;
        xbest(nn) = ii + numions(1,1) - 2 + (1/ny) * (jj-2);
        ybest(nn) = CLUSTERS.composition(ii,jj).bestEnthalpy  - E_rel(ii,jj) - (jj-2)*deltaE;
    end;
    plot(xbest,ybest,':','LineWidth',2);
end;

for j= numions(1,2) : numions(2,2)
    xbest=[]; ybest=[];
    nn=0;
    for i = numions(1,1) : numions(2,1)
        ii = i-numions(1,1)+2;
        jj = j-numions(1,2)+2;
        nn = nn + 1;
        xbest(nn) = ii + numions(1,1) - 2 + (1/ny) * (jj-2);
        ybest(nn) = CLUSTERS.composition(ii,jj).bestEnthalpy  - E_rel(ii,jj) - (jj-2)*deltaE;
    end;
    plot(xbest,ybest,':','LineWidth',2);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting axes, ticks and labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1==0
set(gca,'XTick', 0);
set(gca,'YTick', round(min(E_rel1)+1):round(max(E_rel1)+3));
text(numions(2,1)+0.5,min(E_rel1)-0.2,'numions (A)');
plot([numions(1,1), numions(1,1)+1], [min(E_rel1), min(E_rel1) - (ny+1)*deltaE], 'k', 'LineWidth', 1);
plot([numions(1,1), numions(2,1)+1], [min(E_rel1), min(E_rel1)], 'k', 'LineWidth', 1);
for i = numions(1,2)+1: numions(2,2)
    ttx = numions(1,1)+(i-numions(1,2))/ny;
    tty = min(E_rel1)-(i-numions(1,2))*deltaE*(ny+1)/ny;
    text(ttx-0.02,tty - 0.15, num2str(i));
    plot([ttx,ttx] , [tty,tty+0.1],'k', 'LineWidth', 1);
end
for i = numions(1,1)+1: numions(2,1)
    text(i-0.02,min(E_rel1)-0.15, num2str(i));
    plot([i,i], [min(E_rel1),min(E_rel1)+0.1], 'k', 'LineWidth', 1);
end
text(numions(1,1),min(E_rel1)-0.2, [num2str(numions(1,1)) ',' num2str(numions(1,2))]);
text(numions(1,1)+1, min(E_rel1) - (ny+1)*deltaE - 0.1, 'numions (B)');
text(numions(1,1)+nx/2-1, max(E_rel1)+2.6, ['AxBy      GENERATION ' num2str(POP_STRUC.generation)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=90*3/sqrt(nx*ny);
ccc = 'wrggymck';
operators = ['Random'; 'Heredi'; 'softmu'; 'CoorMu'; 'Permut'; 'TransM'; 'keptBe'; 'convex'];
for k = 1:size(operators,1);
   x=[]; E=[]; n=0;
   for i=1:length(POP_STRUC.POPULATION)
       tmp = [POP_STRUC.POPULATION(i).howCome '1'];
       if ~isempty(findstr(tmp, operators(k,:)))
   	  ix = POP_STRUC.POPULATION(i).numIons(1);
	  iy = POP_STRUC.POPULATION(i).numIons(2);
	  n = n + 1;
	  x(n) = ix + (1/ny) * (iy-numions(1,2));
	  E(n) = POP_STRUC.POPULATION(i).Enthalpies(end) - E_rel(ix-numions(1,1)+2,iy-numions(1,2)+2) - (iy-numions(1,2))*deltaE;
       end 
   end
   scatter(x,E,s,ccc(k),'filled','MarkerEdgeColor','k');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting reference data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fid_ref~=-1)
   xref=[]; Eref=[]; n=0;
   for i = numions(1,1) : numions(2,1)
   for j = numions(1,2) : numions(2,2)
      if E_LJ(i,j)~=0
         n=n+1;
         xref(n) = i + (1/ny) * (j-numions(1,2));
 	 Eref(n) = E_LJ(i,j)  - E_rel(i-numions(1,1)+2,j-numions(1,2)+2) - (j-numions(1,2))*deltaE;
      end
   end
   end
   scatter(xref,Eref,s/3,'k','filled','MarkerEdgeColor','k');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% diapason of axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_ax = numions(1,1);
max_ax = numions(2,1)+1;
nn=0; Ebest = [];
for i = numions(1,1) : numions(2,1)
for j = numions(1,2) : numions(2,2)
    ii = i-numions(1,1)+2;
    jj = j-numions(1,2)+2;
    nn = nn + 1;
    Ebest(nn) = CLUSTERS.composition(ii,jj).bestEnthalpy  - E_rel(ii,jj) - (jj-2)*deltaE;
end
end
min_ay = min(Ebest)-(max(Ebest)-min(Ebest))/30;
max_ay = max(Ebest) + (max(Ebest)-min(Ebest))/2;
if (fid_ref~=-1) && (min(Eref)<min(Ebest))
   min_ay = min(Eref)-(max(Ebest)-min(Eref))/30;
   max_ay = max(Ebest) + (max(Ebest)-min(Eref))/2;
end
axis([min_ax,max_ax,min_ay,max_ay]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variation operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
operators_text = ['Random       '; 'Heredity     '; 'Softmutation '; 'Permutation  '; 'Transmutation'; 'keptBest     '; 'ConvexHull   '];
xxx = max_ax - (max_ax-min_ax)/10;
ccc = 'wrgymck';
for i = 1:size(operators_text,1);
   yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
   scatter(xxx,yyy,s,ccc(i),'filled','MarkerEdgeColor','k');   
   text(xxx+(max_ax-min_ax)/50, yyy, operators_text(i,:));
end
if (fid_ref~=-1)
   scatter(xxx, yyy-(max_ay-min_ay)/30,s/3,'k','filled','MarkerEdgeColor','k');
   text(xxx+(max_ax-min_ax)/50, yyy-(max_ay-min_ay)/30, 'Ref. data');
end

print(h,'-dtiff','-r120',[DebuggingFolder '/Energies_3d_gen' num2str(POP_STRUC.generation) '.tif']);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fitness vs numions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
h = figure;
hold;
set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
xlabel('numIons');
ylabel('Fitness');

%%%%%%%%% Line coonnecting fitnesses of best clusters of all compositions %%%%%%%%%%%%%
xx = []; ff = [];
n=0;
for cx = numions(1,1):numions(2,1)
    for cy = numions(1,2):numions(2,2)
	done = 0;
	for i = 1:numcomp
	    x_ch = CLUSTERS.ConvexHall(i).numIons(1);
	    y_ch = CLUSTERS.ConvexHall(i).numIons(2);
	    if (cx==x_ch) && (cy==y_ch)
		done = 1;
		Eref = CLUSTERS.ConvexHall(i).Enthalpy;
	    end
	end
	
	if ~done
	    for i=1:numcomp
		xch = CLUSTERS.ConvexHall(i).numIons(1);
		ych = CLUSTERS.ConvexHall(i).numIons(2);
		r(i) = sqrt((cx-xch)^2 + (cy-ych)^2);
	    end
       	    [r_upor, n_upor] = sort(r);
	    n_current = 2;
	    while (~done)
		n_current = n_current + 1;
		for i = 1:n_current-2
		    for j = i+1:n_current-1
			x1 = CLUSTERS.ConvexHall(n_upor(i)).numIons(1);
	    		y1 = CLUSTERS.ConvexHall(n_upor(i)).numIons(2);
	    		x2 = CLUSTERS.ConvexHall(n_upor(j)).numIons(1);
	    		y2 = CLUSTERS.ConvexHall(n_upor(j)).numIons(2);
			x3 = CLUSTERS.ConvexHall(n_upor(n_current)).numIons(1);
	    		y3 = CLUSTERS.ConvexHall(n_upor(n_current)).numIons(2);
	    		S = AreaTriangle(x1,y1,x2,y2,x3,y3);
	    		S1 = AreaTriangle(cx,cy,x1,y1,x2,y2);
	    		S2 = AreaTriangle(cx,cy,x1,y1,x3,y3);
	    		S3 = AreaTriangle(cx,cy,x2,y2,x3,y3);
			if (S~=0) && (S1+S2+S3 <= S+0.00000001)
			    E1 = CLUSTERS.ConvexHall(n_upor(i)).Enthalpy;
			    E2 = CLUSTERS.ConvexHall(n_upor(j)).Enthalpy;
			    E3 = CLUSTERS.ConvexHall(n_upor(n_current)).Enthalpy;
			    Eref=(det([x1 y1 E1;x2 y2 E2;x3 y3 E3])-cx*det([y1 E1 1;y2 E2 1;y3 E3 1])+cy*det([x1 E1 1;x2 E2 1;x3 E3 1]))/det([x1 y1 1;x2 y2 1;x3 y3 1]);
			    done = 1;
			end
		    if done break; end
		    end
		if done break; end
		end
	    end

	end
	n=n+1;
	xx(n) = cx + (1/ny) * (cy-numions(1,2));
        ff(n) = CLUSTERS.composition(cx-numions(1,1)+2,cy-numions(1,2)+2).bestEnthalpy - Eref;
    end
end
plot(xx,ff,'--','LineWidth',1);

%%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%
s=90*3/sqrt(nx*ny);
ccc = 'wrggymck';
operators = ['Random'; 'Heredi'; 'softmu'; 'CoorMu'; 'Permut'; 'TransM'; 'keptBe'; 'convex'];
for k = 1:size(operators,1);
   x=[]; F=[]; n=0;
   for i=1:length(POP_STRUC.POPULATION)
       tmp = [POP_STRUC.POPULATION(i).howCome '1'];
       if ~isempty(findstr(tmp, operators(k,:)))
   	  ix = POP_STRUC.POPULATION(i).numIons(1);
	  iy = POP_STRUC.POPULATION(i).numIons(2);
	  n = n + 1;
	  x(n) = ix + (1/ny) * (iy-numions(1,2));
	  F(n) = fitness(i);
        end
   end
   scatter(x,F,s,ccc(k),'filled','MarkerEdgeColor','k');
end

%%%%%%%%%%%%%%%%%%%%%% variation operators %%%%%%%%%%%%%%%%%%%%%%%
min_ax = numions(1,1);
max_ax = numions(2,1)+1;
min_ay = 0;
max_ay = 10;
operators_text = ['Random       '; 'Heredity     '; 'Softmutation '; 'Permutation  '; 'Transmutation'; 'keptBest     '; 'ConvexHull   '];
xxx = max_ax - (max_ax-min_ax)/10;
ccc = 'wrgymck';
for i = 1:size(operators_text,1);
   yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
   scatter(xxx,yyy,60,ccc(i),'filled','MarkerEdgeColor','k');   
   text(xxx+(max_ax-min_ax)/50, yyy, operators_text(i,:));
end
axis([min_ax,max_ax,min_ay,max_ay]);

%%%%%%%%%%%%%%%%%%%%%% ticks & labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0;
max_ticks = 20;
N1_ticks = fix(max_ticks/nx);
delta_tick = fix(ny/N1_ticks) +1;
if (delta_tick==0) delta_tick = 1; end
for ix = numions(1,1):numions(2,1)
for iy = numions(1,2):numions(2,2)
    n=n+1;
    x_tick(n) = ix + (1/ny) * (iy-numions(1,2));
    s_tick(n,:) = '     ';
    if (fix((iy-numions(1,2))/delta_tick)*delta_tick == iy-numions(1,2)) && (numions(2,2)+1-iy >= delta_tick)
        str_tick = [num2str(ix) ',' num2str(iy)];
        if length(num2str(ix))==1 str_tick = [' ' str_tick]; end
        if length(num2str(iy))==1 str_tick = [str_tick ' ']; end
        s_tick(n,:) = str_tick;
    end
end
end

set(gca,'XTick', x_tick,'XTickLabel',s_tick);
print(h,'-dtiff','-r120',[DebuggingFolder '/Fitness_gen' num2str(POP_STRUC.generation) '.tif']);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy/numions vs numions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
h = figure;
hold;
set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
xlabel('numIons');
ylabel('Energy / numIons');

%%%%%%%% Plotting line connecting best clusters of all compositions in current generations %%%%%%%%%%
Ebest = [];
Xbest = [];
n=0;
for ix = numions(1,1):numions(2,1)
for iy = numions(1,2):numions(2,2)
    n=n+1;
    Ebest(n) = CLUSTERS.composition(ix-numions(1,1)+2,iy-numions(1,2)+2).bestEnthalpy / (ix+iy);
    Xbest(n) = ix + (1/ny) * (iy-numions(1,2));
end
end
plot(Xbest,Ebest,'--','LineWidth',1);

%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%
s=90*3/sqrt(nx*ny);
ccc = 'wrggymck';
operators = ['Random'; 'Heredi'; 'softmu'; 'CoorMu'; 'Permut'; 'TransM'; 'keptBe'; 'convex'];
for k = 1:size(operators,1);
   x=[]; EdivN=[]; n=0;
   for i=1:length(POP_STRUC.POPULATION)
       tmp = [POP_STRUC.POPULATION(i).howCome '1'];
       if ~isempty(findstr(tmp, operators(k,:)))
   	  ix = POP_STRUC.POPULATION(i).numIons(1);
	  iy = POP_STRUC.POPULATION(i).numIons(2);
	  n = n + 1;
	  x(n) = ix + (1/ny) * (iy-numions(1,2));
	  EdivN(n) = POP_STRUC.POPULATION(i).Enthalpies(end) / (ix+iy);
        end
   end
   scatter(x,EdivN,s,ccc(k),'filled','MarkerEdgeColor','k');
end

%%%%%%%%%%%%%%%%%%%%%% Plotting reference data %%%%%%%%%%%%%%%%%%%
if (fid_ref~=-1)
   Eref = [];
   Xref = [];
   n=0;
   for ix = numions(1,1):numions(2,1)
   for iy = numions(1,2):numions(2,2)
      if (E_LJ(ix,iy)~=0)
          n=n+1;
          Eref(n) = E_LJ(ix,iy) / (ix+iy);
          Xref(n) = ix + (1/ny) * (iy-numions(1,2));
  	  plot([Xref(n)-1/(2*ny), Xref(n)+1/(2*ny)], [Eref(n),Eref(n)],'k','LineWidth',1);
       end
   end
   end
end

%%%%%%%%%%%%%%%%%%%%%% variation operators %%%%%%%%%%%%%%%%%%%%%%%
min_ax = numions(1,1);
max_ax = numions(2,1)+1;
min_ay = min(Ebest)-(max(Ebest)-min(Ebest))/30;
max_ay = max(Ebest) + (max(Ebest)-min(Ebest))/2;
if (fid_ref~=-1) && (min(Eref)<min(Ebest))
   min_ay = min(Eref)-(max(Ebest)-min(Eref))/30;
   max_ay = max(Ebest) + (max(Ebest)-min(Eref))/2;
end
operators_text = ['Random       '; 'Heredity     '; 'Softmutation '; 'Permutation  '; 'Transmutation'; 'keptBest     '; 'ConvexHull   '];
xxx = max_ax - (max_ax-min_ax)/10;
ccc = 'wrgymck';
for i = 1:size(operators_text,1);
   yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
   scatter(xxx,yyy,60,ccc(i),'filled','MarkerEdgeColor','k');   
   text(xxx+(max_ax-min_ax)/50, yyy, operators_text(i,:));
end
if (fid_ref~=-1)
   plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
   text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
end
axis([min_ax,max_ax,min_ay,max_ay]);

%%%%%%%%%%%%%%%%%%%%%% ticks & labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0;
max_ticks = 20;
N1_ticks = fix(max_ticks/nx);
delta_tick = fix(ny/N1_ticks) +1;
if (delta_tick==0) delta_tick = 1; end
for ix = numions(1,1):numions(2,1)
for iy = numions(1,2):numions(2,2)
    n=n+1;
    x_tick(n) = ix + (1/ny) * (iy-numions(1,2));
    s_tick(n,:) = '     ';
    if (fix((iy-numions(1,2))/delta_tick)*delta_tick == iy-numions(1,2)) && (numions(2,2)+1-iy >= delta_tick)
        str_tick = [num2str(ix) ',' num2str(iy)];
        if length(num2str(ix))==1 str_tick = [' ' str_tick]; end
        if length(num2str(iy))==1 str_tick = [str_tick ' ']; end
        s_tick(n,:) = str_tick;
    end
end
end

set(gca,'XTick', x_tick,'XTickLabel',s_tick);
print(h,'-dtiff','-r120',[DebuggingFolder '/Energy_div_numIons_gen' num2str(POP_STRUC.generation) '.tif']);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------------------------------------------------------------------------------------------------%%
%%------------------------------------ case of one-component clusters -------------------------------------%%
%%---------------------------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(numions,2) == 1	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% printing files best_clusters.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath  = [resFolder '/best_clusters' num2str(POP_STRUC.generation) '.txt'];
fid  = fopen(fpath,  'wt');
fprintf(fid,'composition  energy - all compositions\n');
for i = 1:CLUSTERS.number_compositions
    fprintf(fid,'%9.0f  %9.3f\n',CLUSTERS.composition(i).numIons, CLUSTERS.composition(i).bestEnthalpy);
end

fprintf(fid,'\n');
fprintf(fid,'composition  energy - convex hull\n');
for i = 1:CLUSTERS.ConvexHall_numComp
    fprintf(fid,'%9.0f  %9.3f\n',CLUSTERS.ConvexHall(i).numIons, CLUSTERS.ConvexHall(i).Enthalpy);
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%% printing files best_clusters_POSCARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = numions(1,1) : numions(2,1)
    ii = i-numions(1,1)+1;
    coor     = CLUSTERS.composition(ii).COORDINATES;
    numIons1 = CLUSTERS.composition(ii).numIons;
    lattice  = CLUSTERS.composition(ii).LATTICE;
    symg     = CLUSTERS.composition(ii).symg;
    order    = CLUSTERS.composition(ii).order;
    count    = CLUSTERS.composition(ii).Number;
    Write_POSCAR(atomType, count, symg, numIons1, lattice, coor);
    unix([' cat POSCAR       >>' resFolder '/best_clusters' num2str(POP_STRUC.generation) '_POSCARS']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% reading data from the file "Reference_Data.txt" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_ref = fopen([ORG_STRUC.homePath '/Reference_Data.txt'],'rt');
if (fid_ref~=-1)
   ss=fgets(fid);
   for i=2:75
      N_atoms = fscanf(fid,'%f',1);
      E_RefData(N_atoms) = fscanf(fid,'%f',1);
      ss=fgets(fid);
   end;
   fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% relative energies vs numIons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
h = figure;
hold;
set(gcf,'Visible','off');
xlabel('numIons');
ylabel('Relative energy');

N1 = numions(1,1);
N2 = numions(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting reference data %%%%%%%%%%%%%%%%%%%%%%%%%%
if (fid_ref~=-1)
  x_ref = []; E_ref = [];
  for i = 1:N1-1
    E_ref(i)=10000;
  end
  for i = N1:N2
    E_ref(i) = E_RefData(i) - (E_RefData(N1)+(i-N1)*(E_RefData(N2)-E_RefData(N1))/(N2-N1));
    plot([i-1/4, i+1/4], [E_ref(i),E_ref(i)],'k','LineWidth',1);
  end
end

%%%%%%%% plotting line connecting the best clusters of all compositions %%%%%%%%%%
  x_best = []; E_best = [];
  for i = 1:N1-1
    E_best(i)=10000;
  end
  if (fid_ref~=-1)
    E_Ref1 = E_RefData(N1);
    E_Ref2 = E_RefData(N2);
  else 
    E_Ref1 = CLUSTERS.composition(1).bestEnthalpy;
    E_Ref2 = CLUSTERS.composition(N2-N1+1).bestEnthalpy;
  end
  for i = N1:N2
    x_best(i) = i;
    E_best(i) = CLUSTERS.composition(i-N1+1).bestEnthalpy - (E_Ref1+(i-N1)*(E_Ref2-E_Ref1)/(N2-N1));
  end
  plot(x_best,E_best,'--','LineWidth',1);

%%%%%%%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=600/CLUSTERS.number_compositions;
ccc = 'wrggymck';
operators = ['Random'; 'Heredi'; 'softmu'; 'CoorMu'; 'Permut'; 'TransM'; 'keptBe'; 'convex'];
for k = 1:size(operators,1);
   x=[]; E=[]; n=0;
   for i=1:length(POP_STRUC.POPULATION)
       tmp = [POP_STRUC.POPULATION(i).howCome '1'];
       if ~isempty(findstr(tmp, operators(k,:)))
	  n = n + 1;
	  x(n) = POP_STRUC.POPULATION(i).numIons;
	  E(n) = POP_STRUC.POPULATION(i).Enthalpies(end) - (E_Ref1+(x(n)-N1)*(E_Ref2-E_Ref1)/(N2-N1));
       end
   end
   scatter(x,E,s,ccc(k),'filled','MarkerEdgeColor','k');
end

%%%%%%%%%%%%%%%%%%%%%% variation operators, diapason of axis %%%%%%%%%%%%%%%%%%%%%%%
E = [];
for i=1:length(POP_STRUC.POPULATION)
  E(i) = POP_STRUC.POPULATION(i).Enthalpies(end) - (E_Ref1+(POP_STRUC.POPULATION(i).numIons-N1)*(E_Ref2-E_Ref1)/(N2-N1));
end
min_ax = numions(1,1);
max_ax = numions(2,1);
min_ay = min(E_best)-std(E)/20;
max_ay = min(E_best)+2*std(E);
if (fid_ref~=-1)
  min_ay = min(min(E_ref),min(E_best))-std(E)/20;
  max_ay = min(min(E_ref),min(E_best))+2*std(E);
end

operators_text = ['Random       '; 'Heredity     '; 'Softmutation '; 'Permutation  '; 'Transmutation'; 'keptBest     '; 'ConvexHull   '];
xxx = max_ax - (max_ax-min_ax)/10;
ccc = 'wrgymck';
for i = 1:size(operators_text,1);
   yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
   scatter(xxx,yyy,60,ccc(i),'filled','MarkerEdgeColor','k');   
   text(xxx+(max_ax-min_ax)/50, yyy, operators_text(i,:));
end
if (fid_ref~=-1)
   plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
   text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
end

axis([min_ax,max_ax,min_ay,max_ay]);
set(gca,'XTick', N1:N2)
print(h,'-dtiff','-r120',[resFolder '/Relative_energies_gen' num2str(POP_STRUC.generation) '.tif']);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fittnes vs numIons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
h = figure;
hold;
set(gcf,'Visible','off');
xlabel('numIons');
ylabel('Fitness');

s=600/CLUSTERS.number_compositions;
ccc = 'wrggymck';
operators = ['Random'; 'Heredi'; 'softmu'; 'CoorMu'; 'Permut'; 'TransM'; 'keptBe'; 'convex'];
for k = 1:size(operators,1);
   x=[]; F=[]; n=0;
   for i=1:length(POP_STRUC.POPULATION)
       tmp = [POP_STRUC.POPULATION(i).howCome '1'];
       if ~isempty(findstr(tmp, operators(k,:)))
	  n = n + 1;
	  x(n) = POP_STRUC.POPULATION(i).numIons;
	  F(n) = fitness(i);
        end
   end
   scatter(x,F,s,ccc(k),'filled','MarkerEdgeColor','k');
end

%%%%%%%%%%%%%%%%%%%%%% Fitnesses of best clusters in current generation %%%%%%%%%%%%
xx=[]; FF=[];
for i = N1:N2
    xx(i)=i;
    ii=1;
    while 1
       ii=ii+1;
       if (CLUSTERS.ConvexHall(ii).numIons > i)
          break;
       end
    end
    nn1 = CLUSTERS.ConvexHall(ii-1).numIons;
    nn2 = CLUSTERS.ConvexHall(ii).numIons;
    EE1 = CLUSTERS.ConvexHall(ii-1).Enthalpy;
    EE2 = CLUSTERS.ConvexHall(ii).Enthalpy;
    FF(i) = (CLUSTERS.composition(i-N1+1).bestEnthalpy - (EE1+(i-nn1)*(EE2-EE1)/(nn2-nn1))) / i;
end
plot(xx,FF,'--','LineWidth',1);

%%%%%%%%%%%%%%%%%%%%%% variation operators, diapason of axis %%%%%%%%%%%%%%%%%%%%%%%
F = [];
for i=1:length(POP_STRUC.POPULATION)
  F(i) = fitness(i);
end
min_ax = numions(1,1);
max_ax = numions(2,1);
min_ay = 0;
max_ay = 2*std(F);
operators_text = ['Random       '; 'Heredity     '; 'Softmutation '; 'Permutation  '; 'Transmutation'; 'keptBest     '; 'ConvexHull   '];
xxx = max_ax - (max_ax-min_ax)/10;
ccc = 'wrgymck';
for i = 1:size(operators_text,1);
   yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
   scatter(xxx,yyy,60,ccc(i),'filled','MarkerEdgeColor','k');   
   text(xxx+(max_ax-min_ax)/50, yyy, operators_text(i,:));
end
if (fid_ref~=-1)
   plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
   text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
end

axis([min_ax,max_ax,min_ay,max_ay]);
set(gca,'XTick', N1:N2)
print(h,'-dtiff','-r120',[resFolder '/Fitness_gen' num2str(POP_STRUC.generation) '.tif']);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy/numions vs numions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
h = figure;
hold;
set(gcf,'Visible','off');
xlabel('numIons');
ylabel('Energy / numIons');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting reference data %%%%%%%%%%%%%%%%%%%%%%%%%%
if (fid_ref~=-1)
  E_ref = [];
  for i = 1:N1-1
    E_ref(i)=10000;
  end
  for i = N1:N2
    E_ref(i) = E_RefData(i) / i;
    plot([i-1/4, i+1/4], [E_ref(i),E_ref(i)],'k','LineWidth',1);
  end
end

%%%%%%%% plotting line connecting the best clusters of all compositions %%%%%%%%%%
  x_best = []; E_best = []; 
  for i = N1:N2
    x_best(i) = i;
    E_best(i) = CLUSTERS.composition(i-N1+1).bestEnthalpy / i;
  end
  for i = 1:N1-1
    E_best(i) = E_best(N1);
  end
  plot(x_best,E_best,'--','LineWidth',1);

%%%%%%%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=600/CLUSTERS.number_compositions;
ccc = 'wrggymck';
operators = ['Random'; 'Heredi'; 'softmu'; 'CoorMu'; 'Permut'; 'TransM'; 'keptBe'; 'convex'];
for k = 1:size(operators,1);
   x=[]; E=[]; n=0;
   for i=1:length(POP_STRUC.POPULATION)
       tmp = [POP_STRUC.POPULATION(i).howCome '1'];
       if ~isempty(findstr(tmp, operators(k,:)))
	  n = n + 1;
	  x(n) = POP_STRUC.POPULATION(i).numIons;
	  E(n) = POP_STRUC.POPULATION(i).Enthalpies(end) / x(n);
       end
   end
   scatter(x,E,s,ccc(k),'filled','MarkerEdgeColor','k');
end

%%%%%%%%%%%%%%%%%%%%%% variation operators, diapason of axis %%%%%%%%%%%%%%%%%%%%%%%
E = [];
for i=1:length(POP_STRUC.POPULATION)
  E(i) = POP_STRUC.POPULATION(i).Enthalpies(end) / POP_STRUC.POPULATION(i).numIons;
end
min_ax = numions(1,1);
max_ax = numions(2,1);
min_ay = min(E_best) - std(E)/10;
max_ay = max(E_best) + (max(E_best)-min(E_best));
if (fid_ref~=-1)
  min_ay = min(min(E_ref),min(E_best))-(max(E_best)-min(E_best))/10;
end

operators_text = ['Random       '; 'Heredity     '; 'Softmutation '; 'Permutation  '; 'Transmutation'; 'keptBest     '; 'ConvexHull   '];
xxx = max_ax - (max_ax-min_ax)/10;
ccc = 'wrgymck';
for i = 1:size(operators_text,1);
   yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
   scatter(xxx,yyy,60,ccc(i),'filled','MarkerEdgeColor','k');   
   text(xxx+(max_ax-min_ax)/50, yyy, operators_text(i,:));
end
if (fid_ref~=-1)
   plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
   text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
end

axis([min_ax,max_ax,min_ay,max_ay]);
set(gca,'XTick', N1:N2)
print(h,'-dtiff','-r120',[resFolder '/Energy_div_numIons_gen' num2str(POP_STRUC.generation) '.tif']);
catch
end

end
