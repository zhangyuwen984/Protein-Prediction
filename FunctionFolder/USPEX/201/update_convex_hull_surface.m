function [fitness] = update_convex_hull_surface()
global POP_STRUC
global ORG_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(ORG_STRUC.bulk_ntyp >0) == 1   %substrate is elemental
     if length(ORG_STRUC.numIons) == 1 % only one type of surface atom(eg: C on C(111))
          if POP_STRUC.generation == 1
                 ORG_STRUC.E_AB = POP_STRUC.POPULATION(1).Enthalpies(end);  % for one component
          end
          for i=1:length(POP_STRUC.POPULATION)
                 fitness(i)=POP_STRUC.POPULATION(i).Enthalpies(end)-ORG_STRUC.E_AB*prod(POP_STRUC.POPULATION(i).cell);
              if sum(POP_STRUC.POPULATION(i).Surface_numIons)>0
                 fitness(i)=fitness(i)/sum(POP_STRUC.POPULATION(i).Surface_numIons);
              else
                 fitness(i)=0;
              end
          end
     elseif length(ORG_STRUC.numIons) == 2 % two types of surface (eg: PdO on Pd(100))
          flag = 1;
          E_A = ORG_STRUC.E_A; %-2.524014/2;
          E_B = ORG_STRUC.E_B; %-8.719133/2;
          mu_max = E_B;        %partial pressure 10e-12
          mu_min = E_B - 2;    %partial pressure 10e-12
          POP_STRUC.convex_hull = zeros(2,3);
          POP_STRUC.convex_hull(1,1) = 0;   %upper bound
          POP_STRUC.convex_hull(2,1) = ORG_STRUC.numIons(2);   %lower bound
          POP_STRUC.convex_hull(:,2) = 9999; %initial hull, kill out unreasonable structures.
          for i = 1:length(POP_STRUC.POPULATION)
              N_cell = prod(POP_STRUC.POPULATION(i).cell);
              E_DFT = POP_STRUC.POPULATION(i).Enthalpies(end);
              N_A = POP_STRUC.POPULATION(i).Surface_numIons(1);
              N_B = POP_STRUC.POPULATION(i).Surface_numIons(2);
              E(i) = (E_DFT - N_A*E_A)/N_cell;
              X(i) = N_B/N_cell;
              POP_STRUC.convex_hull = update_hull( POP_STRUC.convex_hull, E(i), X(i), mu_min, mu_max, i);
          end
     end
elseif sum(ORG_STRUC.numIons>0)<3 %substrate is binary system (eg: ZnO)

      flag = 2;
      E_AB = ORG_STRUC.E_AB; %-18.181794/2;
      E_A = ORG_STRUC.E_A; %-2.524014/2;
      E_B = ORG_STRUC.E_B; %-8.719133/2;
      m = ORG_STRUC.bulk_stoi(1); %-8.719133/2;
      n = ORG_STRUC.bulk_stoi(2); %-8.719133/2;
      mu_min = (E_AB - n*E_B)/m ;
      mu_max = E_A;
      POP_STRUC.convex_hull = zeros(2,3);
      POP_STRUC.convex_hull(1,1) =-ORG_STRUC.numIons(2)*m/n;   %upper bound
      POP_STRUC.convex_hull(2,1) = ORG_STRUC.numIons(1);   %lower bound
      POP_STRUC.convex_hull(:,2) = 9999; %initial hull

      for i = 1:length(POP_STRUC.POPULATION)
          cell = POP_STRUC.POPULATION(i).cell;
          E_DFT = POP_STRUC.POPULATION(i).Enthalpies(end);
          N_A = POP_STRUC.POPULATION(i).Surface_numIons(1);
          N_B = POP_STRUC.POPULATION(i).Surface_numIons(2);
          E(i) = (E_DFT - N_B*E_AB/n)/prod(cell);
          X(i) = (N_A-m/n*N_B)/prod(cell);
          POP_STRUC.convex_hull = update_hull( POP_STRUC.convex_hull, E(i), X(i), mu_min, mu_max, i);
      end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(POP_STRUC.convex_hull)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
delta = (mu_max - mu_min)/10;
mu = [mu_min:delta:mu_max];
h = figure;
set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
hold on;
box on;
A = megaDoof(ORG_STRUC.atomType(1));
B = megaDoof(ORG_STRUC.atomType(2));
if flag == 1
title(['E_{formation} = E - N_{' A '}\mu_{' A '} - N_{' B '} \mu_{' B '}']);
xlabel(['\mu(' B ')']);
else
title(['E_{formation} = E - N_{' B '} E_AB/n - (N_{' A '} - m N_{' B '}/n) \mu_{' A '}']);
xlabel(['\mu(' A ')']);
end
ylabel('Enthalpy of formation (eV/cell)');
xlim([mu_min, mu_max]);
count = 0;
for i=1:size(POP_STRUC.convex_hull, 1)
   if POP_STRUC.convex_hull(i,3) > 0
       count = count + 1;
       Ind = round(POP_STRUC.convex_hull(i,3));
       N_cell = prod(POP_STRUC.POPULATION(Ind).cell);
       E_DFT = POP_STRUC.POPULATION(Ind).Enthalpies(end);
       Number = POP_STRUC.POPULATION(Ind).Number;
       N_A = POP_STRUC.POPULATION(Ind).Surface_numIons(1);
       N_B = POP_STRUC.POPULATION(Ind).Surface_numIons(2);
       if flag == 1
          E_form(count,:) = (E_DFT - N_A*E_A - N_B*mu)/N_cell;
       else
          E_form(count,:) = (E_DFT - N_B*E_AB/n - (N_A-m*N_B/n)*mu)/N_cell;
       end
       Fig_title{count} = ['Structure_ ' num2str(Number) ];
   end
end

plot(mu, E_form);
legend(Fig_title);
file_path = [POP_STRUC.resFolder '/Surface_Diagram.pdf'];
print(h,'-dpdf', file_path);
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fitness = toconvexhull(POP_STRUC.convex_hull, E, X);
   disp(['Convex hull updated and fitness calculated ']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  fit = toconvexhull(hull, E, X);

np = size(hull,1);
N = length(E);
for loop = 1:N
   X10 = min(hull(:,1));
   X20 = max(hull(:,1));
   done = 0;
    for i = 1:np
        if abs(X(loop)-hull(i,1))<0.001  %stable point
              fit(loop) = E(loop)- hull(i,2);
              done = 1;
              break
        elseif X(loop) < hull(i,1)
            if hull(i,1) < X20 + 0.001
              X2 = i;
              X20 = hull(i,1);
            end
        elseif X(loop) > hull(i,1)
            if hull(i,1) > X10 - 0.001
              X1 = i;
              X10 = hull(i,1);
            end
        end
    end
    if ~done
       E0= hull(X1,2)+  (X(loop)-X10)*(hull(X1,2)-hull(X2,2))/(X10-X20);
       fit(loop)=E(loop)-E0;
    end
end

