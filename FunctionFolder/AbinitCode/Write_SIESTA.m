function Write_SIESTA(Ind_No)

global POP_STRUC
global ORG_STRUC

Step = POP_STRUC.POPULATION(Ind_No).Step;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;

if ORG_STRUC.molecule==1
   numMols = POP_STRUC.POPULATION(Ind_No).numMols;
   MtypeLIST = POP_STRUC.POPULATION(Ind_No).MtypeLIST;
   typesAList= POP_STRUC.POPULATION(Ind_No).typesAList;
elseif ORG_STRUC.dimension==2
   typesAList= POP_STRUC.POPULATION(Ind_No).typesAList;
else
   item = 0;
   for i = 1:length(numIons)
       for j=1:numIons(i)
           item = item + 1;
           typesAList(item)=ORG_STRUC.atomType(i);
       end
   end
end

try
   [nothing, nothing] = unix(['cat input_' num2str(Step) '.fdf > input.fdf' ]);
catch
   disp(['OOOPPPSSS_input_fdf_is_not_present_for_Step_' num2str(Step) ]);
   quit
end
%%%%%%%%%%Crystal Structure%%%%%%%%%%%%%%
   fid = fopen('STRUC','wt');
   fprintf(fid,'NumberOfAtoms            %4d\n', sum(numIons));
   fprintf(fid,'NumberOfSpecies          %4d\n', length(ORG_STRUC.atomType));
%%%%%%%%LATTCE %%%%%%%%%%%%%%%%%%%%%%%%
   fprintf(fid,'LatticeConstant       1  Ang  \n');
   fprintf(fid,'%%block LatticeVectors\n');
   
   for latticeLoop = 1:3
       fprintf(fid,'%10.6f %10.6f %10.6f\n', lattice(latticeLoop,:));
   end
   
   fprintf(fid,'%%endblock LatticeVectors\n');
%%%%%%%%%%COORDINATES
   if ORG_STRUC.molecule==1
     fprintf(fid,'ZM.UnitsLength Ang\n');
     fprintf(fid,'ZM.UnitsAngle rad\n');
     %write Zmatrix
     fprintf(fid,'%%block Zmatrix\n');
    for ind = 1: sum(numMols)
        format = ORG_STRUC.STDMOL(MtypeLIST(ind)).format;
        ZZZ=POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).ZMATRIX;
        infoFile = cat(2,ORG_STRUC.STDMOL(MtypeLIST(ind)).types, format,ZZZ,ORG_STRUC.STDMOL(MtypeLIST(ind)).optFlags);
        fprintf(fid, 'molecule_cartesian\n');
        fprintf (fid,'%g %g %g %g %f %f %f %g %g %g\n', infoFile');
    end
    fprintf(fid,'%%endblock Zmatrix\n');
   else
    fprintf(fid,'AtomicCoordinatesFormat   Fractional\n');
     fprintf(fid,'%%block AtomicCoordinatesAndAtomicSpecies\n');
     for i=1:sum(numIons)
        for j=1:length(ORG_STRUC.atomType)
            if ORG_STRUC.atomType(j)==typesAList(i)
               label=j;
            end
        end
        fprintf(fid,'%10.6f %10.6f %10.6f %2d\n', COORDINATES(i,:),label);
     end
     fprintf(fid,'%%endblock AtomicCoordinatesAndAtomicSpecies\n');
   end
%KPOINTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     [Kpoints, Error] = Kgrid(lattice, ORG_STRUC.Kresol(Step), ORG_STRUC.dimension);
     if Error == 1  % This LATTICE is extremely wrong, let's skip it from now
         POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
     else
         POP_STRUC.POPULATION(Ind_No).K_POINTS(Step,:)=Kpoints;
     end
     
     fprintf(fid,'%%block kgrid_Monkhorst_Pack\n');
     fprintf(fid,'%2d  0  0 0.0\n', Kpoints(1));
     fprintf(fid,' 0 %2d  0 0.0\n', Kpoints(2));
     fprintf(fid,' 0  0 %2d 0.0\n', Kpoints(3));
     fprintf(fid,'%%endblock kgrid_Monkhorst_Pack\n');
     fclose(fid);
     
     [nothing, nothing] = unix('cat STRUC >> input.fdf');

if ORG_STRUC.dimension==2
     lat_bulk = latConverter(POP_STRUC.POPULATION(Ind_No).Bulk_LATTICE);
     Lattice_par = latConverter(POP_STRUC.POPULATION(Ind_No).LATTICE);

     fid = fopen('Constraint','wt');
     fprintf(fid,'%%block GeometryConstraints\n');
     A1=1;
     A2=-1;
     for type = 1:length(ORG_STRUC.atomType)
         for i=1:numIons(type)
             if POP_STRUC.POPULATION(Ind_No).chanAList(A1) == 0
                if COORDINATES(i,3)*Lattice_par(3) > lat_bulk(3) - ORG_STRUC.thicknessB
                    A2=A2+1;
                end
             end
         end
         if A2 >= 0
            fprintf(fid,'     position  from %4d to %4d \n', A1, A1+A2);
         end
         A1 = A1 + numIons(type);
         A2 = -1;
     end
     fprintf(fid,'%%endblock GeometryConstraints\n');
     fclose(fid)
     [nothing, nothing] = unix('cat Constraint >> input.fdf');
end

