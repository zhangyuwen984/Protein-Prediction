function read_SUBSTRATE()
global ORG_STRUC
    [fid,message] = fopen('POSCAR_SUBSTRATE');
     tmp = fgetl(fid) ;% system description
     factor = str2num(fgetl(fid)); % 1 = numbers in angstroms
%lattice
     LATT = zeros(3);
     for i = 1 : 3
       tmp = fgetl(fid);
       LATT(i,:) = str2num(tmp);
     end
     ORG_STRUC.bulk_lat = LATT*factor;
%%ntyp
     tmp = fgetl(fid);
     ntyp = str2num(tmp);
     if isempty(ntyp)
        tmp=fgetl(fid);
        ntyp = str2num(tmp);
     end
     ORG_STRUC.bulk_ntyp = ntyp;

     natom = sum(ntyp); % number of atoms
%%%
     tmp = fgetl(fid);
     bulk_coordinates = fscanf(fid,'%g',[3,natom]);
     bulk_coordinates = bulk_coordinates';
     ORG_STRUC.bulk_pos(:,1) = bulk_coordinates(:,1);
     ORG_STRUC.bulk_pos(:,2) = bulk_coordinates(:,2);
     ORG_STRUC.bulk_pos(:,3) = bulk_coordinates(:,3);

fclose(fid);
%%% readjust
     T = (1-max(bulk_coordinates(:,3))+min(bulk_coordinates(:,3)))*ORG_STRUC.bulk_lat(3,3);
     if T > 1
        disp('The POSCAR_SUBSTRATE contains too much empty space, Please reduce it ');
        bulk_coordinates(:,3) = bulk_coordinates(:,3) * ORG_STRUC.bulk_lat(3,3);
        ORG_STRUC.bulk_lat(3,3) = ORG_STRUC.bulk_lat(3,3) - T + 0.5;
        bulk_coordinates(:,3) = bulk_coordinates(:,3) - min(bulk_coordinates(:,3)) + 0.1;
        bulk_coordinates(:,3) = bulk_coordinates(:,3)/ORG_STRUC.bulk_lat(3,3);

        disp('Below is our suggested POSCAR_SUBSTRATE ');
        fp = fopen('POSCAR_SUBSTRATE_NEW', 'w');
        fprintf(fp,'POSCAR_SUBSTRATE_NEW\n');
        fprintf(fp,'1.000000\n');
        fprintf(fp,'%8.3f %8.3f %8.3f\n', ORG_STRUC.bulk_lat(1,:));
        fprintf(fp,'%8.3f %8.3f %8.3f\n', ORG_STRUC.bulk_lat(2,:));
        fprintf(fp,'%8.3f %8.3f %8.3f\n', ORG_STRUC.bulk_lat(3,:));
        for i=1:length(ORG_STRUC.bulk_ntyp)
        fprintf(fp,'%4d', ORG_STRUC.bulk_ntyp(i));
        end
        fprintf(fp,'\n');
        fprintf(fp,'Direct\n');
        for i=1:natom
            fprintf(fp,'%9.5f %9.5f %9.5f\n', bulk_coordinates(i,:));
        end
        fclose(fp);
        [nothing, nothing] = unix(['cat POSCAR_SUBSTRATE_NEW']);

        disp('If you believe this new file is good,');
        disp('Please use it to overwrite the old file and resubmit the calculation');
        disp('------------BYE---------------');
        
     quit
     end
     ind = 1;
     ORG_STRUC.bulk_atyp = zeros(sum(ntyp),1);
     for i=1:length(ntyp)
         for j=1:ntyp(i)
             ORG_STRUC.bulk_atyp(ind) = ORG_STRUC.atomType(i);
             ind = ind + 1;
         end
     end
