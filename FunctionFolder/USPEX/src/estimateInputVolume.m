function estimateInputVolume()

global ORG_STRUC

tempFolder = [ORG_STRUC.homePath '/CalcFoldTemp'];
codeFolder = [ORG_STRUC.USPEXPath '/FunctionFolder/Tool'];

ORG_STRUC.latVolume(1)=0;

if ORG_STRUC.molecule == 1
    uff = ORG_STRUC.STDMOL;
    for i = 1:length(uff)
        numIons(1:length(uff(i).types))=1;
        for j = 1:length(uff(i).types)
            atomType(j)=ORG_STRUC.atomType(uff(i).types(j));
        end
        atomVolume(i) =  calcDefaultVolume(numIons, atomType, ORG_STRUC.ExternalPressure, 0, tempFolder, codeFolder);		
        moleVolume(i) =  calcDefaultVolume(numIons, atomType, ORG_STRUC.ExternalPressure, 1, tempFolder, codeFolder);
                
        if ORG_STRUC.dimension == 1
            trueVolume(i)=  0.4*moleVolume(i) + 0.6*atomVolume(i);
        else
            trueVolume(i) = (moleVolume(i) + atomVolume(i))/2;
        end
        
		if ORG_STRUC.varcomp == 1
			ORG_STRUC.latVolume(i) = trueVolume(i);
		else
			ORG_STRUC.latVolume = ORG_STRUC.latVolume + ORG_STRUC.numIons(i)*trueVolume(i);
		end
        clear numIons atomType
    end

%% Print Volume for MOLs
    disp(['']);
    disp(['Default Volume for MOLs : ' ]);
    for i = 1:length(uff)
        disp(['     Volume of MOL_' num2str(i) ' : ' num2str( trueVolume(i) ) '  A^3']);
    end
        if ORG_STRUC.varcomp == 0
                 disp(['     --   Total Volume of MOLs : ' num2str( ORG_STRUC.latVolume ) '  A^3']);
        end
%% Print IonDistances for atoms
    disp([' '])
    fprintf(1,'Default IonDistances for Atoms : ')
    for i=1:length(ORG_STRUC.atomType)
        fprintf(1,'%5s', megaDoof(ceil(ORG_STRUC.atomType(i))));
    end
    fprintf(1,'\n');

    for i=1:length(ORG_STRUC.atomType)
        fprintf(1,'    Minimum distances:    %5s: ', megaDoof(ceil(ORG_STRUC.atomType(i))));
        for j=1:length(ORG_STRUC.atomType)
            fprintf(1,'%4.2f  ', ORG_STRUC.minDistMatrice(i,j));
        end
        fprintf(1,'\n');
    end
    disp([' ']);
    disp(['===============================================================================================']);
    disp(['=   When you cannot generate structures for a long time, for the INPUT.txt                    =']);
    disp(['=    (1) check your "IonDistances" and "MolCenters" parameters, or reduce them;               =']);
    disp(['=    (2) increase the above default volume and set "Latticevalues" parameter.                 =']);
    disp(['===============================================================================================']);
    disp(['']);

    
elseif ORG_STRUC.varcomp == 1
    %001mode%
    if ORG_STRUC.dimension == 0
        for i = 1:length(ORG_STRUC.atomType)
	    numIons4volumecalc = diag(ones(1, length(ORG_STRUC.atomType))); %% we calc ORG_STRUC.latVolue from identity matrix instead ORG_STRUC.numIons
	    ORG_STRUC.latVolume(i)=calcDefaultVolume(numIons4volumecalc(i,:), ORG_STRUC.atomType, ORG_STRUC.ExternalPressure, 0, tempFolder, codeFolder);
	    clear numIons atomType
        end
    %end of 001mode%
    else
        for i = 1:size(ORG_STRUC.numIons,1)
            %ORG_STRUC.latVolume(i)=calcDefaultVolume(ORG_STRUC.numIons(i,:), ORG_STRUC.atomType, ORG_STRUC.ExternalPressure, 0);
	        ORG_STRUC.latVolume(i)=calcDefaultVolume(ORG_STRUC.numIons(i,:), ORG_STRUC.atomType, ORG_STRUC.ExternalPressure, 0, tempFolder, codeFolder);
            clear numIons atomType
        end
    end
else
    %ORG_STRUC.latVolume = calcDefaultVolume(ORG_STRUC.numIons, ORG_STRUC.atomType, ORG_STRUC.ExternalPressure, 0);
	ORG_STRUC.latVolume = calcDefaultVolume(ORG_STRUC.numIons, ORG_STRUC.atomType, ORG_STRUC.ExternalPressure, 0, tempFolder, codeFolder);
end

