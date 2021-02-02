function constraintsOK = distanceCheck(Coordinates, Lattice, composition, minDistMatrice)

% USPEX Version 4.4.3
% Change: only check if actual distance required

global ORG_STRUC

if isempty(find(minDistMatrice>0.0001))
    constraintsOK = 1
elseif size(Coordinates,1) ~= sum(composition)
    constraintsOK = 0
    status = 'Size of coordinates <> sum(numIons)!'
    debug1 = Coordinates
    debug2 = composition
else
    dist = [];
    constraintsOK = 1;
    breakCounter = 0;

    ionChange = zeros(1,length(composition));
    for ind = 1:length(composition)
      ionChange(ind) = sum(composition(1:ind));
    end    

    if length(Lattice)==6
        Lattice = latConverter(Lattice);
    end

    % it is (a lot) easier to calculate the distances between ions, if we first transform the Coordinates 
    % into a matrice, where each colomn represents the coordinates on one lattice vector (basis)
    if ~isempty(find(size(Coordinates)==1))
        Coordinates_temp = Coordinates;
        Coordinates = zeros(3,sum(composition));
        Coordinates(:) = Coordinates_temp(:);
        Coordinates = Coordinates';
    end

    minimalDistance = 10;

    for loop = 1:(size(Coordinates,1)-1)

        if constraintsOK    %%%

            for who = (loop+1):size(Coordinates,1)
               row = find(ionChange>(loop-0.5));
               colomn = find(ionChange>(who-0.5));

               % since the unit cell gets repeated, the maximal distance
               % between two ions (on one lattice vector) is 0.5 (normalized)
               check = (Coordinates(loop,:)-Coordinates(who,:));
               if ORG_STRUC.dimension == 0
                  dist = sqrt(sum((check*Lattice).^2));
                  if dist < minDistMatrice(row(1),colomn(1))
                    constraintsOK = 0;
                    break   %%%
                  end
               else
                for x = 1 : 3
                 for y = 1 : 3
                  for z = 1 : 3
                    dist = sqrt(sum(((check+[x-2,y-2,z-2])*Lattice).^2));
                    if dist<minimalDistance
                       minimalDistance = dist;
                    end
                    if dist < minDistMatrice(row(1),colomn(1))

              % set constraintsOK to zero, so that we know the constraints where violated at least once.
                        constraintsOK = 0;
                        %  % with breakCounter we could count the number of times the
                        %  % constraints have been violated. This information may be
                        %  % of use, depending on the algorithm. In this version it is
                        %  % not used, but if you are changing the code and want
                        %  % this information, you need to delete all lines indicated
                        %  % with  %%% after the code and uncomment the following
                        %  % line
                        %                 breakCounter = breakCounter + 1;

                        break   %%%
                      end
                     end
                    end
                end
               end

               row =[];
               colomn = [];
               dist = [];
               check = [];
               where = [];
            end
        else                %%%
            break           %%%
        end                 %%%

    end
end

if ORG_STRUC.dimension == 0 | ORG_STRUC.dimension == -2
 
 if checkConnectivity(Coordinates, Lattice, composition) == 0
   constraintsOK = 0;
 end
end
