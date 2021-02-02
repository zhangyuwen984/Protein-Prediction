function checkPOTCARS(atomType, abinitioCode, var_comp, homePath)

% checks whether potentials correspond to atom types in INPUT.txt
% quit and display a warning, if not

cd('Specific')
try
 for i = 1 : length(abinitioCode)
  if (abinitioCode(i) == 1) | (abinitioCode(i) == 12) % vasp and mol. vasp
    if var_comp == 1
      for j = 1 : length(atomType)
        if exist(['POTCAR_' megaDoof(atomType(j))]) == 0
          dispError(['ERROR: File POTCAR_' megaDoof(atomType(j)) ' is missing! The calculation has to stop!'], 'Please check manual how to set up the variable composition calculation with VASP.', homePath);
        end
      end
      break;
    else   % fixcomp

      if exist(['POTCAR_' num2str(i)]) == 0
          dispError(['ERROR: File POTCAR_' num2str(i) ' is missing! The calculation has to stop!'], ' ', homePath);
      end      

%      [tmp, Elements] = unix(['./getStuff POTCAR_' num2str(i) ' TITEL 5']);
      Elements=python_uspex([homePath, '/FunctionFolder/getInput.py'],'-f POTCAR -b TITEL -c 4');
      eos = findstr(Elements, sprintf('\n'));
      if isempty(eos)
          dispError(['ERROR: File POTCAR_' num2str(i) ' seems to be empty! The calculation has to stop!'], ' ', homePath);
      end

      if length(eos) ~= length(atomType)
          dispError(['ERROR: File POTCAR_' num2str(i) ' does not have the correct number of atomic potentials! The calculation has to stop!'], ' ', homePath);
      end

      for j = 1 : length(atomType)
        if j == 1
          Element = Elements(1:eos(1)-1);   % read the element name and remove '\n' at the end
        else
          Element = Elements(eos(j-1)+1:eos(j)-1);   % read the element name and remove '\n' at the end
        end
% Note, we need to deal with names like "Li_sv"
        if ~isempty(findstr(Element ,'_'))
          Element = Element(1:findstr(Element ,'_')-1);
        end
        if ~strcmp(Element, megaDoof(atomType(j)))
          dispError(['File POTCAR_' num2str(i) ' has potential for ' megaDoof(atomType(j)) ' missing! The calculation has to stop!'], 'Please check manual how to set up the calculation with VASP.', homePath);
        end
      end      
    end
  end
 end
catch
end

cd(homePath)

function dispError(str1, str2, homePath)
cd(homePath)
disp(str1);
disp(str2);
[nothing, nothing] = unix(['echo ' str1 ' >> Error.txt']);
[nothing, nothing] = unix(['echo ' str2 ' >> Error.txt']);
quit
