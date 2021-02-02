function name = makeLabel(atomType, numIons, position)
name = '';
for loop1 = 1:size(numIons, 1)
    tmp = '';
    for loop2 = 1:size(numIons, 2)
        N_atom = numIons(loop1,loop2);
        if N_atom > 0
            tmp= [tmp megaDoof(atomType(loop2))];
            if N_atom > 1
                label=['{' num2str(N_atom) '}'];
                tmp= [tmp '_' label];
            end
        end
    end
    name{loop1} = tmp;
    text(position(loop1,1), position(loop1,2), tmp, ...
	     'HorizontalAlignment', 'Center','FontSize',15);
end
