function createMake(make_template, angles)
% $Rev: 584 $
% $Author: maxim $
% $Date: 2014-08-30 10:27:47 +0400 (Sat, 30 Aug 2014) $

% Creates input.make file for Tinker execution.

    output_make = 'input_stable.make';

    fin  = fopen(make_template, 'rb');
    fout = fopen(output_make,'w' );

    tline = fgets(fin);
    counter = 0;
    format = '%*s %f %f';


    while ischar(tline)
        num_list = zeros();
        %disp(tline);
        outline = tline;
        if findstr(tline, '0.001') > 0 & findstr(tline, '-0.001') > 0
            counter = 1;
        end

        if counter >= 1
            num_list = sscanf(tline, format, [1, inf]);
            if ~isempty(num_list)
                phi_psi = angles(counter, :);
                phi = phi_psi(1);
                psi = phi_psi(2);
                % First replace psi since it's with '-' sign:
                outline = strrep(tline  , strcat('-', num2str(sprintf('%.3f', counter/1000)), ' '), num2str(sprintf('%8.3f', psi)));
                % Then replace phi:
                outline = strrep(outline, strcat(' ', num2str(sprintf('%.3f', counter/1000)), ' '), num2str(sprintf('%8.3f', phi)));
                counter = counter + 1;
            end
        end
        fprintf(fout,'%s', outline);
        tline = fgets(fin);     % Go to the next row
    end
    fclose(fin);    % End of file reading
    fclose(fout);

    return
end