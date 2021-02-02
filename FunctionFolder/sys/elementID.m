function ID = elementID(name)

name = strtrim(name);  %remove blank ending if necessary
ID = [];
if strcmp(name, 'H.5') | strcmp(name, 'H.75') | strcmp(name, 'H1.25') | strcmp(name, 'H1.5')
   ID = 1;
else
   for i=1:105  %to find the index for the element
      if strcmp(lower(name), lower(megaDoof(i))) | strcmp(lower(name), lower(elementFullName(i)))
         ID = i;
         break;
      end
   end
end

if isempty(ID)
   warningStr = ['Element ' name ' does not exist. Please go back to check it again.\n'];
   USPEXmessage(0, warningStr, 0);
   quit
end
