function result = callAWK(scriptname, filename, varargin)


global ORG_STRUC

option = [];
if ~isempty(varargin)
   option=varargin{:};
end


USPEXPath = ORG_STRUC.homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
if findstr(uspexmode,'exe')
   [nothing,USPEXPath]= unix('echo -e $USPEXPATH');
   USPEXPath(end)=[];
end


awkScript=[USPEXPath, '/FunctionFolder/Tool/', scriptname];

%disp(['awk -f ' awkScript '  ' option '  '  filename ]);
[nothing, resultStr] = unix(['awk -f ' awkScript '  ' option '  '  filename ]);
try
   result=str2num(resultStr);
catch
   result=[];
   USPEXmesage(1010,resultStr,0);
   disp(resultStr);
end
