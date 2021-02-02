function [] = safesave(filename0,varargin)

% FUNCTION:
%   [] = safesave(filename,varargin)
%
% DESCRIPTION:
%   Saves specified input arguments to a mat-file
%   and creates a backup of the original file if necessary
%
% INPUT ARGUMENTS:
%   filename = string containing desired filename (either with or without '.mat' extension)
%   varargin = comma-separated list of variables to save (NOT strings: the function uses the original variable names)
%

% Process input arguments
if nargin < 2
    error('safesave requires at least two input arguments: a file name and one or more variables to save to that file...')
end

% Make sure that extension is included in filename (*.mat)

if ~strcmp( filename0(end-3:end), '.mat' )
    fileFullName = [filename0 '.mat']; % If extension not present, add it.
else
    fileFullName= filename0;
end

% Assign caller's variable names
nonamecount = 0;
varstring = '';

    for i = 1:length(varargin)
        if ~isempty(inputname(1+i))
           varname = inputname(1+i);
        else
           nonamecount = nonamecount + 1;
           varname = ['noname' num2str(nonamecount,'%03.0f')]
        end
         varstring = [varstring ',''' varname ''''];
        eval([varname ' = varargin{i};']);
     end


% Assign caller's variable names

mark=strfind(fileFullName, '/');
if mark
   filename=fileFullName(mark(end)+1:end);
   pathname=fileFullName(1:mark(end)-1);
else
   filename=fileFullName;
   pathname=['./'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([pathname '/matfileLocker'])
	disp([' '])
	disp(['ERROR: '])
	disp(['    ', fileFullName, ' mat file is locked by USPEX! Please check your calculation !']);
	disp(['    For continuing the USPEX calculation, please delete the file :', pathname, 'matfileLocker'  ]);
	exit;
end


loadFailed    =0;
loadBackFailed=0;
if exist(fileFullName)
   try
      whoIsThere = whos('-file',fileFullName);
   catch
      loadFailed=1;
   end
   try
      whoIsThere = whos('-file',[fileFullName,'.backup']);
   catch
      loadBackFailed=1;
   end


   if (0==loadFailed) && (1==loadBackFailed)
      [nothing, nothing] = unix( ['cp ', fileFullName '        ' , fileFullName '.' '`date +"%b%d-%T"`' ]) ;
      [nothing, nothing] = unix( ['cp ', fileFullName '        ' , fileFullName '.backup']);
      USPEXmessage(1001, '', 0);      
   end
   if (1==loadFailed) && (0==loadBackFailed)
      [nothing, nothing] = unix( ['cp ', fileFullName '.backup ' , fileFullName '.' '`date +"%b%d-%T"`' ]) ;
      [nothing, nothing] = unix( ['cp ', fileFullName '.backup ' , fileFullName]);
      USPEXmessage(1002, '', 0);
   end
   if (loadFailed==1) && (loadBackFailed==1)
      [nothing, nothing] = unix(['echo -e "Interupt mat files " >> ERROR'])
      exit;
   end
end


tryTimes=10;

[nothing, nothing] = unix(['touch ', pathname '/matfileLocker']);
for iC = 1:tryTimes
   saveFailed=0;

   % Save data to new mat-file
   %filename0
   savestring = ['save(''' fileFullName '''' varstring ')']; % String with save() command
   %savestring
   eval(savestring);

   try
     whoIsThere = whos('-file',fileFullName);
   catch
     saveFailed=1;
   end

   if saveFailed==1
      if iC < 10
         USPEXmessage(1003, '', 0);
         pause(3);
      else
        [nothing, nothing] =  unix(['echo -e "Ooops, your file system met a serious problem! USPEX exsiting..." >> Warning.log']);
	% [nothing, nothing] = unix(['cp ', fileFullName, '.backup ' , fileFullName ]) ;
         exit;
      end
   else
      [nothing, nothing] = unix(['cp ', fileFullName, ' ' , fileFullName, '.backup' ]);
      [nothing, nothing] = unix(['rm ' pathname '/matfileLocker'] );
      return;
   end
end
