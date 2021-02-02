function [result status] = python_uspex(varargin)
%PYTHON Execute Perl command and return the result.
%   PYTHON(PYTHONFILE) calls python script specified by the file PYTHONFILE
%   using appropriate python executable.
%
%   PYTHON(PYTHONFILE,ARG1,ARG2,...) passes the arguments ARG1,ARG2,...
%   to the python script file PYTHONFILE, and calls it by using appropriate
%   python executable.
%
%   RESULT=PYTHON(...) outputs the result of attempted python call.  If the
%   exit status of python is not zero, an error will be returned.
%
%   [RESULT,STATUS] = PYTHON(...) outputs the result of the python call, and
%   also saves its exit status into variable STATUS.
%


  outputFormat=0;   %% string
% outputFomat=1;   %% numerical
    
cmdString = '';


% Add input to arguments to operating system command to be executed.
% (If an argument refers to a file on the MATLAB path, use full file path.)
for i = 1:1
    thisArg = varargin{i};
    if ~ischar(thisArg)
        error(message('MATLAB:python:InputsMustBeStrings'));
    end
    if i==1
        if exist(thisArg, 'file')==2
            % This is a valid file on the MATLAB path
            if isempty(dir(thisArg))
                % Not complete file specification
                % - file is not in current directory
                % - OR filename specified without extension
                % ==> get full file path
                thisArg = which(thisArg);
            end
        else
            % First input argument is PythonFile - it must be a valid file
            error(message('MATLAB:python:FileNotFound', thisArg));
        end
    end

    % Wrap thisArg in double quotes if it contains spaces
    if isempty(thisArg) || any(thisArg == ' ')
        thisArg = ['"', thisArg, '"']; %#ok<AGROW>
    end

    % Add argument to command string
    cmdString = [cmdString, ' ', thisArg]; %#ok<AGROW>
end

for i = 2:nargin
    if ischar(varargin{i})
       option = varargin{i};
       cmdString = [cmdString, ' ', option]; %#ok<AGROW>
	elseif isnumeric(varargin{i}) && (i==nargin) && (varargin{i}==1)
	   outputFormat=1;
	elseif isnumeric(varargin{i}) && (i==nargin) && (varargin{i}==0)
	   outputFormat=0;
    else
       error(message('MATLAB:python:WrongInputFormat'));
    end       
    % Add argument to command string
end
% Execute Perl script
if isempty(cmdString)
    error(message('MATLAB:python:NoPythonCommand'));
elseif ispc
    % PC
    pythonCmd = fullfile(matlabroot, 'sys\python\win32\bin\');
    cmdString = ['python -W ignore' cmdString];
    pythonCmd = ['set PATH=',pythonCmd, ';%PATH%&' cmdString];
    [status, result0] = dos(pythonCmd);
else
    % UNIX
    [status ignore] = unix('which python'); %#ok
    if (status == 0)
        cmdString = ['python -W ignore', cmdString];
        [status, result0] = unix(cmdString);
    else
        error(message('MATLAB:python:NoExecutable'));
    end
end

% Check for errors in shell command
if nargout < 2 && status~=0
   error('MATLAB:python:ExecutionError', ...
        'System error: %s \nCommand executed: %s', result0, cmdString);
end

% To remove some warnings
% Such as: "python: /home/qiang/source/matlab/bin/glnxa64/libz.so.1: no version information available (required by python)"
if ~isempty(findstr(result0,'python'))
%   disp([' ']);
%   disp(['Warning : some complains for the python.m call']);
%   disp(['  ', result0])
end

startPlace=findstr(result0,'<CALLRESULT>');
result0=result0(startPlace+13:end);


if outputFormat==0
   result=result0;
else
   result=str2num(result0);
   if ~isreal(result)
      error('MATLAB:python:WrongInputFormat', ...
                    '  Output Format error:  \n      %s  is not numeric. \n  Command executed: \n      %s', result0, cmdString);
   end
end
