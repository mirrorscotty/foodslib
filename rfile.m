     function[filem] = rfile()  

%  RFILE     Creates user-specified results file for design simulation
%
%            Synopsis: MATLAB function which sets the filename for the 
%            FOODS-LIB model design results.
% 
%            Call: rfile() 
%
%            Returned Argument:
%             filem - user-specified file name 

%  *************************************************************************  

%   Print to results to specified file name
    global filem;
    if ~isstr(filem)

%   global filem;
%   filehere = isstr(filem);
%   if (filehere == 0)
      % CHANGE MADE 4/21/01: ALLOWS STRING INPUT
      filem = input('Enter filename for design results:   ','s');
      % END CHANGE
      fid = fopen(filem,'a+');
    else
      fid = fopen(filem,'a+');
    end
 


