     function[fileuo] = uofile(UOutil,UOsize,UOnum,UOecon,UOname)  

%  UOFILE    Creates user-specified unit operation cost data file 
%
%            Synopsis: MATLAB function which write the unit operation
%            cost data to user-specified file. Data is used in 
%            economic analysis.   
% 
%            Call: uofile(UOutil,UOsize,UOnum,UOecon,UOname) 
%
%            Passed Arguments: 
%             UOutil - unit operation utilities demand array
%             UOsize - unit operation size  
%             UOnum - number of unit operations of UOsize 
%             UOecon - unit operation economic array 
%             UOname - unit operation name  
%
%            Returned Argument:
%             fileuo - user-specified unit operation cost data file name 

%  *************************************************************************  

%   Print to results to specified file name
    global fileuo;
    if ~isstr(fileuo)
      % CHANGE MADE 4/21/01: ALLOWS STRING INPUT  
      fileuo = input('Enter filename for unit operation costs data:   ','s');
      % END CHANGE 4/21/01
      
      fid = fopen(fileuo,'a+');
    else
      fid = fopen(fileuo,'a+');
    end

    fprintf(fid,'%12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %12.11e %s %2.0f\n', ...
     UOutil(1),UOutil(2),UOutil(3),UOutil(4),UOutil(5),UOutil(6),UOsize,UOnum, ...
     UOecon(1),UOecon(2),UOecon(3),UOecon(4),UOecon(5),UOecon(6), ...
     UOecon(7),UOecon(8),UOecon(9),UOecon(10),UOname,-1);

 


