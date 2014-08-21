     function[status] = turn(descript)  

%  TURN      String interpreter for on/off status 
%
%            Synopsis: MATLAB function which turns a particular
%             sectional of FOODS-LIB on or off
% 
%            Call: turn(descript) 
%  
%            Passed Argument: 
%             descript - 'on' or 'off' (string)  
%                   
%            Returned Argument:
%             status - numerical representation of on/off (1,0)  

%  *************************************************************************  

%    Convert string to number

     if (strcmp(descript,'off'))
       status = 0;                 

     elseif (strcmp(descript,'on'))
       status = 1;          

     end  

