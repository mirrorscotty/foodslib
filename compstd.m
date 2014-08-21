    function [Ecompnew] = compstd(Ecompold) 

%  COMPSTD   Standardize length of component and component type arrays 
%  
%            For use in unit operation S-functions 
%
%            Synopsis: MATLAB function which lengthens or shortens stream 
%             stream component or component type arrays to standard "compmax" 
%             size. 
%
%            Called from unit operation S-functions 
% 
%            Call: compstd(Ecompold)  
%  
%            Passed Arguments: 
%             Ecompold - current stream component array; 
%                        may be in Fcomp or Fctype 
% 
%            Returned Argument:  
%             Ecompnew - stream component array updated to standard length 
%    

% ***************************************************************************

%     Maximum number of components   
      csize = compmax;   

%     Actual length of Ecompold  
      Elength = length(Ecompold); 

%     For stream component array length less than compmax 
      if (Elength < csize)  

%       Extend array to compmax length 
        Ecompnew = Ecompold;
        Ecompnew(csize) = 0;  

%     For stream component array length greater than compmax  
      elseif (Elength > csize) 

        if (Ecompold(1) < 1)      % Signal for component matrix
%       Sort array so that largest component mass fractions are retained  
%       Truncate u-array to usizemax  
        Fuarray = Fuarray(1:usize,1);   
        end 
      end   
  


 
