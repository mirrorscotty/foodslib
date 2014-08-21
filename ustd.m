    function [Fuarray] = ustd(Fuarray) 

%  USTANDARD Standardize length of stream u-array 
%  
%            For use in unit operation S-functions 
%
%            Synopsis: MATLAB function which lengthens or shortens stream 
%             u-array to standard "usizemax" size. Converts stream 
%             temperature in K to C.  
%
%            Called from unit operation S-functions 
% 
%            Call: ustandard(Fuarray)  
%  
%            Passed Arguments: 
%             Fuarray - current stream array; may be in Fmass or Fme format
%             usizemax - Current maximum u-array size 
% 
%            Returned Argument:  
%             Fuarray - stream array updated to standard length and with 
%                       temperature variable in (C)  

% ***************************************************************************

%     Maximum u-array size
      usize = usizemax;   

%     Actual length of Fuarray
      Flength = length(Fuarray); 

%     Location in u-array of temperature 
      Ftemploc = ulocate('temp');  

%     If temperature is given in Fuarray 
      if (Flength >= Ftemploc)    
%       Convert temperature in K to C and reset u-array 
        Tref = 273.15;            % Reference temperature 
        FtempK = Fuarray(ulocate('temp'));  % Stream temperature (K)  

%       For positive temperature (K); convert K -> C 
        if (FtempK > 0) FtempC = FtempK - Tref; end   
        Fuarray(ulocate('temp')) = FtempC;
      end  

%     For stream u-array length less than usizemax 
      if (length(Fuarray) < usize)  

%       Extend array to usizemax length 
        Fuarray(usize,1) = 0;  

%     For stream u-array length greater than usizemax  
      elseif (length(Fuarray) > usize) 

%       Truncate u-array to usizemax  
        Fuarray = Fuarray(1:usize,1);   
      end   
  


 
