
       function [index] = findval1(array,value)

%  FINDVAL1  Find index of all elements in 1D array with value specified
%
%            Synopsis: MATLAB function that determines the index locations 
%            that contains that value specified. Returns a vector;
%            the vector is empty if value is not found.
% 
%            Call: findval1(array,value) 
%   
%            Passed Argument:
%             array - 1-Dimensional array to be searched
%             value - value to be found in array 
%
%            Returned Argument:  
%             index - index of value
%

%  ************************************************************************

       J = 1;  

%      Search length of array for value 
       for I = 1:length(array)
         if (array(I) == value)
%          Set index to location of value 
           index(J) = I; 
           J = J + 1;  
         end
       end 

       if (J == 1) index = []; end    


