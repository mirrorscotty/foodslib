   function[NAMESNEW] = namesize(NAMES,NSIZE)

%  NAMES     Standard name length for printing tables  
%
%            Synopsis: MATLAB function which standardizes the size of
%            all name strings in a matrix.  
%
%            Call: namesize(NAMES,NSIZE)  
%
%            Passed Arguments:
%             NAMES - matrix of names 
%             NSIZE - standard size of all names  
%
%           Returned Arguement: 
%             NAMENEW - matrix of names all in standardized size             

% ******************************************************************************

%   Limiting or sizing all stream names to standard number of characters
    NAMES(1,(NSIZE+1)) = ' ';  % Extend the names array to 1+ standard 
                               % characters by adding a space

    for I = 1:size(NAMES,1)    % For all names less than 15 characters, extend
      for J = 1:NSIZE          %   name length to 15 char by adding spaces
        if (abs(NAMES(I,J)) == 0) NAMES(I,J) = ' '; end
      end
    end 
 
    NAMESNEW = NAMES(:,1:NSIZE); % Truncate the names array to standard
                                 % number of characters




