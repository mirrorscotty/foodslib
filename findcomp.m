
       function [Fcompi] = findcomp(Fcomp,Fctype,ctype)

%  FINDCOMP  Find mass fraction of specified component in component array 
%
%            Synopsis: MATLAB function that determines the mass fraction 
%            for the component specified. Returns a mass fraction of zero 
%            if the type is not found.
% 
%            Call: findcomp(Fcomp,Fctype,ctype) 
%   
%            Passed Argument:
%             Fcomp - component mass fractions (w/w) 
%             Fctype - component types
%             ctype - component for which a mass fraction is being sought  
%
%            Returned Argument:  
%             Fcompi - mass fraction (w/w) of component being sought  
%

%  ************************************************************************

%         Find index for component type  
          I = findval1(Fctype,ctype);  

          Fcompi = 0; 

%         If component is not found; return zero value 
          if isempty(I)
            Fcompi = 0;

%         If multiple incidences of component are found; sum mass fractions
          elseif (length(I) > 1) 
            for J = 1:length(I) 
              Fcompi = Fcompi + Fcomp(J);  
            end   

%         Else return mass fraction for component  
          else 
            Fcompi = Fcomp(I);   
          end  

   


