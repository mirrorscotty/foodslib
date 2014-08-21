     function [Fts] =  tsolids(Fcomp,Fctype)

%  TSOLIDS  Total solids of stream based on composition  
%  
%           Synopsis: MATLAB function which sums all the non-water  
%           component mass fractions (w/w) to determine the total 
%           solids (w/w).    
%  
%           Call: tsolids(Fcomp,Fctype)
%   
%           Passed Variables: 
%            Fcomp - array of component mass fractions (w/w)  
%            Fctype - array of component type identifiers  
%    
%           Returned Variable:
%            Fts - total solids (w/w) 
 
%  *********************************************************************** 

%    Initialize food total solids (w/w)   
     Fts = 0;

%    For each component      
     for I = 1:length(Fctype)

%      If the component is not water  
       if (Fctype(I) ~= comptype('water')) & (Fctype(I) ~= comptype('air')) 

%        Sum total solids content of food (w/w)
         Fts = Fts + Fcomp(I); 

       end
     end
 


