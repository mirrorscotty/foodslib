
    function[flagmc] = purechk(Fcomp,Fctype) 

%  PURECHK    Checks stream composition for multiple components 
%
%             Synopsis: MATLAB function to determine if a stream contains
%               only water or only air-and-water 
%
%             Call:  purechk(Fcomp,Fctype)  
% 
%             Passed Arguments: 
%              Fcomp - array of component mass fractions (w/w)
%              Fctype - array of component types 
%
%             Returned Argument:  
%              flagmc - flag for pure water stream (= 1)
%                       flag for air-water stream (= 2)

% ************************************************************************ 

%    Flag for multi-component food stream
     flagmc = 0;

%    Index for water component  
     xh2o = [];  % initialize
     xh2o = findval1(Fctype,comptype('water'));

%    Check for pure water stream  
     if ~isempty(xh2o)   
%      Water fraction (w/w) equals 1
       if (abs(Fcomp(xh2o) - 1) < 1e-10)
%        Flag for pure water stream
         flagmc = 1;  
       end  
     end  

%    Index for air component 
     xair = [];  % initialize
     xair = findval1(Fctype,comptype('air'));  

%    Check for pure air-water stream
     if ~isempty(xair) & ~isempty(xh2o)  
%      Water and air fractions (w/w) sum to 1
       if (abs(Fcomp(xh2o)+Fcomp(xair)-1) < 1e-10)
%        Flag for air-water stream
         flagmc = 2;  
       end  
     end  
 




