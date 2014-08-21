     function[uindex] = ulocate(info)  

%  ULOCATE   Locate index in u-array for specified information  
%
%            Synopsis: MATLAB function which locates the u-array index  
%            which holds the information type specified.
% 
%            Call: ulocate(info) 
%  
%            Passed Argument: 
%             info = information type for which index is being saught
%                    (string)  
%                   
%            Returned Argument:
%             uindex - u-array index  

%  *************************************************************************  

%    Locate uarray index for specified information    

%  ************************** MASS BALANCE LEVEL ****************************

     if (strcmp(info,'mass'))
       uindex = 1;          % Mass flow rate (kg/h)        

     elseif (strcmp(info,'comp1st'))
       uindex = 2;          % First stream component (w/w) 

     elseif (strcmp(info,'complast'))
       uindex = 11;         % Last stream component (w/w) 

     elseif (strcmp(info,'ctype1st'))
       uindex = 12;         % First stream component type 

     elseif (strcmp(info,'ctypelast'))  
       uindex = 21;         % Last stream component type  

%  *********************** ENERGY BALANCE LEVEL ***************************

     elseif (strcmp(info,'temp')) 
       uindex = 22;         % Temperature (C)

     elseif (strcmp(info,'pressure'))
       uindex = 23;         % Pressure (MPa) 

     elseif (strcmp(info,'quality'))  
       uindex = 24;         % Quality/State 

%  ***********************  MICROBIAL DEATH KINETICS *********************

     elseif (strcmp(info,'microDval'))   
       uindex = 25;          % D value (min)

     elseif (strcmp(info,'microDtemp'))   
       uindex = 26;           % D value temperature basis (C) 

     elseif (strcmp(info,'microZval'))   
       uindex = 27;           % Z value (C)

     elseif (strcmp(info,'microEa'))   
       uindex = 28;            % activation energy (kcal/mol) 

     elseif (strcmp(info,'microload'))   
       uindex = 29;           % microbial load (user defined units)

%  ***************************  QUALITY KINETICS *************************

     elseif (strcmp(info,'qualDval'))    
       uindex = 30;           % D value (min) 
 
     elseif (strcmp(info,'qualDtemp'))   
       uindex = 31;           % D value temperature basis (C) 
 
     elseif (strcmp(info,'qualZval'))    
       uindex = 32;           % Z value (C) 
 
     elseif (strcmp(info,'qualEa'))    
       uindex = 33;           % activation energy (kcal/mol) 

     elseif (strcmp(info,'qualload'))    
       uindex = 34;            % quality load (user defined units) 

%  ************************* ECONOMICS LEVEL  ****************************

     elseif (strcmp(info,'cost'))     
       uindex = 35;            % Cost ($/kg)
 
     elseif (strcmp(info,'inflation'))     
       uindex = 36;            % Annual cost inflation rate
 
     end  

