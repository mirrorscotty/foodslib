     function [sys,x0] = feeds(t,x,u,flag,Fname)  

%  FEEDS     System input feed variables u-array re-assignment
%
%            For use with FOODS-LIB models   
%  
%            Synopsis: S-Function builds standard feed stream array. 
%            Extends u_array from 36 inputs to "usizemax" variables.  
%            U_array elements 37-"usizemax" are available for passing 
%            new information. Checks for pure stream (ie. water) and
%            verifies constant values.
%  
%            Call: feeds(t,x,u,flag,Fname)
%  
%            Passed Arguments:
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 information
%             Fname - stream name 
% 
%            Returned Arguments:
%             sys - SIMULINK S-function variable: for flag = 0,
%                   sys = sizes of parameters and initial conditions.
%                   For flag = 3, sys = system outputs which are the
%                   u-arrays for stream entering the model.
  
%  ********************* SIMULINK INTERFACE INFORMATION *********************
 
usize = usizemax;                   % Current u_array size for streams
noinputs = 36;                      % Current number of constant inputs
 
if flag == 0 
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, usize, noinputs, 0, 1];   
   x0 = ones(size(u));  

elseif abs(flag) == 1
   % if flag=1, return state derivatives, xDot
     disp('Set Simulation-Parameters Stop Time to Zero.')
     break;

elseif flag == 3 | flag == 9  
   % if flag=3, return system outputs, y

%    Feed Stream Variables 

%  ************************** MASS BALANCE LEVEL ****************************
 
     F = u(1);               % Flow rate (kg/h) 
     for I = 1:10
      Fctype(I,1) = u(2*I);    % Component fraction types   
      Fcomp(I,1) = u(2*I+1);   % Component fractions (w/w)  
     end  

     if (any(Fcomp < 0)) 
       errortxt = [10 'FEED COMPOSITION ERROR: At least one feed' 10 ...
       'component has a negative mass fraction. Check values entered ' 10 ...
       'in FOOD COMPOSITION workspaces.' 10];
       disp(errortxt); 
       return;
     end 

     Fh2o = 1 - tsolids(Fcomp,Fctype);      % Calculated water fraction (w/w)  

%    For case where the entered water mass fraction and calculated water
%    fraction do not match.  
     if ((Fh2o ~= Fcomp(1)) & (Fcomp(1) > 0)) 
       warntxt = [10 'FOOD COMPOSITION WARNING: The water mass fraction entered' 10 ...
       'for a feed stream did not bring the sum of the component' 10 ...
       'mass fractions to 1. A new value has been calculated.' 10]; 
       disp(warntxt);
     end 

     Fcomp(1) = Fh2o;

%    For case where sum of mass fractions (not water) exceeds 1 
     if (Fh2o < 0.0)
       errortxt = [10 'FOOD COMPOSITION WARNING: The total solids of a feed stream' 10 ...
       'is greater than 1. Check values entered in FOOD COMPOSITION' 10 ...
       'workspaces.' 10];  
       disp(errortxt); 
       return;
     end 

%    Check for multiple-component stream (pure water = 1) 
     flagmc = purechk(Fcomp,Fctype); 

%  ************************* ENERGY BALANCE LEVEL ***************************

     FtempC = u(22);           % Temperature (C)  
     Fpres = u(23);            % Pressure (MPa)  
     Fx    = u(24);            % Quality/State   

%    For pure water stream; verify energy balance level inputs    
     if (flagmc == 1)   
%      Energy and cost information u-array required for steams.m call  
       Fec = [FtempC;Fpres;Fx;0;0]; 
       Fuarray = steams(0,0,Fec,3); % call format for S-function   
       Fctype = Fuarray(ulocate('ctype1st'):ulocate('ctypelast'));
                                             % Component types  
       FtempC = Fuarray(ulocate('temp'));    % Temperature (C)
       Fpres = Fuarray(ulocate('pressure')); % Pressure (MPa)
       Fx    = Fuarray(ulocate('quality'));  % Quality/State
     end   

%  ********************  MICROBIAL DEATH KINETICS ********************* 

     FmDval = u(25);           % D value (min) 
     FmDtempC = u(26);         % D value temperature basis (C) 
     FmZval = u(27);           % Z value (C) 
     FmEa = u(28);             % activation energy (kcal/mol) 
     Fmload = u(29);           % microbial load (user defined units)   

%  ************************  QUALITY KINETICS ************************* 

     FqDval = u(30);           % D value (min)
     FqDtempC = u(31);         % D value temperature basis (C)
     FqZval = u(32);           % Z value (C)
     FqEa = u(33);             % activation energy (kcal/mol)  
     Fqload = u(34);           % quality load (user defined units)   

%  ************************* ECONOMICS LEVEL **************************

     Fcost = u(35);            % Cost ($/kg)  
     Finfl = u(36);            % Annual cost inflation rate 

%  ************** Re-ORGANIZATION & EXTENSION OF U-ARRAY ************** 

%    Re-organization of u-array  
     u(ulocate('mass')) = F;   
     Fc1st = ulocate('comp1st');  
     Fclast = ulocate('complast');  
     u(Fc1st:Fclast) = Fcomp;  
     Fct1st = ulocate('ctype1st');
     Fctlast = ulocate('ctypelast');
     u(Fct1st:Fctlast) = Fctype;  
     u(ulocate('temp')) = FtempC;
     u(ulocate('pressure')) = Fpres;
     u(ulocate('quality')) = Fx; 
     u(ulocate('microDval')) = FmDval; 
     u(ulocate('microDtemp')) = FmDtempC;
     u(ulocate('microZval')) = FmZval;
     u(ulocate('microEa')) = FmEa;
     u(ulocate('microload')) = Fmload; 
     u(ulocate('qualDval')) = FqDval;
     u(ulocate('qualDtemp')) = FqDtempC; 
     u(ulocate('qualZval')) = FqZval;
     u(ulocate('qualEa')) = FqEa;
     u(ulocate('qualload')) = Fqload;
     u(ulocate('cost')) = Fcost;
     u(ulocate('inflation')) = Finfl;

%    Extension of u-array to specified "usizemax" elements    
     for I = (length(u)+1):usize  
       u(I) = 0;       
     end  

     sys = [u;];   

%  Send feed cost information to stream cost/value data file 
   if flag == 9
     if (Fcost > 0)
       strmtype = 1;
       scfile(strmtype,F,Fcost,Finfl,Fname);
     end 
     sys = [];   
   end

else
   % Otherwise, no need to return anything
   % since this is continuous system
   sys = [];   

end
 
