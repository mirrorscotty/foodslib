     function [sys,x0] = products(t,x,u,flag,Pname)

%  PRODUCTS  System output product variables u-array re-assignment 
% 
%            For use with FOODS-LIB models
%    
%            Synopsis: S-Function adds economic information to standard 
%            product stream array.
%    
%            Call: products(t,x,u,flag,Pname)
%    
%            Passed Arguments:
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                  and product economic information
%             Pname - Product stream name; entered as Parameter 
%    
%            Returned Arguments:
%             sys - SIMULINK S-function variable: for flag = 0,
%                   sys = sizes of parameters and initial conditions.
%                   For flag = 3, sys = system outputs which are the
%                   u-arrays for stream entering the model.
 
%  ********************* SIMULINK INTERFACE INFORMATION *********************

usize = usizemax;                   % Maximum stream u-array size  
noinputs = 2;                       % Number of product income inputs  
 
if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, usize, usize+noinputs, 0, 1];
   x0 = ones(size(u));

elseif abs(flag) == 1
   % if flag=1, return state derivatives, xDot
     disp('Set Simulation-Parameters Stop Time to Zero.')
     break;

elseif flag == 3 | flag == 9 
   % if flag=3, return system outputs, y

%  **************************** ECONOMICS LEVEL ****************************

     P = u(ulocate('mass'));       % Product mass flow rate (kg/h)  
     Pincome = u(51);              % Income ($/kg)
     Pinfl = u(52);                % Annual income inflation rate (d %)

%  ***************** Re-ORGANIZATION & EXTENSION OF U-ARRAY ****************
 
%    Re-organization of u-array
     u(35) = Pincome; 
     u(36) = Pinfl;   

     sys = [u(1:usize)];

%  Send product value information to stream cost/value data file 
   if flag == 9
     if (Pincome > 0)
       strmtype = 2;
       scfile(strmtype,P,Pincome,Pinfl,Pname);
     end 
     sys = []; 
   end

 
else
   % Otherwise, no need to return anything
   % since this is continuous system
   sys = []; 
 
end

 
