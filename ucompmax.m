     function[ucomp] = ucompmax()  

%  UCOMPMAX  Current number of components compositions in stream u-array  
%
%            Synopsis: MATLAB function which determines the maximum 
%            number of components compositions which can be defined 
%            for each stream u-array. Currently this value is 10.
% 
%            Call: ucompmax() 
%
%            Returned Argument:
%             ucomp - maximum number of components   

%  *************************************************************************  

%    First component composition index in stream u-array   
     ucomp1st = ulocate('comp1st');    

%    Last component composition index in stream u-array   
     ucomplast = ulocate('complast');  

%    Total number of components
     ucomp = ucomplast - ucomp1st + 1;   

