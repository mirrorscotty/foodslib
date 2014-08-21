     function[valuey] = inflate(value1,rate,year,dorc)

%  INFLATE  Inflates value to year N value  
%
%           Synopsis: MATLAB function which determines the value
%           in year Y based on the year 1 value and the rate of
%           inflation  
%
%           Call: inflate(value1,rate,year,dorc)
% 
%           Passed Arguments: 
%            value1 - year 1 value 
%            rate - rate of inflation 
%            year - year to which year 1 value is inflated 
%            dorc - discrete (0) or compounded (1) inflation 
%

%  **********************************************************************

%    Determines the value in year Y based on the value in year 1
%    and the rate.  
 
     if (dorc == 0) 
%      Discrete inflation of stream cost
       valuey = value1*(1+rate)^(year-1);  
     else  
%      Compounded inflation of stream cost
       valuey = value1*exp(rate*(year-1));   
     end   

     
