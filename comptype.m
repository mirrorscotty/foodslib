     function[cindex] = comptype(info)  

%  COMPTYPE  Index for component type   
%
%            Synopsis: MATLAB function which returns index of stream
%            component type  
% 
%            Call: comptype(info) 
%  
%            Passed Argument: 
%             info = information type for which index is being saught
%                    (string)  
%                   
%            Returned Argument:
%             cindex - component type index  

%  *************************************************************************  

%    Assign component type index for specified information    

     if (strcmp(info,'water'))
       cindex = 1;                  

     elseif (strcmp(info,'protein'))
       cindex = 2;           

     elseif (strcmp(info,'carbohydrate'))
       cindex = 3;  

     elseif (strcmp(info,'fiber'))
       cindex = 4;      

     elseif (strcmp(info,'fat'))  
       cindex = 5;     

     elseif (strcmp(info,'ash')) 
       cindex = 6;  

     elseif (strcmp(info,'air'))
       cindex = 7;  

     end  
