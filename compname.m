     function[cname] = compname(cindex)  

%  COMPNAME  Stream component name    
%
%            Synopsis: MATLAB function which returns index of stream
%            component type  
% 
%            Call: compname(cindex) 
%  
%            Passed Argument: 
%             cindex - component type index  
%                   
%            Returned Argument:
%             cname - component type name (string)  

%  *************************************************************************  

%    Assign component type index for specified information    

     if (cindex == 1)
       cname = 'Water        ';                  

     elseif (cindex == 2)
       cname = 'Protein      ';           

     elseif (cindex == 3)
       cname = 'Carbohydrate ';  

     elseif (cindex == 4)
       cname = 'Fiber        ';      

     elseif (cindex == 5)  
       cname = 'Fat          ';     

     elseif (cindex == 6) 
       cname = 'Ash          ';  

     elseif (cindex == 7)
       cname = 'Air          ';  

     end  
