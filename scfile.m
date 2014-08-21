     function[filesc] = scfile(Stype,S,Scost,Sinfl,Sname)  

%  SCFILE    Creates user-specified stream cost/value data file 
%
%            Synopsis: MATLAB function which writes stream cost/value 
%            data to user-specified filename. Data is for economic 
%            analysis.
% 
%            Call: scfile(Stype,S,Scost,Sinfl,Sname)  
%
%            Passed Arguments:
%             Stype - stream type: raw material (1); product (2);
%                       waste (3) 
%             S - stream flow rate 
%             Scost - stream cost ($/kg) 
%             Sinfl - stream  cost inflation rate (d percent)
%             Sname - stream name  
%
%            Returned Argument:
%             filesc - user-specified stream cost/value data file name 

%  *************************************************************************  

%   Print cost/value data to specified file 
    global filesc;
    if ~isstr(filesc)
      % CHANGE MADE 4/21/01: ALLOWS STRING INPUT
      filesc = input('Enter filename for stream costs/value data:   ','s');
      % END CHANGE 4/21/01
      
      fid = fopen(filesc,'a+');
    else
      fid = fopen(filesc,'a+');
    end

    fprintf(fid,'%2.0f %12.11f %12.11f %12.11f %s %2.0f\n',Stype,S,Scost, ...
      Sinfl,Sname,-1);

 


