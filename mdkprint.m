 
   function[] = mdkprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY)
 
%  MDKPRINT  Thermal death-rate kinetics of microbes results printout   
%
%            Synopsis: MATLAB function which prints the microbial thermal
%            death-rate kinetics results table for a unit operation block; 
%            called only by unit operation S-functions.  
% 
%            Call: mdkprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY)
%
%            Passed Arguments: 
%             UNAME -  unit operation name 
%             SINAMES - matrix of names of all system inputs being assessed 
%             SONAMES - matrix of names of all system outputs being assessed 
%             SNNAMES - matrix of names of all system intermediates being assessed 
%             SIARRAY - system input array; each column of the array represents 
%              a stream specification which includes: D value, D temperature basis, 
%              Z value, Ea value, and microbial load (Fmdk array format) 
%             SOARRAY - system output array; each column of the array represents 
%              a stream specification which includes: D value, D temperature basis, 
%              Z value, Ea value, and microbial load (Emdk array format), plus 
%              processing time and bulk temperature  
%             SNARRAY - system inermediate array; each column of the array represents
%              a stream specification which includes: D value, D temperature basis,
%              Z value, Ea value, and microbial load (Nmdk array format), plus 
%              estimated minimum processing time and bulk temperature  

% ******************************************************************************

%   Print results to user specified filename which was determined in rfile.m 
    global filem;
    % THIS IS WHERE THE CHANGE WAS MADE 4/21/01
    fid = fopen(filem,'a+');
    %fid = filem;
    % END CHANGE 4/21/01

%   Limiting or sizing all stream names to 15 characters 
    SINAMES = namesize(SINAMES,15); 
    SONAMES = namesize(SONAMES,15);  
    if (size(SNNAMES,1) >= 1)
      INT = 1;
      SNNAMES = namesize(SNNAMES,15);
    else INT = 0;
    end
 
    fprintf(fid,'     \n');
    fprintf(fid,'     \n');
    fprintf(fid,'%s MICROBIAL KINETICS RESULTS\n',UNAME);
    fprintf(fid,'     \n');
    fprintf(fid,'Stream               D        D Temp       Z          Ea         Time      Temp       Load       \n');
    fprintf(fid,'Name                           Base \n');
    fprintf(fid,'                   (min)       (C)        (C)     (kJ/g-mole)    (min)      (C)   \n');
    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
    fprintf(fid,'SYSTEM INPUTS:\n');
    ninputs = size(SIARRAY,2); 
    for I = 1:ninputs   


      Dvalm = SIARRAY(1,I);             % D value (min)
      DtempCm =  SIARRAY(2,I);          % D value temperature basis (C)      
      Zvalm =  SIARRAY(3,I);            % Z value (C)
      Eam =  SIARRAY(4,I);              % activation energy (kcal/mol)
      Noin = SIARRAY(5,I);              % microbial load
  
      fprintf(fid,'%s  %8.4f     %6.2f   %6.2f     %8.2f       ---       ---      %8.3e \n', ...
       SINAMES(I,:),Dvalm,DtempCm,Zvalm,Eam,Noin);
    end   

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM OUTPUTS:\n');
    noutputs = size(SOARRAY,2);     
    for I = 1:noutputs      

      Dvalm = SOARRAY(1,I);             % D value (min)
      DtempCm = SOARRAY(2,I);          % D value temperature basis (C)
      Zvalm =  SOARRAY(3,I);            % Z value (C)
      Eam =  SOARRAY(4,I);              % activation energy (kcal/mol)
      Nout = SOARRAY(5,I);              % microbial load
      tmin = SOARRAY(6,I);              % time at bulk temperature (min)
      Tbulk = SOARRAY(7,I);             % bulk temperature (C)  

      fprintf(fid,'%s  %8.4f     %6.2f   %6.2f     %8.2f   %8.4f   %8.2f    %8.3e \n', ...
       SONAMES(I,:),Dvalm,DtempCm,Zvalm,Eam,tmin,Tbulk,Nout);
    end  
    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');

    if (INT == 1) 

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM INTERMEDIATES:\n');
    noutputs = size(SNARRAY,2);    
    for I = 1:noutputs     
 
      Dvalm = SNARRAY(1,I);             % D value (min)
      DtempCm = SNARRAY(2,I);          % D value temperature basis (C)
      Zvalm =  SNARRAY(3,I);            % Z value (C)
      Eam =  SNARRAY(4,I);              % activation energy (kcal/mol)
      Nout = SNARRAY(5,I);              % microbial load
      tmin = SNARRAY(6,I);              % time at bulk temperature (min)
      Tbulk = SNARRAY(7,I);             % bulk temperature (C) 
 
      fprintf(fid,'%s  %8.4f     %6.2f   %6.2f     %8.2f   %8.4f   %8.2f    %8.3e \n', ...
       SNNAMES(I,:),Dvalm,DtempCm,Zvalm,Eam,tmin,Tbulk,Nout);
    end % for I
    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');

    end % intermediates print


    fprintf(fid,'     \n');
    fprintf(fid,'     \n');

    endtxt = [10 'END OF (',UNAME,') MICROBIAL KINETICS PRINT TO FILE' 10];
    disp(endtxt);


