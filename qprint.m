 
   function[] = qprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY)
 
%  QPRINT    Quality kinetics results printout   
%
%            Synopsis: MATLAB function which prints the quality kinetics 
%            results table for a unit operation block; called only by unit 
%            operation S-functions.  
% 
%            Call: qprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY)
%
%            Passed Arguments: 
%             UNAME -  unit operation name 
%             SINAMES - matrix of names of all system inputs being assessed 
%             SONAMES - matrix of names of all system outputs being assessed 
%             SONAMES - matrix of names of all system intermediates being assessed 
%             SIARRAY - system input array; each column of the array represents 
%              a stream specification which includes: D value, D temperature basis, 
%              Z value, Ea value, and quality load (Fqual array format) 
%             SOARRAY - system output array; each column of the array represents 
%              a stream specification which includes: D value, D temperature basis, 
%              Z value, Ea value, and quality load (Equal array format), plus
%              estimated maximum processing time and temperature 
%             SNARRAY - system intermediate array; each column of the array represents
%              a stream specification which includes: D value, D temperature basis,
%              Z value, Ea value, and quality load (Nqual array format), plus 
%              estimated maximum processing time and temperature  

% ******************************************************************************

%   Print results to user specified filename which was determined in rfile.m 
    global filem;
    % THIS IS WHERE THE CHANGE WAS MADE 4/21/01
    fid = fopen(filem,'a+');
    %fid = filem;
    % END CHANGE 4/32/01

%   Limiting or sizing all stream names to 15 characters 
    SINAMES = namesize(SINAMES,15); 
    SONAMES = namesize(SONAMES,15);   
    SONAMES = namesize(SONAMES,15); 
    if (size(SNNAMES,1) >= 1)
      INT = 1;
      SNNAMES = namesize(SNNAMES,15);
    else INT = 0;
    end

    fprintf(fid,'     \n');
    fprintf(fid,'     \n');
    fprintf(fid,'%s QUALITY KINETICS RESULTS\n',UNAME);
    fprintf(fid,'     \n');
    fprintf(fid,'Stream               D        D Temp       Z          Ea         Time      Temp       Load       \n');
    fprintf(fid,'Name                           Base \n');
    fprintf(fid,'                   (min)       (C)        (C)     (kJ/g-mole)    (min)      (C)   \n');
    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
    fprintf(fid,'SYSTEM INPUTS:\n');
    ninputs = size(SIARRAY,2); 
    for I = 1:ninputs   


      Dvalq = SIARRAY(1,I);             % D value (min)
      DtempCq =  SIARRAY(2,I);          % D value temperature basis (C)      
      Zvalq =  SIARRAY(3,I);            % Z value (C)
      Eaq =  SIARRAY(4,I);              % activation energy (kcal/mol)
      Qoin = SIARRAY(5,I);              % quality load
 
      fprintf(fid,'%s  %8.4f     %6.2f   %6.2f     %8.2f       ---       ---      %8.3e \n', ...
       SINAMES(I,:),Dvalq,DtempCq,Zvalq,Eaq,Qoin);
    end   

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM OUTPUTS:\n');
    noutputs = size(SOARRAY,2);     
    for I = 1:noutputs      

      Dvalq = SOARRAY(1,I);             % D value (min)
      DtempCq = SOARRAY(2,I);          % D value temperature basis (C)
      Zvalq =  SOARRAY(3,I);            % Z value (C)
      Eaq =  SOARRAY(4,I);              % activation energy (kcal/mol)
      Qout = SOARRAY(5,I);              % quality load
      tmin = SOARRAY(6,I);              % time at bulk temperature (min)
      Tbulk = SOARRAY(7,I);             % bulk temperature (C)
 
      fprintf(fid,'%s  %8.4f     %6.2f   %6.2f     %8.2f   %8.4f   %8.2f    %8.3e \n', ...
       SONAMES(I,:),Dvalq,DtempCq,Zvalq,Eaq,tmin,Tbulk,Qout);
    end  
    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');

    if (INT == 1) 

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM INTERMEDIATES:\n');
    noutputs = size(SNARRAY,2);
    for I = 1:noutputs

      Dvalq = SNARRAY(1,I);             % D value (min)
      DtempCq = SNARRAY(2,I);          % D value temperature basis (C)
      Zvalq =  SNARRAY(3,I);            % Z value (C)
      Eaq =  SNARRAY(4,I);              % activation energy (kcal/mol)
      Qout = SNARRAY(5,I);              % quality load
      tmin = SNARRAY(6,I);              % time at bulk temperature (min)
      Tbulk = SNARRAY(7,I);             % bulk temperature (C)

      fprintf(fid,'%s  %8.4f     %6.2f   %6.2f     %8.2f   %8.4f   %8.2f    %8.3e \n', ...
       SNNAMES(I,:),Dvalq,DtempCq,Zvalq,Eaq,tmin,Tbulk,Qout);
    end % for I  

    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
    end % intermediates print

    fprintf(fid,'     \n');
    fprintf(fid,'     \n');

    endtxt = [10 'END OF (',UNAME,') QUALITY KINETICS PRINT TO FILE' 10];
    disp(endtxt);


