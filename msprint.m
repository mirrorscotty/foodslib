   function[] = msprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY)
 
%  MSPRINT   Mass balance results printout   
%
%            Synopsis: MATLAB function which prints the mass balance results 
%            table for a unit operation block; called only by unit operation 
%            S-functions.  
% 
%            Call: mprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,
%                          SNARRAY)
%
%            Passed Arguments: 
%             UNAME -  unit operation name 
%             SINAMES - matrix of names of all system inputs 
%             SONAMES - matrix of names of all system outputs 
%             SNNAMES - matrix of names of all system intermediates 
%             SIARRAY - system input array; each column of the array represents 
%              a stream specification which includes: mass flow rate; component 
%              fractions, and component types (Fmass array format) 
%             SOARRAY - system output array; each column of the array represents 
%              a stream specification which includes: mass flow rate; component 
%              fractions, and component types (Emass array format)
%             SNARRAY - system intermediate array; each column of the array 
%              represents a stream specification which includes: mass flow rate; 
%              component fractions, and component types (Imass array format)  

% ******************************************************************************

%   Print results to user specified filename which was determined in rfile.m 
    global filem 
    fid = fopen(filem,'a+');
%   fid = filem 

%   Printable component types 
    Pctype = [comptype('water');comptype('protein');comptype('carbohydrate');comptype('fiber');comptype('fat');comptype('ash');0;0;0;0]; 

%   Limiting or sizing all stream names to 15 characters 
    SINAMES = namesize(SINAMES,15);    
    SONAMES = namesize(SONAMES,15);
    if (size(SNNAMES,1) >= 1)  
      INT = 1;
      SNNAMES = namesize(SNNAMES,15);
    else INT = 0; 
    end  

    fprintf(fid,'     \n');
    fprintf(fid,'%s RESULTS\n',UNAME);
    fprintf(fid,'     \n');
    fprintf(fid,'%s MASS BALANCE RESULTS\n',UNAME);
    fprintf(fid,'     \n');
    fprintf(fid,'Stream             Mass Flow     Total    Water   Protein   Carb.    Fiber     Fat      Ash\n');
    fprintf(fid,'Name                 Rate        Solids\n');
    fprintf(fid,'                    (kg/h)       (w/w)    (w/w)    (w/w)    (w/w)    (w/w)    (w/w)    (w/w)\n');
    fprintf(fid,'----------------------------------------------------------------------------------------------\n');
    fprintf(fid,'SYSTEM INPUTS:\n');
    ninputs = size(SIARRAY,2); 
    for I = 1:ninputs   

      F = SIARRAY(1,I);                % Input stream mass flow rate  
      Fcomp = SIARRAY(2:11,I);         % Input stream component mass fractions (w/w)
      Fctype = SIARRAY(12:21,I);       % Input stream component types  
      Fts = tsolids(Fcomp,Fctype);     % Input stream total solids   

      if (any(abs(Fctype-Pctype) > 0))
        Foldct = Fctype;  
        Foldcomp = Fcomp;
        Fctype = zeros(size(Fctype,1),1);  
        Fcomp = zeros(size(Fctype,1),1);  
        for J = 1:6  
          Fcomp(J) = findcomp(Foldcomp,Foldct,J);
          Fctype(J) = J;  
        end  
      end  

      fprintf(fid,'%s  %11.2f    %6.5f  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f\n', ...
       SINAMES(I,:),F,Fts,Fcomp(1),Fcomp(2),Fcomp(3),Fcomp(4),Fcomp(5),Fcomp(6));
    end   

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM OUTPUTS:\n');
    noutputs = size(SOARRAY,2);     
    for I = 1:noutputs      

      E = SOARRAY(1,I);                % Output stream mass flow rate  
      Ecomp = SOARRAY(2:11,I);         % Output stream component mass fractions (w/w)
      Ectype = SOARRAY(12:21,I);       % Output stream component types 
      Ets = tsolids(Ecomp,Ectype);     % Output stream total solids    

      if (any(abs(Ectype-Pctype) > 0))
        Eoldct = Ectype;  
        Eoldcomp = Ecomp;  
        Ectype = zeros(size(Ectype,1),1);  
        Ecomp = zeros(size(Ectype,1),1);  
        for J = 1:6
          Ecomp(J) = findcomp(Eoldcomp,Eoldct,J);  
          Ectype(J) = J;  
        end
      end  

      if (E > 0) 
      fprintf(fid,'%s  %11.2f    %6.5f  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f\n', ...
       SONAMES(I,:),E,Ets,Ecomp(1),Ecomp(2),Ecomp(3),Ecomp(4),Ecomp(5),Ecomp(6));
      end 
    end  
    fprintf(fid,'----------------------------------------------------------------------------------------------\n');

%   Perform an overall mass balance on the system     
    Fsum = sum(SIARRAY(1,:));          % Sum all system input mass flow rates  
    Esum = sum(SOARRAY(1,:));          % Sum all system output mass flow rates 
    mbal = Fsum - Esum;                % Overall mass balance   

    fprintf(fid,'OVERALL MASS\n');
    fprintf(fid,'BALANCE:         %11.2f\n',mbal);
    fprintf(fid,'----------------------------------------------------------------------------------------------\n');

    if (INT == 1)

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM INTERMEDIATES:\n');
    nints = size(SNARRAY,2);    
    for I = 1:nints      

      N = SNARRAY(1,I);              % Intermediate stream mass flow rate
      Ncomp = SNARRAY(2:11,I);       % Intermediate stream component mass fractions (w/w)
      Nctype = SNARRAY(12:21,I);      % Intermediate stream component types
      Nts = tsolids(Ncomp,Nctype);     % Intermediate stream total solids

      if (any(abs(Nctype-Pctype) > 0))
        Noldct = Nctype;
        Noldcomp = Ncomp;
        Nctype = zeros(size(Nctype,1),1);
        Ncomp = zeros(size(Nctype,1),1);
        for J = 1:6
          Ncomp(J) = findcomp(Noldcomp,Noldct,J);
          Nctype(J) = J;
        end
      end  
 
      fprintf(fid,'%s  %11.2f    %6.5f  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f\n', ...
       SNNAMES(I,:),N,Nts,Ncomp(1),Ncomp(2),Ncomp(3),Ncomp(4),Ncomp(5),Ncomp(6));
    end 
    fprintf(fid,'----------------------------------------------------------------------------------------------\n');
 
    end % intermediates print 


    fprintf(fid,'     \n');
    fprintf(fid,'     \n');

    endtxt = [10 'END OF (',UNAME,') MASS BALANCE PRINT TO FILE' 10];
    disp(endtxt);


