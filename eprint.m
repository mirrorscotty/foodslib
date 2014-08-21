
   function[] = eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY)  

%  EPRINT    Energy balance results printout    
%
%            Synopsis: MATLAB function which prints the energy balance results 
%            table for a unit operation block; called only by unit operation 
%            S-functions.  
%
%            Call: eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY) 
%
%            Passed Arguments: 
%             UNAME -  unit operation name 
%             SINAMES - matrix of names of all system inputs 
%             SONAMES - matrix of names of all system outputs
%             SNNAMES - matrix of names of all system intermediates
%             SIARRAY - system input array; each column of the array represents 
%              a stream specification which includes: mass flow rate; component 
%              fractions, component types, temperature, pressure, and quality 
%              (Fme array format) 
%             SOARRAY - system output array; each column of the array represents 
%              a stream specification which includes: mass flow rate; component 
%              fractions, component types, temperature, pressure, and quality 
%              (Eme array format) 
%             SNARRAY - system intermediate array; each column of the array 
%              represents a stream specification which includes: mass flow rate; 
%              component fractions, component types, temperature, pressure, 
%              and quality (Nme array format)

% **********************************************************************************

%   Print results to user specified filename which was determined in rfile.m 
    global filem;
    % THIS IS WHERE THE CHANGE WAS MADE 4/21/01
    fid = fopen(filem,'a+');
    %fid = filem;
    % END CHANGE 4/21/01

%   Limiting or sizing all stream names to 15 characters 
    SINAMES = namesize(SINAMES,15);   
    SONAMES = namesize(SONAMES,15);   
    if (size(SNNAMES,1) >=1)
      INT = 1;
      SNNAMES = namesize(SNNAMES,15);
    else INT = 0;
    end

%   Reference temperature (K)
    Tref = tref;

    fprintf(fid,'     \n');
    fprintf(fid,'%s ENERGY BALANCE RESULTS\n',UNAME);
    fprintf(fid,'     \n');
    fprintf(fid,'Stream                   Total       Temp.      Press.     Quality       Heat       Enthalpy\n');
 
    fprintf(fid,'Name                    Enthalpy                                       Capacity\n');
    fprintf(fid,'                         (kJ/h)       (C)       (MPa)                  (kJ/kg-K)    (KJ/kg)\n');
 
    fprintf(fid,'----------------------------------------------------------------------------------------------\n');
    fprintf(fid,'SYSTEM INPUTS:\n');
    ninputs = size(SIARRAY,2); 
    for I = 1:ninputs   

      F = SIARRAY(1,I);                % Input stream mass flow rate (kg/h)  
      Fcomp = SIARRAY(2:11,I);         % Input stream component mass fractions (w/w)
      Fctype = SIARRAY(12:21,I);       % Input stream component types  
      flagmc = purechk(Fcomp,Fctype);  % Flag for pure water stream (=1) 
      FtempK = SIARRAY(22,I);          % Input stream temperature (K) 
      FtempC = FtempK - Tref;          % Input stream temperature (C)  
      Fpres = SIARRAY(23,I);           % Input stream pressure (MPa) 
      Fx = SIARRAY(24,I);              % Input stream quality  
      Fh = enthalpy(SIARRAY(:,I));     % Input stream enthalpy (kJ/kg)      
      Fhtot(I) = F*Fh;                    % Input stream total enthalpy (kJ/h)  

      if (flagmc == 0)                 % Multi-component stream...  
        Fcp = cp(Fcomp,Fctype,FtempK); %   Heat capacity (kJ/kg-K)   
  
        fprintf(fid,'%s  %15.2f    %6.2f       ---         ---        %6.3f     %8.2f\n',...
        SINAMES(I,:),Fhtot(I),FtempC,Fcp,Fh);
      elseif (flagmc == 1)               % Pure water stream       
        fprintf(fid,'%s  %15.2f    %6.2f     %6.4f       %4.3f         ---      %8.2f\n',...
        SINAMES(I,:),Fhtot(I),FtempC,Fpres,Fx,Fh);
      elseif (flagmc == 2)               % Pure water stream
        fprintf(fid,'%s  %15.2f    %6.2f     %6.4f         ---         ---      %8.2f\n',...
        SINAMES(I,:),Fhtot(I),FtempC,Fpres,Fh);
      end  

    end   

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM OUTPUTS:\n');
    noutputs = size(SOARRAY,2);     
    for I = 1:noutputs      

      E = SOARRAY(1,I);                % Output stream mass flow rate (kg/h)  
      Ecomp = SOARRAY(2:11,I);         % Output stream component mass fractions (w/w)
      Ectype = SOARRAY(12:21,I);       % Output stream component types
      flagmc = purechk(Ecomp,Ectype);  % Flag for pure water stream (=1)
      EtempK = SOARRAY(22,I);          % Output stream temperature (K)
      EtempC = EtempK - Tref;          % Output stream temperature (C)
      Epres = SOARRAY(23,I);           % Output stream pressure (MPa)
      Ex = SOARRAY(24,I);              % Output stream quality
      Eh = enthalpy(SOARRAY(:,I));     % Output stream enthalpy (kJ/kg)
      Ehtot(I) = E*Eh;                 % Output stream total enthalpy (kJ/h)
 
      if (flagmc == 0)                 % Multi-component stream...
        Ecp = cp(Ecomp,Ectype,EtempK); %   Heat capacity (kJ/kg-K)
 
        fprintf(fid,'%s  %15.2f    %6.2f       ---         ---        %6.3f     %8.2f\n',...
        SONAMES(I,:),Ehtot(I),EtempC,Ecp,Eh);
      elseif (flagmc == 1)               % Pure water stream
        fprintf(fid,'%s  %15.2f    %6.2f     %6.4f       %4.3f         ---      %8.2f\n',...
        SONAMES(I,:),Ehtot(I),EtempC,Epres,Ex,Eh);
      elseif (flagmc == 2)               % Pure water stream
        fprintf(fid,'%s  %15.2f    %6.2f     %6.4f         ---         ---      %8.2f\n',...
        SONAMES(I,:),Ehtot(I),EtempC,Epres,Eh);
      end

    end  
    fprintf(fid,'----------------------------------------------------------------------------------------------\n');

%   Perform an overall energy balance on the system     
    Fhtotsum = sum(Fhtot);          % Sum all system input total enthalpies  
    Ehtotsum = sum(Ehtot);          % Sum all system output total enthalpies 
    ebal = Fhtotsum - Ehtotsum;     % Overall energy balance   

    fprintf(fid,'OVERALL ENERGY\n');
    if (ebal > 0.01)
      fprintf(fid,'BALANCE:         %15.2f (Heat Gain By System)\n',ebal);
    elseif (ebal < -0.01) 
      fprintf(fid,'BALANCE:         %15.2f (Heat Loss From System)\n',ebal);
    else  
      fprintf(fid,'BALANCE:         %15.2f\n',ebal);
    end 

    fprintf(fid,'----------------------------------------------------------------------------------------------\n');

    if (INT == 1)

    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM INTERMEDIATES:\n');
    nints = size(SNARRAY,2);
    for I = 1:nints
 
      N = SNARRAY(1,I);                % Intermed. stream mass flow rate (kg/h)
      Ncomp = SNARRAY(2:11,I);         % Intermed. stream component mass fractions (w/w)
      Nctype = SNARRAY(12:21,I);       % Intermed. stream component types
      flagmc = purechk(Ncomp,Nctype);  % Flag for pure water stream (=1)
      NtempK = SNARRAY(22,I);          % Intermed. stream temperature (K)
      NtempC = NtempK - Tref;          % Intermed. stream temperature (C)
      Npres = SNARRAY(23,I);           % Intermed. stream pressure (MPa)
      Nx = SNARRAY(24,I);              % Intermed. stream quality
      Nh = enthalpy(SNARRAY(:,I));     % Intermed. stream enthalpy (kJ/kg)
      Nhtot(I) = N*Nh;                 % Intermed. stream total enthalpy (kJ/h)
 
      if (flagmc == 0)                 % Multi-component stream...
        Ncp = cp(Ncomp,Nctype,NtempK); %   Heat capacity (kJ/kg-K)
 
        fprintf(fid,'%s  %15.2f    %6.2f       ---         ---        %6.3f     %8.2f\n',...
        SNNAMES(I,:),Nhtot(I),NtempC,Ncp,Nh);
      elseif (flagmc == 1)               % Pure water stream
        fprintf(fid,'%s  %15.2f    %6.2f     %6.4f       %4.3f         ---      %8.2f\n',...
        SNNAMES(I,:),Nhtot(I),NtempC,Npres,Nx,Nh);
      elseif (flagmc == 2)               % Pure water stream
        fprintf(fid,'%s  %15.2f    %6.2f     %6.4f         ---         ---      %8.2f\n',...
        SNNAMES(I,:),Nhtot(I),NtempC,Npres,Nh);
      end
 
    end
    fprintf(fid,'----------------------------------------------------------------------------------------------\n');

    end % intermediates print


    fprintf(fid,'     \n');
    fprintf(fid,'     \n');

    endtxt = [10 'END OF (',UNAME,') ENERGY BALANCE PRINT TO FILE' 10]; 
    disp(endtxt);
 


