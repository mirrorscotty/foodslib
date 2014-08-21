      function [sys,x0] = ufs(t,x,u,flag,Uname,Fname,Rname,Pname)

%  UFS       Ultrafiltration model   
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for ultrafiltration steady-state 
%            model. Determines mass and energy balances. 
%   
%            Call: splitms(t,x,u,flag,Uname,Fname,Rname,Pname)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Uname - mass splitter name
%             Fname - feed stream name 
%             Rname - retentate stream name 
%             Pname - permeate stream name 
%  
%            Returned Arguments: 
%             sys - SIMULINK S-function variable: for flag = 0, 
%                   sys = sizes of parameters and initial conditions. 
%                   For flag = 3, sys = system outputs which are the 
%                   u-arrays for each stream exiting the unit operation 
% 

%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = 1;                            % Number of system input streams
SO = 2;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = 49;                     % Number of ultrafiltration inputs 

if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, (SO*usize), ((SI*usize)+UOinputs), 0, 1];
   x0 = ones(size(u));

elseif flag == 3 | flag == 9
   % if flag=3, return system outputs
  
%     Breakup of SIMULINK constructed u-array into stream arrays and 
%      operating information array  
 
      Fuarray = u(1:usize);          % full feed stream  array 
      UOuarray = u(SI*usize+1:SI*usize+UOinputs); 
                                     % full operating inputs array

%  *********************** USER-DEFINED DESIGN LEVEL ************************

      dlevel = UOuarray(1);           % design level
      mandq = UOuarray(2);            % micro and quality assessment
      econ = UOuarray(3);             % economic evaluation  

%  ******************* MASS BALANCE LEVEL VARIABLE ASSIGNMENTS **************

%     Feed stream information  
      F = Fuarray(ulocate('mass')); % mass flow rate (kg/h)
      Fcomp = Fuarray(ulocate('comp1st'):ulocate('complast'));
                                    % component mass fractions (w/w)
      Fctype = Fuarray(ulocate('ctype1st'):ulocate('ctypelast')); 
                                    % component types  
      Fmass = [F;Fcomp;Fctype];     % mass information array 

%     Ultrafiltration information for mass balances  
      Rts = UOuarray(4);            % retentate total solids 
      for I = 1:9
        RJctype(I) = UOuarray(3+2*I);  % rejection components   
        RJrat(I) = UOuarray(4+2*I);    % rejection ratio by component  
      end  

%  ***************************  MASS BALANCE LEVEL ************************** 

%     Total solids of feed (w/w)
      Fts = tsolids(Fcomp,Fctype);

%     User error message:
      if (Rts < Fts)  % Retentate total solids less than feed solids 
        errortxt = [10 'ULTRAFILTRATION (',Uname,') ERROR: Specified total solids of' 10 ...
        'retentate is less than total solids of feed stream. ' 10];
        disp(errortxt);
        return;        % Stop execution of algorithm
      end
 
%     Given the component rejection ratios and the retentate total solids;
%     a system of non-linear equations must be solved to obtain the
%     retentate composition and flow rate. Newton's Method is used.
 
%     Find all non-zero component mass fraction indexes in feed
      nonzero = find(Fcomp);  
%     Number of non-zero component mass fractions
      Nnonzero = nnz(Fcomp);  
%     Set new composition array
      J = 0;
      for I = 1:Nnonzero
        if (Fctype(nonzero(I)) ~= comptype('water'))
          J = J+1;
          Ccomp(J) = Fcomp(nonzero(I));  
          Cctype(J) = Fctype(nonzero(I));  
        else 
          Nnonzero = Nnonzero-1; 
        end 
      end 

%     Place rejection ratios in same order as components 
      for I = 1:Nnonzero
        J = findval1(RJctype,Cctype(I));  
        if ~isempty(J)
          reject(I) = RJrat(J);
        end  
      end  

      if (Nnonzero > 1)
%       Initial guess at weight concentration ratio 
%        (WCR = feed rate/retentate rate
        unknown(Nnonzero+1) = Rts/Fts;  
 
%       Initial guesses at retentate composition based on guess of weight
%        concentration ratio
        for I = 1:Nnonzero
          unknown(I) = Ccomp(I)*unknown(Nnonzero+1)^reject(I);
        end
 
%       Initialize function array
        f = ones(size(1,Nnonzero+1));
%       Convergence criteria; WCR yields desired retentate solids
        while (any(abs(f)>1e-5) == 1)

%         Evaluate function: f(unknown)
          f(1:Nnonzero) = unknown(1:Nnonzero) - ... 
            Ccomp(1:Nnonzero) .* unknown(Nnonzero+1) .^reject(1:Nnonzero);
          f(Nnonzero+1) = sum(unknown(1:Nnonzero)) - Rts;
          B = -1 .*f;

%         Evaluate Jacobian matrix; filled with partial derivatives of f function
          Apd = zeros(Nnonzero+1);
          for I = 1:Nnonzero
            Apd(I,I) = 1;
            Apd(Nnonzero+1,I) = 1;
          end
          for I = 1:Nnonzero
            Apd(I,Nnonzero+1) = -reject(I) .*Ccomp(I) .*(unknown(I)) .^(reject(I)-1);
          end
          Apd(Nnonzero+1,Nnonzero+1) = 0;

%         Solve linear equation; D values are change from base point 
%          values of unknowns
          D = (inv(Apd)*B')';
          for I = 1:Nnonzero
            if (abs(unknown(I)) > 1e-10) & (abs(D(I)/unknown(I)) > 0.5)
              D(I) = 0.5*unknown(I)*abs(D(I))/D(I);
            end
          end
 
%         Improved estimate of retentate composition and WCR
          unknown = unknown + D;
        end

      elseif (Nnonzero == 1)
        unknown(1) = Rts;  
        unknown(Nnonzero+1) = (unknown(1)/Fts)^(1/reject(nonzero(1)));  
      end

%     Final value of WCR
      WCR = unknown(Nnonzero+1);  
       
%     Retentate composition (w/w) and flow rate (kg/h)
      Rcomp = unknown(1:Nnonzero);  
      Rctype = Cctype;
      Rcomp(Nnonzero+1) = 1 - sum(Rcomp);  
      Rctype(Nnonzero+1) = comptype('water');  
      if (length(Rcomp) < compmax) 
        Rcomp(compmax) = 0;
        Rctype(compmax) = 0;
      end  
      R = F/unknown(Nnonzero+1);  
       
%     Permeate flow rate (kg/h)
      P = F - R;
       
%     Permeate composition (w/w)
      for I = 1:Nnonzero 
%       Pcomp(I) = abs((F*Ccomp(nonzero(I)) - R*Rcomp(nonzero(I)))/P)
        Pcomp(I) = abs((F*Ccomp(I) - R*Rcomp(I))/P);  
      end  
      Pcomp(Nnonzero+1) = 1 - sum(Pcomp);
      if (length(Pcomp) < compmax) 
        Pcomp(compmax) = 0;
      end  
      Pctype = Rctype; 
       
%     Exit stream mass arrays; used for print results table  
      Rmass = [R;Rcomp';Rctype'];  
      Pmass = [P;Pcomp';Pctype'];

%     System output arrays 
      Ruarray = Rmass;  
      Puarray = Pmass;  

      endtxt = [10 'END OF ULTRAFILTRATION (',Uname,') MASS BALANCE' 10];  
      disp(endtxt);


%  **************** ENERGY BALANCE LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 2)  

%     Reference temperature (K)
      Tref = tref;

%     Feed stream information  
      FtempC = Fuarray(ulocate('temp'));% temperature (C) 
      FtempK = FtempC + Tref;            % temperature (K)
      Fpres = Fuarray(ulocate('pressure'));  % pressure (MPa)
      Fx = Fuarray(ulocate('quality'));    % quality    
      Fme = [Fmass;FtempK;Fpres;Fx]; 
                                         % mass and energy information array  

%  *********************** ENERGY BALANCE LEVEL *****************************

%     Retentate and Permeate temperatures (K)
      RtempK = FtempK;
      PtempK = FtempK;
     
%     Retentate and permeate pressure (MPa) 
      Rpres = Fpres;
      Ppres = Fpres;  

%     Retentate and permeate quality 
      Rx = Fx;
      Px = Fx;  

%     Update mass and energy info arrays 
      Rme = [Rmass;RtempK;Rpres;Rx]; 
      Pme = [Pmass;PtempK;Ppres;Px]; 

%     System output arrays
      Ruarray = Rme;
      Puarray = Pme;

      endtxt = [10 'END OF ULTRAFILTRATION (',Uname,') ENERGY BALANCE' 10];
      disp(endtxt);

      end % Energy Balance Level  

%  *************** BASIC TRANSPORT LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 3)

%     Unit operation information for basic transport (sizing)
      flux = UOuarray(23);          % flux through membrane (L/h-m^2)

      end % Basic Transport Level Assignments

%  ************************** BASIC TRANSPORT LEVEL *************************

      if (dlevel >= 3)

%     Density of permeate (kg/m^3)
      Frho = rho(Fcomp,Fctype,FtempK);  

%     Density of retentate (kg/m^3)
      Rrho = rho(Rcomp,Rctype,RtempK);
 
%     Density of permeate (kg/m^3)
      Prho = rho(Pcomp,Pctype,PtempK);  

%     Membrane surface area (m^2) - Note: flux (L/m^2-h)
      Area = P*1000/(Prho*flux);  

      endtxt = [10 'END OF ULTRAFILTRATION (',Uname,') BASIC TRANSPORT' 10];
      disp(endtxt);

      end  % Basic Transport Level

%  ************* ADVANCED TRANSPORT LEVEL VARIABLE ASSIGNMENTS **************

      if (dlevel >= 4)

%     Unit operation inforamtion for advanced transport
      mumodel = UOuarray(24);     % viscosity model
      fluxmodl = UOuarray(25);    % flux model 

%     Hagen-Poiseuille law model inputs 
      Ptrans = UOuarray(26);      % transmembrane pressure (MPa)
      Dpore = UOuarray(27);       % mean pore diameter (microns)
      delx = UOuarray(28);        % skin thickness (microns)
      Npores = UOuarray(29);           % pore density (pores/cm^2)

%     Mass transfer (film theory) model inputs 
      Dh = UOuarray(30);          % hydraulic (tube/fiber) diameter (m)
      Tlength = UOuarray(31);     % tube length (m)
      Ntubes = UOuarray(32);      % number of tubes/fibers
      diff = UOuarray(33);        % diffusivity (cm^2/s)
      Cg = UOuarray(34);          % gel layer concentration (w/w)
      default = UOuarray(35);     % default to Grober, Leveque, or
                                    %  Coburn-Chilton mass transfer models
      A = UOuarray(36);           % mass transfer constant
      alpha = UOuarray(37);       % mass transfer Reynold's number exponent
      beta = UOuarray(38);        % mass transfer Schmidt number exponent
      w = UOuarray(39);           % mass transfer boundary layer 
                                    % development exponent
 
%  ************************ ADVANCED TRANSPORT LEVEL ************************

%     Power law values for viscosity at bulk temperature
      powerlaw = viscos(Fcomp,Fctype,FtempK,mumodel);
      Kbulk = powerlaw(1);       % Consistency Coefficient (Pa-s^n)
      nbulk = powerlaw(2);       % Flow behavior index

      viscbulk = Kbulk;

      if (fluxmodl == 1) 
%   -----------------------------------------------------------------------
%     Hagen-Poiseuille law for streamline flow through channels - predicting
%      flux under pressure dependent conditions

%     Porosity (percent of membrane surface occupied by pores)
%     Pore diameter (microns --> cm)
      eta = Npores*pi*(Dpore*1e-4)^2/4;  

%     Flux (L/m^2-h); time (s --> h); (m^3 --> L)
      fluxpd = eta*Dpore^2*Ptrans*1000*3600/(32*delx*viscbulk);  

%     Membrane surface area (m^2) - Note: flux (L/m^2-h)
      Area = P*1000/(Prho*fluxpd);  

      elseif (fluxmodl == 2) 
%   -----------------------------------------------------------------------
%     Mass transfer (film theory) model for predicting flux under pressure
%      independent conditions

%     Velocity (m/s)
      vel = F*4/(Ntubes*Frho*pi*Dh^2*3600);  
 
%     Reynolds Number  
      NRe = reynolds(Dh,Dh,vel,Frho,Kbulk,nbulk,0); 
 
%     Schmidt Number
      NSc = viscbulk*100^2/(Frho*diff);  
 
%     For default models
      if (default == 1)
        if NRe > 4000
%         Chilton-Colburn membrane system
          A = 0.023;  
          alpha = 0.8;  
          beta = 0.33;  
          w = 0;  
        elseif NRe < 1800
%         Channel length to attain fully developed flow (m)
          Lv = 0.0575*NRe*Dh;  
          if (Lv > Tlength)
%           Grober membrane system
            A = 0.664;  
            alpha = 0.5;  
            beta = 0.33;  
            w = 0.5;  
          else
%           Leveque membrane system
            A = 1.86;  
            alpha = 0.33;  
            beta = 0.33;
            w = 0.33;  
          end
        end
      end
 
%     Sherwood Number
      NSh = A*(NRe^alpha)*(NSc^beta)*(Dh/Tlength)^w;  
 
%     Mass transfer coefficient (cm^2/m-s --> L/m^2-h)
      k = NSh*diff*1000*3600/(100^2*Dh);  
 
%     Flux (L/m^2-h)
      fluxpid = k*log(Cg/Fts);  
 
%     Membrane surface area (m^2) - Note: flux (L/m^2-h)
      Area = P*1000/(Prho*fluxpid);  

      end  
 
      endtxt = [10 'END OF ULTRAFILTRATION (',Uname,') ADVANCED TRANSPORT' 10];
      disp(endtxt);

      end  % Advanced Transport Level

%  *************** MICROBIAL & QUALITY VARIABLE ASSIGNMENTS *****************

%     Feed Microbial Death Kinetics
      Dvalm = Fuarray(ulocate('microDval'));   % D value (min)
      DtempCm =  Fuarray(ulocate('microDtemp'));% D value temperature basis (C)
      DtempKm = DtempCm+Tref; % D value temperature basis (K)
      Zvalm =  Fuarray(ulocate('microZval'));   % Z value (C)
      Eam =  Fuarray(ulocate('microEa'));       % activation energy (kcal/mol)
      Noin =  Fuarray(ulocate('microload'));     % microbial load
      Fmdk = [Dvalm;DtempCm;Zvalm;Eam;Noin];  % microbial array 

%     Feed Quality Kinetics
      Dvalq = Fuarray(ulocate('qualDval'));    % D value (min)
      DtempCq = Fuarray(ulocate('qualDtemp')); % D value temperature basis (C)
      DtempKq = DtempCq+Tref; % D value temperature basis (K)
      Zvalq = Fuarray(ulocate('qualZval'));     % Z value (C)
      Eaq = Fuarray(ulocate('qualEa'));       % activation energy (kcal/mole)
      Qoin =  Fuarray(ulocate('qualload'));     % quality load
      Fqual = [Dvalq;DtempCq;Zvalq;Eaq;Qoin];  % quality array 
 
%  ********************* MICROBIAL & QUALITY ASSESSMENT *********************

%     Pass microbial and quality assessment variables to exit streams 
      Rmdk = Fmdk;
      Pmdk = zeros(size(Fmdk,1));
      Rqual = Fqual;
      Pqual = zeros(size(Fqual,1));

      if (mandq == 1)   

%     Bulk temperature for microbial assessment (K)
      Tbulk = FtempK;   

%     Minimum Residence Time (min)
      tmin = 0;   
 
%     Final microbial load 
      Nout = mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,Tbulk,tmin);  
 
%     Log reduction in microbial load
      DlogNred = log10(Noin/Nout);

%     Temperature for quality assessment (K)
      Tqual = FtempK;

%     Maximum processing time = Minimum processing time (min) 
      tqmax = 0;

%     Final quality load  
      Qout =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,tqmax); 

%     Log reduction in quality 
      DlogQred = log10(Qoin/Qout);

%     Exit stream microbial and quality arrays
      Rmdk = Fmdk;  
      Rmdk(size(Fmdk,1)) = Nout;
      Rqual = Fqual;
      Rqual(size(Fqual,1)) = Qout;

      endtxt = [10 'END OF ULTRAFILTRATION (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

      end  % microbial and quality assessment 

%     Update system output array
%     Ruarray = [Rme;Rmdk;Rqual];
%     Puarray = [Pme;Pmdk;Pqual];

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(40:49); % economic analysis array 
      end  

%  ********************* SYSTEM OUTPUT ARRAY ASSIGNMENT *********************

%     Size stream info arrays to full system output u-arrays  
      Ruarray = ustd(Ruarray);
      Puarray = ustd(Puarray);

      sys =[Ruarray;Puarray;];    
 
%  ****************************** PRINTED OUTPUT ****************************

  if flag == 9    % print only after a solution is found 

%   Print results to specified file name; 
%    call filename manager function rfile.m 
    fid = rfile;

%   Name of unit operation; system input streams; system output streams  
    UNAME = Uname; 
    SINAMES = str2mat(Fname);
    SONAMES = str2mat(Rname,Pname);

%   System inputs array 
    SIARRAY = [Fmass]; 
%   System outputs array 
    SOARRAY = [Rmass Pmass];

%   No intermediate streams
    SNNAMES = [];
    SNARRAY = []; 

%   Print mass balance results table
    msprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    if (dlevel >= 2) 
%     System inputs array  
      SIARRAY = [Fme]; 
%     System outputs array 
      SOARRAY = [Rme Pme]; 

%     Print energy balance results table
      eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    end  

    fprintf(fid,'ULTRAFILTRATION MODEL INPUTS\n');     
    fprintf(fid,'     \n');
    fprintf(fid,' Weight concentration ratio \n');
    fprintf(fid,'  (WCR = feed rate/retentate rate)      %8.5e \n',WCR);
    for I = 1:length(RJctype) 
    index = RJctype(I);
    if (index > 0)  
    fprintf(fid,' %s Rejection Ratio                     %8.5f \n',compname(index),RJrat(I));
    end  
    end
    fprintf(fid,'     \n');

    if (dlevel == 3) 
    fprintf(fid,'  Flux                                 %8.3f L/m^2-h\n',flux);
    end
    fprintf(fid,'     \n');

    if (dlevel == 4) 
      if (fluxmodl == 1)  
    fprintf(fid,' Hagen-Poiseuille Model - Pressure Dependent:\n');
    fprintf(fid,'  Transmembrane Pressure                %8.5f MPa\n',Ptrans);
    fprintf(fid,'  Mean Pore Diameter                    %8.5f microns\n',Dpore);
    fprintf(fid,'  Skin Thickness                        %8.5f microns\n',delx);
    fprintf(fid,'  Pore Density                           %8.4e pores/cm^2\n',Npores);
    fprintf(fid,'  Porosity                             %8.4f \n',eta);

    fprintf(fid,'     \n');
      elseif (fluxmodl == 2) 
    fprintf(fid,' Mass Transfer (Thin Film) Model - Pressure Independent:\n');
    fprintf(fid,'  Hydraulic (tube/fiber) Diameter       %8.5f m\n',Dh);
    fprintf(fid,'  Tube/Fiber Length                     %8.5f m\n',Tlength);
    fprintf(fid,'  Number of Tubes/Fibers                 %5.0f \n',Ntubes);
    fprintf(fid,'  Diffusivity                            %8.5e cm^2/s\n',diff);
    fprintf(fid,'  Gel Layer Concentration               %8.5f w/w\n',Cg);
    fprintf(fid,'  Model Coefficient and Exponents\n');
    fprintf(fid,'   A - Coefficient                      %8.5f w/w\n',A);
    fprintf(fid,'   alpha - NRe exponent                 %8.5f w/w\n',alpha);
    fprintf(fid,'   beta - NSc exponent                  %8.5f w/w\n',beta);
    fprintf(fid,'   w - (diameter/length) exponent       %8.5f w/w\n',w);
      end 
    end 

    if (dlevel >= 3)
    fprintf(fid,'     \n');
    fprintf(fid,'ULTRAFILTRATION MODEL OUTPUTS\n');
    fprintf(fid,'     \n');

    fprintf(fid,'  Feed Density                         %8.2f kg/m^3\n',Frho); 
    fprintf(fid,'  Retentate Density                    %8.2f kg/m^3\n',Rrho); 
    fprintf(fid,'  Permeate Density                     %8.2f kg/m^3\n',Prho); 
    fprintf(fid,'  Feed Volumetric Flow                 %8.2f L/h\n',F/Frho*1000); 
    fprintf(fid,'  Retentate Volumetric Flow            %8.2f L/h\n',R/Rrho*1000); 
    fprintf(fid,'  Permeate Volumetric Flow             %8.2f L/h\n',P/Prho*1000); 
    end 
    if (dlevel == 3)
    fprintf(fid,'     \n');
    fprintf(fid,'  Membrane Surface Area                %8.2f m^2\n',Area); 
    end

    if (dlevel == 4)
    fprintf(fid,'  Bulk Consistency Coefficient     %8.3e Pa-s^n\n',Kbulk);
    fprintf(fid,'  Bulk Flow Behavior Index (n)      %8.3f \n',nbulk);

    fprintf(fid,'     \n');

      if (fluxmodl == 1) 
    fprintf(fid,'  Flux                                  %8.3f L/m^2-h\n',fluxpd);                         
    fprintf(fid,'  Membrane Surface Area                %8.2f m^2\n',Area);
    fprintf(fid,'     \n');

      elseif (fluxmodl == 2) 

    fprintf(fid,'  Fluid velocity                        %8.1f m/s\n',vel);
    fprintf(fid,'  Reynolds Number (NRe)                 %8.1f \n',NRe);
    fprintf(fid,'  Schmidt Number (NSc)                  %8.1f \n',NSc);
    fprintf(fid,'  Sherwood Number (NSc)                 %8.1f \n',NSh);
    fprintf(fid,'  Mass Transfer Coefficient            %8.3f L/m^2-h\n',k); 
    fprintf(fid,'  Flux                                 %8.3f L/m^2-h\n',fluxpid);
    fprintf(fid,'  Membrane Surface Area                %8.2f m^2\n',Area);
      end  
    fprintf(fid,'     \n');

    endtxt = [10 'END OF (',Uname,') TRANSPORT RESULTS PRINT TO FILE' 10];
    disp(endtxt); 
 
    end 

%   Microbial and Quality Assessment results printout 
    if (mandq == 1) & (dlevel >= 2)

    fprintf(fid,'     \n');
    fprintf(fid,'  Microbial (Worst Case) Assessment:   \n');
    fprintf(fid,'   Minimum Residence Time                %5.2f min \n',tmin);
    fprintf(fid,'   Minimum Processing Temperature        %5.2f C \n',Tbulk-Tref);
    fprintf(fid,'   Log Reduction in Microbial Load       %5.2f \n',DlogNred);
    fprintf(fid,'     \n'); 
    fprintf(fid,'  Quality (Worst Case) Assessment:   \n');
    fprintf(fid,'   Maximum Residence Time                %5.2f min \n',tqmax);
    fprintf(fid,'   Maximum Processing Temperature        %5.2f C \n',Tqual-Tref);  
    fprintf(fid,'   Log Reduction in Quality Load         %5.2f \n',DlogQred);

%     System input and exit stream names for microbial and quality assessment 
      SIMQNAME = str2mat(Fname);
      SOMQNAME = str2mat(Rname,Pname); 

%     System input microbial array
      SIMARRAY = [Fmdk];
%     System outputs microbial array
      SOMARRAY = [Rmdk Pmdk]; 

%     No intermediate streams
      SNMQNAME = [];
      SNMARRAY = [];
      SNQARRAY = [];

      mdkprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIMARRAY,SOMARRAY,SNMARRAY);

%     System input quality array
      SIQARRAY = [Fqual];
%     System outputs quality array
      SOQARRAY = [Rqual Pqual];  

      qprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIQARRAY,SOQARRAY,SNQARRAY);
    end 

    if (econ >= 1) 

      Uutil = [0 0 0 0 0 0]; 
      Usize = Area;
      Unum = 1; 
      Uecon = ECON;  
      uofile(Uutil,Usize,Unum,Uecon,Uname);

    end  


    fprintf(fid,'     \n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'     \n');

    sys = []; 
  end % flag 9 print 

else
   % Otherwise, no need to return anything since this is continuous system
   sys = [];

end
 
