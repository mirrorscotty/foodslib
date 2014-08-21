      function [sys,x0] = pumps(t,x,u,flag,Uname,Fname,Ename)  

%  PUMPS     Pump model   
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for pump steady-state model. 
%            Determine mass and energy balances; performs design at 
%            basic and advanced transport level. 
%   
%            Call: pumps(t,x,u,flag,Uname,Fname,Ename)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Uname - fitting name
%             Fname - feed stream name 
%             Ename - exit stream name 
%  
%            Returned Arguments: 
%             sys - SIMULINK S-function variable: for flag = 0, 
%                   sys = sizes of parameters and initial conditions. 
%                   For flag = 3, sys = system outputs which are the 
%                   u-arrays for each stream exiting the unit operation 
%
%            References: 
%             Geankoplis, C.J. 1983.  "Transport Processes and Unit
%               Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
% 

%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = 1;                            % Number of system input streams
SO = 1;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = 20;                     % Number of pump inputs 

if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, (SO*usize), ((SI*usize)+UOinputs), 0, 1];
   x0 = ones(size(u));

elseif flag == 3 | flag == 9
   % if flag=3, return system outputs
  
%     Breakup of SIMULINK constructed u-array into stream arrays and 
%      operating information array  
 
      Fuarray = u(1:usize);          % full feed stream array 
      UOuarray = u(SI*usize+1:SI*usize+UOinputs); 
                                     % full operating inputs array

%  *********************** USER-DEFINED DESIGN LEVEL ************************

      dlevel = UOuarray(1);           % design level
      mandq = UOuarray(2);            % micro and quality assessment
      econ = UOuarray(3);             % economic evaluation  

%  ******************* MASS BALANCE LEVEL VARIABLE ASSIGNMENTS **************

%     Feed stream 1 information  
      F = Fuarray(ulocate('mass')); % mass flow rate (kg/h)
      Fcomp = Fuarray(ulocate('comp1st'):ulocate('complast'));
                                    % component mass fractions (w/w)
      Fctype = Fuarray(ulocate('ctype1st'):ulocate('ctypelast')); 
                                    % component types  
      Fmass = [F;Fcomp;Fctype];     % mass information array 

%  ***************************  MASS BALANCE LEVEL ************************** 

%     Overall Mass Balance (kg/h) 
      E = F; 

%     Component balances (w/w) and types 
      Ecomp = Fcomp;
      Ectype = Fctype;  

%     Exit stream mass arrays; used for print results table  
      Emass = [E;Ecomp;Ectype];

%     System output arrays 
      Euarray = Emass;  

      endtxt = [10 'END OF PUMP (',Uname,') MASS BALANCE' 10];  
      disp(endtxt); 

%  **************** ENERGY BALANCE LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 2)  

%     Reference temperature (K)
      Tref = tref;

%     Feed stream information  
      FtempC = Fuarray(ulocate('temp'));% temperature (C) 
      FtempK = FtempC + Tref;            % temperature (K) 
      Fpres = Fuarray(ulocate('pressure'));  
      Fx = Fuarray(ulocate('quality')); 
      Fme = [Fmass;FtempK;Fpres;Fx]; 
                                         % mass and energy information array  

%  *********************** ENERGY BALANCE LEVEL *****************************

%     Exit stream information 
      EtempK = FtempK;
      Epres = Fpres; 
      Ex = Fx;

%     Update mass and energy info arrays 
      Eme = [Emass;EtempK;Epres;Ex]; 

%     System output arrays
      Euarray = Eme;

      endtxt = [10 'END OF PUMP (',Uname,') ENERGY BALANCE' 10];
      disp(endtxt);

      end % Energy Balance Level  

%  ************* BASIC TRANSPORT LEVEL VARIABLE ASSIGNMENTS **************

      if (dlevel >= 3) 

%     Pump inforamtion for advanced transport  
      diamin = UOuarray(4);     % inlet diameter (m)
      diamout = UOuarray(5);    % outlet diameter (m)
      delP = UOuarray(6);      % pressure increase (MPa)
      delH = UOuarray(7);       % height differential (m)
      eta = UOuarray(8);        % pump efficiency
      etae = UOuarray(9);       % motor efficiency
      mumodel = UOuarray(10);   % viscosity model 

%  ************************ BASIC TRANSPORT LEVEL ************************

%     Bulk temperature (K)
      Tbulk = (FtempK+EtempK)/2;  

%     Bulk density (kg/m^3)
      FEme = Fme;
      xtemp = ulocate('temp'); 
      FEme(xtemp) = Tbulk;
      rhobulk = density(FEme);

%     Inlet cross sectional area of tube (m^2)
      Axin = pi*(diamin/2)^2;  
%     Volumetric flowrate (m^3/h)
      Vflowin = F/rhobulk;
%     Velocity (m/s)
      velin = Vflowin/(Axin*3600);  

%     Outlet cross sectional area of tube (m^2)
      Axout = pi*(diamout/2)^2;  
%     Volumetric flowrate (m^3/h)
      Vflowout = E/rhobulk;
%     Velocity (m/s)
      velout = Vflowout/(Axout*3600);  

%     Bulk viscosity (Pa-s); estimating as water
      mu = viscos(Fcomp,Fctype,FtempK,mumodel);  
      K = mu(1);
      n = mu(2);
 
%     Reynolds Number
      NRein = reynolds(diamin,diamin,velin,rhobulk,K,n,0);
      NReout = reynolds(diamout,diamout,velout,rhobulk,K,n,0);
      alpha = 1.0;
      if ((NRein < 2100) & (NReout < 2100)) alpha = 1/2; end
 
%     Exit pressure (MPa)
      Epres = Fpres + delP; 
 
%     Gravitational constant (m/s^2)
      g = gravity;
 
%     Shaft Work (kJ/kg)
      Ws = - (1/(2*alpha)*(velout^2-velin^2)*1e-3+1e3*(delP)/rhobulk+ ...
           g*delH*1e-3); 
 
%     Pump Work (kJ/kg)
      Wp = -Ws/eta;  
 
%     Energy loss due to inefficiencies of pump (kJ/h)
      Q = F*(-Ws/(1-eta));

%     Feed enthalpy (kJ/kg)  
      Fh = enthalpy(Fme);  

%     Reset exit temperature (K) as if all energy loss is converted to
%      heat which is transfered to the fluid. Use Bisection method.
%     Tlow = EtempK;
%     Thigh = EtempK + 50;
%     Qdiff = 1
%     while (Qdiff > 1e-2)
%       EtempK = (Tlow+Thigh)/2
%       Exit stream heat capacity (kJ/kg-K)
%       Ecp =  cp(Ecomp,Ectype,EtempK);
%       Enthalpy of exit stream (kJ/kg)
%       Eh = Ecp*(EtempK-Tref);
%       Heat difference
%       Qdiff = Q + F*Fh - E*Eh   
%       if (Qdiff > 0)
%         Tlow = EtempK 
%       elseif (Qdiff < 0)
%         Thigh = EtempK 
%       end
%       Qdiff = abs(Qdiff)
%     end
 
%     Pump Power with units conversion (kJ/h --> kW)
      ppower = F*Wp/3600;  
%     Units conversion (kW --> hp)
      hp = ppower/0.74570;  
 
%     Electrical power input (kW)
      epower = ppower/etae;  
 
%     Proposed motor efficiency based on pump power
%     Pump power data (kW)
      size = [.5,2,5,15,150];  
%     Matching pump efficiency data
      effic = [0.75,0.8,0.84,0.87,0.93];  
      if (ppower <= 0.5)
        etaprop = 0.75;  
      elseif (ppower >= 150)
        etaprop = 0.93;  
      else
        etaprop = interp1(size,effic,ppower);  
      end
 
%     Proposed electrical power input (kW)
      epowerp = ppower/etaprop;  

%     Update mass and energy info arrays
      Eme = [Emass;EtempK;Epres;Ex];

%     System output arrays
      Euarray = Eme;

      endtxt = [10 'END OF PUMP (',Uname,') BASIC TRANSPORT' 10]; 
      disp(endtxt); 

      end  % Basic Transport Level 

%  ************** ADVANCED TRANSPORT LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 4)
 
%     Unit operation information for advanced transport (sizing)
 
      end % Advanced Transport Level Assignments
 
%  *********************** ADVANCED TRANSPORT LEVEL *************************
 
      if (dlevel == 4) 
 
      endtxt = [10 'END OF PUMP (',Uname,') ADVANCED TRANSPORT' 10];
      disp(endtxt);
 
      end  % Advanced Transport Level

%  *************** MICROBIAL & QUALITY VARIABLE ASSIGNMENTS *****************

%     Feed Microbial Death Kinetics
      Dvalm = Fuarray(ulocate('microDval'));   % D value (min)
      DtempCm = Fuarray(ulocate('microDtemp'));% D value temperature basis (C)
      DtempKm = DtempCm+Tref; % D value temperature basis (K)
      Zvalm = Fuarray(ulocate('microZval'));   % Z value (C)
      Eam = Fuarray(ulocate('microEa'));       % activation energy (kcal/mol)
      Noin = Fuarray(ulocate('microload'));     % microbial load
      Fmdk = [Dvalm;DtempCm;Zvalm;Eam;Noin];  % microbial array 

%     Feed Quality Kinetics
      Dvalq = Fuarray(ulocate('qualDval'));    % D value (min)
      DtempCq = Fuarray(ulocate('qualDtemp')); % D value temperature basis (C)
      DtempKq = DtempCq+Tref; % D value temperature basis (K)
      Zvalq = Fuarray(ulocate('qualZval'));     % Z value (C)
      Eaq = Fuarray(ulocate('qualEa'));       % activation energy (kcal/mole)
      Qoin = Fuarray(ulocate('qualload'));     % quality load
      Fqual = [Dvalq;DtempCq;Zvalq;Eaq;Qoin];  % quality array 
 
%  ********************* MICROBIAL & QUALITY ASSESSMENT *********************

%     Pass microbial and quality assessment variables to exit streams
      Emdk = Fmdk;  
      Equal = Fqual;  

      if (mandq == 1) & (dlevel >= 3)   

%     Bulk temperature for microbial assessment (K)
      Tbulk = FtempK;   

%     Minimum Residence Time (min)
      tmin = 0;  
 
%     Final microbial load 
      Nout = mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,Tbulk,tmin)
 
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
      Emdk = Fmdk;  
      Emdk(length(Emdk)) = Nout;
      Equal = Fqual;
      Equal(length(Fqual)) = Qout;

      endtxt = [10 'END OF PUMP (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

      end  % microbial and quality assessment 

%     Update system output array
      Euarray = [Eme;Emdk;Equal];


%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(11:20);        % economic analysis array 
      end  

%  ********************* SYSTEM OUTPUT ARRAY ASSIGNMENT *********************

%     Size stream info arrays to full system output u-arrays  
      Euarray = ustd(Euarray);

      sys =[Euarray];    
 
%  ****************************** PRINTED OUTPUT ****************************

  if flag == 9    % print only after a solution is found 

%   Print results to specified file name; 
%    call filename manager function rfile.m 
    fid = rfile;

%   Name of unit operation; system input streams; system output streams  
    UNAME = Uname; 
    SINAMES = str2mat(Fname);
    SONAMES = str2mat(Ename);

%   System inputs array 
    SIARRAY = [Fmass]; 
%   System outputs array 
    SOARRAY = [Emass];

%   No intermediate streams
    SNNAMES = [];
    SNARRAY = []; 

%   Print mass balance results table
    msprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    if (dlevel >= 2) 
%     System inputs array  
      SIARRAY = [Fme]; 
%     System outputs array 
      SOARRAY = [Eme];

%     Print energy balance results table
      eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    end  

    if (dlevel >= 3)   
    fprintf(fid,'PUMP MODEL INPUTS\n');     
    fprintf(fid,'     \n');

    fprintf(fid,'  Inlet Diameter:                   %8.2f m\n',diamin);
    fprintf(fid,'  Outlet Diameter:                  %8.2f m\n',diamout);
    fprintf(fid,'  Pressure Differential:           %8.3e MPa\n',delP);
    fprintf(fid,'  Height Differential:              %8.2f m\n',delH);

    end  

    if (dlevel >= 3) 
    fprintf(fid,'     \n');
    fprintf(fid,'PUMP MODEL OUTPUTS\n');     
    fprintf(fid,'     \n');
    fprintf(fid,'  Inlet Velocity:                   %8.2f m/s\n',velin);
    fprintf(fid,'  Outlet Velocity:                  %8.2f m/s\n',velout);
    fprintf(fid,'  Inlet Reynold Number:            %8.3e \n',NRein);
    fprintf(fid,'  Outlet Reynold Number:           %8.3e \n',NReout);

    fprintf(fid,'     \n');
    fprintf(fid,'  Shaft Work:                     %8.3e kJ/kg\n',Ws);
    fprintf(fid,'  Pump Efficiency:                  %8.3f \n',eta);
    fprintf(fid,'  Pump Power:                      %8.3e kW\n',ppower);
    fprintf(fid,'                                   %8.3e hp\n',hp);
    fprintf(fid,'  Motor Efficiency:                 %8.3f \n',etae);
    fprintf(fid,'  Electrical Power Input:          %8.3e kW\n',epower);
    fprintf(fid,'     \n');
    fprintf(fid,'  Proposed Motor Efficiency Based   \n');
    fprintf(fid,'   on Pump Power:                   %8.3f \n',etaprop);
    fprintf(fid,'  Proposed Electrical Power Input  %8.3e kW\n',epowerp);
    fprintf(fid,'     \n');
    fprintf(fid,'  Density                           %8.2f kg/m^3\n',rhobulk);
    fprintf(fid,'  Consistency Coefficient          %8.3e Pa-s^n\n',K);
    fprintf(fid,'  Flow Behavior Index (n)           %8.3f \n',n);
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
      SOMQNAME = str2mat(Ename); 

%     System input microbial array
      SIMARRAY = [Fmdk];
%     System outputs microbial array
      SOMARRAY = [Emdk;tmin;Tbulk-Tref;];

%     No intermediate streams
      SNMQNAME = [];
      SNMARRAY = [];
      SNQARRAY = [];

      mdkprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIMARRAY,SOMARRAY,SNMARRAY);

%     System input quality array
      SIQARRAY = [Fqual];
%     System outputs quality array
      SOQARRAY = [Equal;tqmax;Tqual-Tref;];

      qprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIQARRAY,SOQARRAY,SNQARRAY);
    end 

    if (econ >= 1) 

      Uutil = [0 0 0 0 0 0]; 
      Usize = diami;
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
 
