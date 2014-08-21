      function [sys,x0] = pipings(t,x,u,flag,Uname,Fname,Ename)

%  PIPEINGS  Straight pipe model   
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for straight pipe steady-state 
%            model. Determine mass and energy balances; performs design at 
%            basic and advanced transport level. 
%   
%            Call: unitops(t,x,u,flag,Uname,Fname,Ename)
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
UOinputs = 28;                     % Number of straight pipe inputs 

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

      endtxt = [10 'END OF STRAIGHT PIPE (',Uname,') MASS BALANCE' 10];  
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

      endtxt = [10 'END OF STRAIGHT PIPE (',Uname,') ENERGY BALANCE' 10];
      disp(endtxt);

      end % Energy Balance Level  

%  *************** BASIC TRANSPORT LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 3)

%     Unit operation information for basic transport (sizing) 

      end % Basic Transport Level Assignments  

%  ************************** BASIC TRANSPORT LEVEL *************************

      if (dlevel == 3) 

      endtxt = [10 'END OF UNIT (',Uname,') BASIC TRANSPORT' 10];
      disp(endtxt);

      end  % Basic Transport Level 

%  ************* ADVANCED TRANSPORT LEVEL VARIABLE ASSIGNMENTS **************

      if (dlevel >= 4) 

%     Pipe fitting inforamtion for advanced transport  
      diami = UOuarray(4);      % pipe inner diameter (m)
      diamo = UOuarray(5);      % pipe outer diameter (m)
      length = UOuarray(6);     % pipe length (m)
      vhorient = UOuarray(7);   % vertical (0) or horizontal (1) orientation
      rough = UOuarray(8);      % pipe roughness (m)
      k(1) = UOuarray(9);       % pipe thermal conductivity (kJ/h-m-K)
      insul(1) = UOuarray(10); % 1st insulation layer thinkness (m)
      k(2) = UOuarray(11);     % 1st insulation thermal conductivity (kJ/h-m-K)
      insul(2) = UOuarray(12); % 2nd insulation layer thinkness (m)
      k(3) = UOuarray(13);     % 2nd insulation thermal conductivity (kJ/h-m-K)
      insul(3) = UOuarray(14); % 3rd insulation layer thinkness (m)
      k(4) = UOuarray(15);     % 3rd insulation thermal conductivity (kJ/h-m-K)
      mumodel = UOuarray(16);  % viscosity model 
      TambC = UOuarray(17);     % ambient air temperature (C) 
      Tamb = TambC + Tref;     % (K)  
      Pamb = UOuarray(18);     % ambient air pressure (MPa)  

%  ************************ ADVANCED TRANSPORT LEVEL ************************

%     Resitances to heat transfer due to metal and insulating layers
%     Assume maximum layers is 3
      lmax = 3;
%     Radii (m)
      r(1) = diami/2;
      r(2) = diamo/2;
      layers = 1;  % Metal layer
      for I = 1:lmax
        if (insul(I) > 0)
          r(I+2) = r(I+1) + insul(I);  
          layers = layers + 1;   % Additional insulating layers
        end
      end
 
%     Surface area of each layer (m^2)
      A = 2*pi*length*r;  
 
%     Log mean areas for layers (m^2)
      for I = 1:layers
        Alm(I) = (A(I+1)-A(I))/log(A(I+1)/A(I));  
      end
 
%     Heat transfer resistances (h-K/kJ)
%     Metal and insulating layers resistance
      for I = 1:layers
        resist(I) = (r(I+1)-r(I))/(k(I)*Alm(I));  
      end
%     Total resistance  of layers (h-K/kJ)
      Rlayers = sum(resist);  
 
%     Convergence criteria: Iteration difference between last and current
%      bulk temperature (K)
      delTbulk = 1;  
      IT = 0; 
      Q = 0;  
      Ro = 0;  
      while (delTbulk > 1e-8)
        Tbulk = (FtempK+EtempK)/2; 
 
%       Bulk density (kg/m^3)
        FEme = Fme;  
        xtemp = ulocate('temp');  
        FEme(xtemp) = Tbulk;  
        rhobulk = density(FEme);
 
%       Inner cross sectional area of tube (m^2)
        Axi = pi*(diami/2)^2;   
%       Volumetric flowrate (m^3/h)
        Volflow = F/rhobulk;
%       Velocity (m/s)
        velocity = Volflow/(Axi*3600);  

%       Viscosity 
        powerlaw = viscos(Fcomp,Fctype,FtempK,mumodel);  
        Kbulk = powerlaw(1);  % Consistency Coeff (Pa-s^n)   
        nbulk = powerlaw(2);  % Flow Behavior Index  

%       Bulk viscosity (Pa-s); taken at max velocity in tube
        mubulk = Kbulk;
 
%       Reynolds Number
        NRe = reynolds(diami,diami,velocity,rhobulk,Kbulk,nbulk,0);
 
%       Heat capacity of exit stream (kJ/kg-K)
        cpbulk = cp(Fcomp,Fctype,Tbulk);
 
%       Thermal conductivity (kJ/h-m-K)
        kbulk = thermc(Fcomp,Fctype,Tbulk);
 
%       Prandtl Number
        NPr = prandtl(cpbulk,mubulk,kbulk);
 
%       Viscosity at wall  
        if (IT == 0)
          muwall = mubulk;  
        else
%         Inside wall temperature (K)
          Tis = Tamb - Q*(Rlayers+Ro);  
%         Power law values for viscosity at wall temperature
          powerlaw = viscos(Fcomp,Fctype,Tis,mumodel);
          Kwall = powerlaw(1);       % Consistency Coefficient (Pa-s^n)
          nwall = powerlaw(2);       % Flow behavior index

%         Maximum shear rate (1/s)
          shearmax = shrwall(nwall,Volflow,diami);

%         Viscosity at tube wall (Pa-s)
          muwall = Kwall*shearmax^(nwall-1);
        end
 
%       Nusselt Number
        if (NRe < 2100)  % Laminar Flow
          NNu = 1.86*(NRe*NPr*D/L)^(1/3)*(mubulk/muwall)^0.14;
        elseif (NRe > 1e4)  % Turbulent Flow
          NNu = 0.023*NRe^0.8*NPr^(1/3)*(mubulk/muwall)^0.14;
        else   % Transition flow; Hausen's Equation (Perry's)
          NNu = 0.116*(NRe^(2/3)-125)*NPr^(1/3)*(1+(diami/length)^(2/3))*(mubulk/muwall)^0.14;
        end
 
%       Inside heat transfer coefficient (kJ/m^2-h-K)
        hi = NNu*kbulk/diami;
 
%       Outside film temperature (K)
        if (IT == 0)
          Tfilm = Tamb;  
        else
          Tos = Tbulk + Q*(Rlayers+Ri);  
          Tfilm = (Tamb + Tos)/2;  
        end
 
%       Air density (kg/m^3) at film temperature
        denair = rhoair(Tfilm,Pamb);  
 
%       Volumetric coefficient of expansion of air (1/K) at film temp.
        Beta = 1/Tfilm;  
 
%       Air viscosity (Pa-s) at film temperature
        muairf = muair(Tfilm);  
 
%       Positive temperature difference between the wall and the air (K)
%       delTair = abs(Tfilm - Tamb)/2;  
        delTair = abs(Tfilm - Tamb);  
 
%       Grashof Number
        NGr = (r(layers+1))^3*denair^2*9.80665*Beta*delTair/muairf^2;  
 
%       Air heat capacity (kJ/kg-K) at film temperature
        cpairf = cpair(Tfilm);  
 
%       Air thermal conductivity (kJ/h-m-K) at film temperature
        kairf = kair(Tfilm);  
 
%       Prandtl Number for Natural Convection of air
        NPrair = cpairf*muairf*3600/kairf;  
 
%       Product of Grashof and Prandtl Numbers
        NGrNPr = NGr*NPrair;  
 
%       Selection of constants for natural convection
%       Good for horizontal pipes with D < 0.2 m
        if (vhorient == 0)
          if (NGrNPr < 1e-5)
            a = 0.49;
            m = 0;
          elseif (NGrNPr >= 1e-5) & (NGrNPr < 1e-3)
            a = 0.71;
            m = 1/25;
          elseif (NGrNPr >= 1e-3) & (NGrNPr < 1)
            a = 1.09;
            m = 1/10;
          elseif (NGrNPr >= 1) & (NGrNPr < 1e4)
            a = 1.09;
            m = 1/5;
          elseif (NGrNPr >= 1e4) & (NGrNPr < 1e9)
            a = 0.53;
            m = 1/4;
          elseif (NGrNPr >= 1e9)
            a = 0.13;
            m = 1/3;
          end
%       Good for vertical pipes with L < 1 m
        elseif (vhorient == 1)
          if (NGrNPr(S) < 1e4)
            a = 1.36;
            m = 1/5;
          elseif (NGrNPr(S) >= 1e4) & (NGrNPr(S) < 1e9)
            a = 0.59;
            m = 1/4;
          elseif (NGrNPr(S) >= 1e9)
            a = 0.13;
            m = 1/3;
          end
        end
 
%       Outside Nusselts Number
        NNuair = a*(NGrNPr)^m;  
 
%       Outside heat transfer coefficient (kJ/m^2-h-K)
        ho = NNuair*kairf/(r(layers+1));  
 
%       Inside resistance
        Ri = 1/(hi*A(1));  
%       Outside resistance
        Ro = 1/(ho*A(layers+1));  
 
%       Heat Transfer Rate (kJ/h)
        Q =  (Tamb - Tbulk)/(Ri+Rlayers+Ro);  

%       Enthalpy of feed stream (kJ/kg) 
        Fh = enthalpy(Fme);

%       Enthalpy of exit stream (kJ/kg) 
        Eh = (F*Fh + Q)/F; 
 
%       Find exit temperature (K) by Bisection method 
        if (Q > 0) 
          EtempKlow = FtempK;  
          EtempKhigh = Tamb;  
        else
          EtempKlow = Tamb;  
          EtempKhigh = FtempK;   
        end  
        Ehdiff = 1;  
        while (Ehdiff > 1e-4)
          EtempK = (EtempKhigh+EtempKlow)/2;
          Eme = [Emass;EtempK;Epres;Ex];
          Ehnew = enthalpy(Eme);  
          Ehdiff = Eh - Ehnew;  
          if (Ehdiff > 0)  
            EtempKlow = EtempK;
          else
            EtempKhigh = EtempK;  
          end 
        end  
        
%       EtempK = (F*Fh + Q)/(F*cpbulk) + Tref
 
%       Converegence Criteria
        delTbulk = abs(Tbulk - (FtempK+EtempK)/2);  
        IT = IT + 1; 
      end
 
%     Fanning friction factor
      if (NRe < 2100)  % Laminar Flow
        f = 16/NRe;
      else     % Turbulent flow: Colebrook equation (Perry's)
        friclow = 0.001;
        frichigh = 0.25;
        fdiff = 1;
        while (fdiff > 1e-5)
          f = (friclow+frichigh)/2;
          fdiff = -4*log10(rough/(3.7*diami)+1.256/(NRe*sqrt(f))) - 1/sqrt(f);
          if (fdiff > 0)
            frichigh = f;
          else
            friclow = f;
          end
          fdiff = abs(fdiff);
        end
      end
 
%     Pressure drop (MPa)
      delP = (4*f*rhobulk*length/diami*velocity^2/2)/1e6;
%     Exit pressure (MPa)
      Epres = Fpres - delP;

%     Update mass and energy info arrays
      Eme = [Emass;EtempK;Epres;Ex];

%     System output arrays
      Euarray = Eme;
 
      endtxt = [10 'END OF STRAIGHT PIPE (',Uname,') ADVANCED TRANSPORT' 10]; 
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

%     Pipe volume (m^3)
      pipevol = Axi*length;  
%     Average time in piping (h)
      timeh = pipevol/Volflow;  
%     Time units conversion (h --> min)
      time = timeh*60;  

%     Bulk temperature for microbial assessment (K)
      Tbulk = Tbulk;   

%     Minimum Residence Time (min)
      tmin = time;  
 
%     Final microbial load 
      Nout = mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,Tbulk,tmin);  
 
%     Log reduction in microbial load
      DlogNred = log10(Noin/Nout);

%     Temperature for quality assessment (K)
      Tqual = FtempK;

%     Maximum processing time = Minimum processing time (min) 
      tqmax = time;

%     Final quality load  
      Qout =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,tqmax); 

%     Log reduction in quality 
      DlogQred = log10(Qoin/Qout);

%     Exit stream microbial and quality arrays
      Emdk = Fmdk;  
      Emdk(size(Fmdk,1)) = Nout;
      Equal = Fqual;
      Equal(size(Fqual,1)) = Qout;

      endtxt = [10 'END OF STRAIGHT PIPE (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

%     Update system output array
      Euarray = [Eme;Emdk;Equal];

      end  % microbial and quality assessment 

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(19:28);        % economic analysis array 
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

    if (dlevel >= 4)   
    fprintf(fid,'UNIT MODEL INPUTS\n');     
    fprintf(fid,'     \n');

    fprintf(fid,'  Inner Diameter                    %8.4f m\n',diami);
    fprintf(fid,'  Outer Diameter                    %8.4f m\n',diamo);
    fprintf(fid,'  Length                            %8.4f m\n',length);
    if (vhorient == 0)
    fprintf(fid,'  Pipe Orientation                  %s \n','vertical');
    else 
    fprintf(fid,'  Pipe Orientation                  %s \n','horizontal');
    end
    fprintf(fid,'  Pipe Roughness                   %8.3e m\n',rough);
    fprintf(fid,'  Viscosity Model                   %8.0f \n',mumodel);
    fprintf(fid,'  Ambient Temperature               %8.2e C\n',TambC);
    fprintf(fid,'  Ambient Pressure                 %8.3e MPa\n',Pamb);

    end  

    if (dlevel >=4) 
    fprintf(fid,'     \n');
    fprintf(fid,'UNIT MODEL OUTPUTS\n');     
    fprintf(fid,'     \n');
    fprintf(fid,'  Velocity                          %8.4f m/s\n',velocity);
    fprintf(fid,'  Reynold Number                   %8.3e \n',NRe);
    fprintf(fid,'  Pipe Roughness                   %8.3e m\n',rough);
    fprintf(fid,'  Fanning Friction Factor           %8.4f \n',f);
    fprintf(fid,'  Prandtl Number                    %8.2f \n',NPr);
    fprintf(fid,'  Head Loss                        %8.3e MPa\n',delP);
    fprintf(fid,'  L/Di ratio                        %8.2f \n',length/diami);
    fprintf(fid,'  Inside Heat Transfer Coeff.      %8.3e kJ/h-m^2-K \n',hi);
 
    fprintf(fid,'     \n');
    fprintf(fid,'  NGr*NPr                          %8.3e  \n',NGrNPr);
    fprintf(fid,'  Outside Heat Transfer Coeff.     %8.3e kJ/h-m^2-K \n',ho);
 
    fprintf(fid,'     \n');
    fprintf(fid,'  Heat Transfer Rate              %8.3e kJ/h \n',Q);
    fprintf(fid,'     \n');
    fprintf(fid,'  Mean Temperature Profile \n');
    fprintf(fid,'   Bulk Temperature                 %8.3f C\n',Tbulk-Tref);
    fprintf(fid,'   Inside Surface Temperature       %8.3f C\n',Tis-Tref);
    fprintf(fid,'   Outside Surface Temperature      %8.3f C\n',Tos-Tref);
    fprintf(fid,'   Outside Film Temperature         %8.3f C\n',Tfilm-Tref);
    fprintf(fid,'   Ambient Temperature              %8.3f C\n',TambC);
 
    fprintf(fid,'     \n');
    fprintf(fid,'  Density                           %8.2f kg/m^3\n',rhobulk);
    fprintf(fid,'  Bulk Consistency Coefficient     %8.3e Pa-s^n\n',Kbulk);
    fprintf(fid,'  Bulk Flow Behavior Index (n)      %8.3f \n',nbulk);
    fprintf(fid,'  Wall Consistency Coefficient     %8.3e Pa-s^n\n',Kwall);
    fprintf(fid,'  Wall Flow Behavior Index (n)      %8.3f \n',nwall);
    fprintf(fid,'     \n');

    fprintf(fid,'Piping Resistance to Heat Transfer \n');
    fprintf(fid,'     \n');
    fprintf(fid,'Layer         Thickness   Inner    Outer   Log Mean   Thermal   Resistance\n');
    fprintf(fid,'                          Radius   Radius    Area      Cond   \n');
    fprintf(fid,'                 (m)       (m)      (m)     (m^2)    (kJ/h-m-K)   (h-K/kJ) \n');
    fprintf(fid,'----------------------------------------------------------------------------\n');
    fprintf(fid,'pipe           %5.4f     %5.4f   %5.4f    %5.4f   %8.4f   %8.4e   \n',(diamo-diami),r(1),r(2),Alm(1),k(1),resist(1));
    for I = 1:layers-1
    fprintf(fid,'insulation %1.0f   %5.4f     %5.4f   %5.4f    %5.4f   %8.4f   %8.4e   \n',I,insul(I),r(I+1),r(I+2),Alm(I+1),k(I+1),resist(I+1));
    end 
    fprintf(fid,'----------------------------------------------------------------------------\n');
 
    fprintf(fid,'     \n');
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
 
