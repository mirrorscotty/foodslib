      function [sys,x0] = mix2s(t,x,u,flag,Uname,Fname1,Fname2,Ename)

%  MIX2S     Point mix for 2 streams   
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for 2 stream point mix steady-state 
%            model. Determines mass and energy balances. 
%   
%            Call: unitops(t,x,u,flag,Uname,Fname1,Fname2,Ename)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Uname - 2 stream point mix name
%             Fname1 - feed stream 1 name 
%             Fname2 - feed stream 2 name 
%             Ename - exit stream name 
%  
%            Returned Arguments: 
%             sys - SIMULINK S-function variable: for flag = 0, 
%                   sys = sizes of parameters and initial conditions. 
%                   For flag = 3, sys = system outputs which are the 
%                   u-arrays for each stream exiting the unit operation 
%
% 

%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = 2;                            % Number of system input streams
SO = 1;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = 13;                      % Number of point mix inputs 

if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, (SO*usize), ((SI*usize)+UOinputs), 0, 1];
   x0 = ones(size(u));

elseif flag == 3 | flag == 9
   % if flag=3, return system outputs
  
%     Breakup of SIMULINK constructed u-array into stream arrays and 
%      operating information array  

      for I = 1:SI  
        Fuarray(:,I) = u((usize*(I-1)+1):(I*usize)); % full feed stream arrays  
      end  
      UOuarray = u(SI*usize+1:SI*usize+UOinputs); 
                                     % full operating inputs array

%  *********************** USER-DEFINED DESIGN LEVEL ************************

      dlevel = UOuarray(1);           % design level
      mandq = UOuarray(2);            % micro and quality assessment
      econ = UOuarray(3);             % economic evaluation  

%  ******************* MASS BALANCE LEVEL VARIABLE ASSIGNMENTS **************

%     Feed streams information  
      for I = 1:SI 
        F(I) = Fuarray(ulocate('mass'),I); % mass flow rate (kg/h)
        Fcomp(:,I) = Fuarray((ulocate('comp1st'):ulocate('complast')),I);
                                    % component mass fractions (w/w)
        Fctype(:,I) = Fuarray((ulocate('ctype1st'):ulocate('ctypelast')),I); 
                                    % component types  
      end 
      Fmass = [F;Fcomp;Fctype];     % mass information array 

%  ***************************  MASS BALANCE LEVEL ************************** 

%     Overall mass balance (kg/h)
      E = sum(F);

%     Composition of exit stream sub-algorithm 
      K = 1;  
      maxcomp = size(Fcomp(:,1),1);
      Ecomp = zeros(maxcomp,1);  
      Ectype = zeros(maxcomp,1);  
%     Lowest component (index assignment) type 
      mincomp = min(Fctype);  
      lowcomp = mincomp(1);  
      if (lowcomp == 0) lowcomp = 1; end  
%     Highest component (index assignment) type 
      maxcomp = max(Fctype);  
      highcomp = maxcomp(1);  
%     For each component type  
      for J = lowcomp:highcomp 
%       For each feed stream 
        for I = 1:SI 
%         Find mass fraction of specified component type J 
          comp = findcomp(Fcomp(:,I),Fctype(:,I),J);
%         Mass flow of component type J in feed stream I 
          Mcomp(I) = F(I)*comp;
        end  
%       Sum of all mass flows of component type J 
        Etotcomp = sum(Mcomp);  
        if (Etotcomp > 0) 
          Ecomp(K) = Etotcomp/E;
          Ectype(K) = J;  
          K = K + 1;
        end 
        if (K > 10) & (J < highcomp) 
          errortxt = [10 'POINT MIX (',Uname,') ERROR: Number of component types being mixed' 10 ...
          'exceeds 10.' 10];  
          disp(errortxt); 
          return;
        end   
      end % for J components 

%     Exit stream mass arrays; used for print results table  
      Emass = [E;Ecomp;Ectype];

%     System output arrays 
      Euarray = Emass;  

      endtxt = [10 'END OF POINT MIX (',Uname,') MASS BALANCE' 10];  
      disp(endtxt); 

%  **************** ENERGY BALANCE LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 2)  

%     Reference temperature (K)
      Tref = tref;

%     Feed stream information  
      FtempC = Fuarray(ulocate('temp'),:);% temperature (C) 
      FtempK = FtempC + Tref;            % temperature (K) 
      Fpres = Fuarray(ulocate('pressure'),:); % pressure (MPa) 
      Fx = Fuarray(ulocate('quality'),:);    % quality  
      Fme = [Fmass;FtempK;Fpres;Fx]; % mass and energy information array  

%  *********************** ENERGY BALANCE LEVEL *****************************

%     Enthalpy of feed streams (kJ/kg)
      for I = 1:SI 
        Fh(I) = enthalpy(Fme(:,I));
      end  

%     Total enthalpy of exit stream (kJ/kg) 
      Ehtot = sum(F.*Fh); 

%     Exit temperature (K) determined by bisection method
      Tlow = min(FtempK);
      Thigh = max(FtempK);
      if (Tlow == Thigh) 
        EtempK = Tlow; 
      else  
        Tdiff = 1;
        while (abs(Tdiff) > 1e-8)
          EtempK = (Tlow+Thigh)/2;
          Ecp = cp(Ecomp,Ectype,EtempK);
          EtempKc = Ehtot/(E*Ecp) + Tref;
          Tdiff = EtempK - EtempKc;
          if (Tdiff > 0)
            Thigh = EtempK;
          else
            Tlow = EtempK;
          end
        end
      end  

%     Exit stream pressure (MPa) and quality  
      Epres = mean(Fpres);
      Ex = 0; 

%     Update mass and energy info arrays 
      Eme = [Emass;EtempK;Epres;Ex]; 

%     System output arrays
      Euarray = Eme;

      endtxt = [10 'END OF POINT MIX (',Uname,') ENERGY BALANCE' 10];
      disp(endtxt);

      end % Energy Balance Level  

%  *************** MICROBIAL & QUALITY VARIABLE ASSIGNMENTS *****************

%     Feed Microbial Death Kinetics
      Dvalm = Fuarray(ulocate('microDval'),:);   % D value (min)
      DtempCm =  Fuarray(ulocate('microDtemp'),:);% D value temperature basis (C)
      DtempKm = DtempCm+Tref; % D value temperature basis (K)
      Zvalm =  Fuarray(ulocate('microZval'),:);   % Z value (C)
      Eam =  Fuarray(ulocate('microEa'),:);       % activation energy (kcal/mol)
      Noin =  Fuarray(ulocate('microload'),:);     % microbial load
      Fmdk = [Dvalm;DtempCm;Zvalm;Eam;Noin];  % microbial array 

%     Feed Quality Kinetics
      Dvalq = Fuarray(ulocate('qualDval'),:);    % D value (min)
      DtempCq = Fuarray(ulocate('qualDtemp'),:); % D value temperature basis (C)
      DtempKq = DtempCq+Tref; % D value temperature basis (K)
      Zvalq = Fuarray(ulocate('qualZval'),:);     % Z value (C)
      Eaq = Fuarray(ulocate('qualEa'),:);       % activation energy (kcal/mole)
      Qoin =  Fuarray(ulocate('qualload'),:);     % quality load
      Fqual = [Dvalq;DtempCq;Zvalq;Eaq;Qoin];  % quality array 

%  ********************* MICROBIAL & QUALITY ASSESSMENT *********************

%     The microbe with the largest D-value will continue to be tracked
%     Find stream with microbial kinetics with largest D value (min)
      Jguess = find(Dvalm == max(Dvalm));
      J = Jguess(1);
%     Pass kinetics for that microbe
      Dvalmax = Dvalm(J);
      Dtempmax = DtempKm(J);
      Zvalmax = Zvalm(J);
      Eamax = Eam(J);
%     Sum microbial loads for all streams tracking the desired microbe
      sumNin = 0;
      for I = 1:SI  
        if (Dvalm(I) == Dvalmax) & (Zvalm(I) == Zvalmax) & (Eam(I) == Eamax)
          sumNin = sumNin + F(I)*Noin(I);
        end
      end  
%     Exit stream microbial load
      Nout = sumNin/E;
 
%     The quality parameter with the smallest D-value will continue to be
%      tracked
%     Find stream with quality having smallest D value (min)
      Jguess = find(Dvalq == min(Dvalq));  
      J = Jguess(1);  
%     Pass kinetics for that microbe
      Dvalqmin = Dvalq(J);
      Dtqmin = DtempKq(J);
      Zvalqmin = Zvalq(J);
      Eaqmin = Eaq(J);
%     Sum quality for all streams tracking the desired quality parameter
      sumQin = 0;
      for I = 1:SI  
        if (Dvalq(I) == Dvalqmin) & (Zvalq(I) == Zvalqmin) & (Eaq(I) == Eaqmin)
          sumQin = sumQin + F(I)*Qoin(I);
        end
      end
%     Exit stream quality
      Qout = sumQin/E;

%     Pass microbial and quality assessment variables to exit streams
      Emdk = [Dvalmax;Dtempmax-Tref;Zvalmax;Eamax;Nout];  
      Equal = [Dvalqmin;Dtqmin-Tref;Zvalqmin;Eaqmin;Qout];  

      if (mandq == 1) & (dlevel >= 2)   

%     Mixed microbial and quality loads 
      Nmix = Nout;
      Qmix = Qout;

%     Bulk temperature for microbial assessment (K)
      Tbulk = EtempK;   

%     Minimum Residence Time (min)
      tmin = 0;  
 
%     Final microbial load 
      Nout = mdeathk(Dvalmax,Dtempmax,Zvalmax,Eamax,Nmix,Tbulk,tmin);  
 
%     Log reduction in microbial load
      DlogNred = log10(Nmix/Nout);

%     Temperature for quality assessment (K)
      Tqual = EtempK;

%     Maximum processing time = Minimum processing time (min) 
      tqmax = 0;

%     Final quality load  
      Qout =  qualdegr(Dvalqmin,Dtqmin,Zvalqmin,Eaqmin,Qmix,Tqual,tqmax); 

%     Log reduction in quality 
      DlogQred = log10(Qmix/Qout);

%     Exit stream 1 microbial and quality arrays
      Emdk(size(Emdk,1)) = Nout;
      Equal(size(Equal,1)) = Qout;

      endtxt = [10 'END OF UNIT (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

      end  % microbial and quality assessment 

%     Update system output array
      E1uarray = [Eme;Emdk;Equal];

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(19:28); % economic analysis array 
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
    for I = 1:SI  
      FnameF = ['Fname' int2str(I)];
      if exist(FnameF) 
        Fnsize = size(eval([FnameF]),2);
        SINAMES(I,1:Fnsize) = eval([FnameF]);  
      end  
    end  
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
      SIMQNAME = SINAMES;
      SOMQNAME = str2mat(Ename); 

%     System input microbial array
      SIMARRAY = [Fmdk];
%     System outputs microbial array
      SOMARRAY = [Emdk;tmin;Tbulk-Tref];

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
      Usize = 1;
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
 
