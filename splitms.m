      function [sys,x0] = splitms(t,x,u,flag,Uname,Fname,Ename1,Ename2,Ename3,Ename4,Ename5)

%  SPLITMS   Stream mass split model   
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for stream mass split steady-state 
%            model. Determines mass and energy balances. 
%   
%            Call: splitms(t,x,u,flag,Uname,Fname,E1name,E2name,E2name,
%                          E3name,E4name,E5name)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Uname - mass splitter name
%             Fname - feed stream name 
%             Ename1 - exit stream 1 name 
%             Ename2 - exit stream 2 name 
%             Ename3 - exit stream 3 name 
%             Ename4 - exit stream 4 name 
%             Ename5 - exit stream 5 name 
%  
%            Returned Arguments: 
%             sys - SIMULINK S-function variable: for flag = 0, 
%                   sys = sizes of parameters and initial conditions. 
%                   For flag = 3, sys = system outputs which are the 
%                   u-arrays for each stream exiting the unit operation 
%
% 

%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = 1;                            % Number of system input streams
SO = 5;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = 18;                      % Number of mass split inputs 

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

%     Mass split information for mass balances  
%     Maximum number of outlet streams
      maxouts = 5;
      for I = 1:maxouts
        Fsplit(I) = UOuarray(3+I);   % mass flowrate split fractions
      end

%  ***************************  MASS BALANCE LEVEL ************************** 

%     Number of non-zero mass fraction splits
      scount = nnz(Fsplit);  
 
%     Check that fractions do not exceed 1
      chkfrac = 1 - sum(Fsplit);
%     Send remaining fraction to next available stream
      if (chkfrac > 0) & (scount < maxouts)
        Jguess = find(Fsplit == 0); 
        J = Jguess(1);  
        Fsplit(J) = chkfrac;  
        scount = scount + 1;  
      elseif (chkfrac > 0) & (scount == maxouts)
        errortxt = [10 'MASS SPLIT (',Uname,') ERROR: Split fractions do not sum to 1, and there' 10 ...
        'is no available stream to which the remaining fraction can to be passed.' 10]; 
        disp(errortxt);  
        return;
      elseif (chkfrac < 0)
        errortxt = [10 'MASS SPLIT (',Uname,') ERROR: Split fractions sum exceeds 1.' 10];
        disp(errortxt);  
        return;
      end
 
%     Overall mass balance (kg/h)
      E = F .* Fsplit;  

%     Exit stream composition and types (w/w)
      Ecomp = zeros(size(Fcomp,1),maxouts);
      Ectype = zeros(size(Fcomp,1),maxouts);
      for I = 1:maxouts
        if (E(I) > 0)
          Ecomp(:,I) = Fcomp;
          Ectype(:,I) = Fctype;
        end       
      end 

%     Exit stream mass arrays; used for print results table  
      Emass = [E;Ecomp;Ectype];

%     System output arrays 
      Euarray = Emass;  

      endtxt = [10 'END OF MASS SPLIT (',Uname,') MASS BALANCE' 10];  
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

%     Enthalpy of feed stream (kJ/kg)
      Fh = enthalpy(Fme);  

%     Exit temperature (K), pressure(MPa), quality
      for I = 1:maxouts   
        if (Fsplit(I) > 0) 
          EtempK(I) = FtempK;
          Epres(I) = Fpres;
          Ex(I) = Fx; 
        else 
          EtempK(I) = 0;
          Epres(I) = 0;
          Ex(I) = 0;  
        end  
      end  

%     Exit pressure (MPa) 

%     Update mass and energy info arrays 
      Eme = [Emass;EtempK;Epres;Ex]; 

%     System output arrays
      Euarray = Eme;

      endtxt = [10 'END OF MASS SPLIT (',Uname,') ENERGY BALANCE' 10];
      disp(endtxt);

      end % Energy Balance Level  

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
      Emdk = zeros(size(Fmdk,1),maxouts);
      Equal = zeros(size(Fqual,1),maxouts);
      for I = 1:maxouts 
        if (Fsplit(I) > 0) 
          Emdk(:,I) = Fmdk;  
          Equal(:,I) = Fqual;  
        end 
      end  

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
      for I = 1:maxouts 
        if (Fsplit(I) > 0) 
          Emdk(:,I) = Fmdk;  
          Emdk(size(Fmdk,1),I) = Nout;
          Equal(:,I) = Fqual;
          Equal(size(Fqual,1),I) = Qout;
        end  
      end  

      endtxt = [10 'END OF UNIT (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

      end  % microbial and quality assessment 

%     Update system output array
      Euarray = [Eme;Emdk;Equal];

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(9:18); % economic analysis array 
      end  

%  ********************* SYSTEM OUTPUT ARRAY ASSIGNMENT *********************

%     Size stream info arrays to full system output u-arrays  
      Euarray = ustd(Euarray);

      sys =[Euarray(:,1);Euarray(:,2);Euarray(:,3);Euarray(:,4);Euarray(:,5);];    
 
%  ****************************** PRINTED OUTPUT ****************************

  if flag == 9    % print only after a solution is found 

%   Print results to specified file name; 
%    call filename manager function rfile.m 
    fid = rfile;

%   Name of unit operation; system input streams; system output streams  
    UNAME = Uname; 
    SINAMES = str2mat(Fname);
    SONAMES = str2mat(Ename1,Ename2,Ename3,Ename4,Ename5);

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
      SOARRAY = []; 
      for I = 1:maxouts 
        if (Fsplit(I) > 0) 
          SOARRAY = [SOARRAY Eme(:,I)];
        end 
      end

%     Print energy balance results table
      eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    end  

    fprintf(fid,'UNIT MODEL INPUTS\n');     
    fprintf(fid,'     \n');
    for I = 1:maxouts 
      if(Fsplit(I) > 0)   
        fprintf(fid,'  Mass Flow Split to Exit %1.0f      %5.4f w/w \n',I,Fsplit(I)); 
      end
    end  
    fprintf(fid,'     \n');

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
      SOMQNAME = str2mat(Ename1,Ename2,Ename3,Ename4,Ename5); 

%     System input microbial array
      SIMARRAY = [Fmdk];
%     System outputs microbial array
      SOMARRAY = []; 
      for I = 1:maxouts
        if (Fsplit(I) > 0) 
          Emdkpr = [Emdk(:,I);tmin;Tbulk-Tref;];
          SOMARRAY = [SOMARRAY Emdkpr];
        end 
      end  

%     No intermediate streams
      SNMQNAME = [];
      SNMARRAY = [];
      SNQARRAY = [];

      mdkprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIMARRAY,SOMARRAY,SNMARRAY);

%     System input quality array
      SIQARRAY = [Fqual];
%     System outputs quality array
      SOQARRAY = [];  
      for I = 1:maxouts
        if (Fsplit(I) > 0) 
          Equalpr = [Equal(:,I);tqmax;Tqual-Tref;]; 
          SOQARRAY = [SOQARRAY Equalpr];
        end  
      end   

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
 
