      function [sys,x0] = template(t,x,u,flag,Uname,F1name,F2name,Sname,E1name,E2name)

%  UNITOPS   Unit operation template file  
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function TEMPLATE for unit operation steady-state 
%            model. May determine mass and energy balances, perform design at 
%            basic and advanced transport level. 
%   
%            Call: unitops(t,x,u,flag,Uname,F1name,F2name,Sname,E1name,E2name)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Uname - unit operation name
%             F1name - feed stream 1 name 
%             F2name - feed stream 2 name 
%              ...
%             Sname - water/steam stream name 
%             E1name - exit stream 1 name 
%             E2name - exit stream 2 name 
%  
%            Returned Arguments: 
%             sys - SIMULINK S-function variable: for flag = 0, 
%                   sys = sizes of parameters and initial conditions. 
%                   For flag = 3, sys = system outputs which are the 
%                   u-arrays for each stream exiting the unit operation 
%
%            References:  
% 

%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = ?;                            % Number of system input streams
SO = ?;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = ?;                      % Number of unit operating inputs 

if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, (SO*usize), ((SI*usize)+UOinputs), 0, 1];
   x0 = ones(size(u));

elseif flag == 3 | flag == 9
   % if flag=3, return system outputs
  
%     Breakup of SIMULINK constructed u-array into stream arrays and 
%      operating information array  
 
      F1uarray = u(1:usize);          % full feed stream 1 array 
      F2uarray = u(usize+1:2*usize);  % full feed stream 2 array 
      UOuarray = u(SI*usize+1:SI*usize+UOinputs); 
                                     % full operating inputs array

%  *********************** USER-DEFINED DESIGN LEVEL ************************

      dlevel = UOuarray(1);           % design level
      mandq = UOuarray(2);            % micro and quality assessment
      econ = UOuarray(3);             % economic evaluation  

%  ******************* MASS BALANCE LEVEL VARIABLE ASSIGNMENTS **************

%     Feed stream 1 information  
      F1 = F1uarray(ulocate('mass')); % mass flow rate (kg/h)
      F1comp = F1uarray(ulocate('comp1st'):ulocate('complast'));
                                    % component mass fractions (w/w)
      F1ctype = F1uarray(ulocate('ctype1st'):ulocate('ctypelast')); 
                                    % component types  
      F1mass = [F1;F1comp;F1ctype];     % mass information array 

%     REPEAT feed stream information read-ins for each feed stream 

%     Unit operation information for mass balances  


%  ***************************  MASS BALANCE LEVEL ************************** 


%     Exit stream mass arrays; used for print results table  
      E1mass = [E1;E1comp;E1ctype];
      E2mass = [E2;E2comp;E2ctype];

%     System output arrays 
      E1uarray = E1mass;  
      E2uarray = E2mass;   

      endtxt = [10 'END OF UNIT (',Uname,') MASS BALANCE' 10];  
      disp(endtxt); 

%  **************** ENERGY BALANCE LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 2)  

%     Reference temperature (K)
      Tref = tref;

%     Feed stream information  
      F1tempC = F1uarray(ulocate('temp'));% temperature (C) 
      F1tempK = F1tempC + Tref;            % temperature (K) 
      F1me = [F1mass;F1tempK;F1uarray(ulocate('pressure'):ulocate('quality'))]; 
                                         % mass and energy information array  
%     Steam supply information
      S = Suarray(ulocate('mass'));      % mass flow rate (kg/h)
      Scomp = Suarray(ulocate('comp1st'):ulocate('complast'));   
                                         % component mass fractions (w/w)
      Sctype = Suarray(ulocate('ctype1st'):ulocate('ctypelast')); 
                                         % component types  
      Smass = [S;Scomp;Sctype];          % mass information array 
      StempC = Suarray(ulocate('temp')); % temperature (C) 
      StempK = StempC + Tref;            % temperature (K)  
      Spres = Suarray(ulocate('pressure'));    % pressure (MPa)   
      Sx = Suarray(ulocate('quality'));   % quality
      Sme = [Smass;StempK;Spres;Sx];
                                         % mass and energy information array 

%     Unit operation information for energy balances  

%  *********************** ENERGY BALANCE LEVEL *****************************

%     Update mass information arrays
      Smass = [S;Scomp;Sctype];
      E1mass = [E1;E1comp;E1ctype];   
      E2mass = [E2;E2comp;E2ctype];

%     Update mass and energy info arrays 
      Sme = [Smass;StempK;Spres;Sx]; % if energy balance used for mass flow rate 
      E1me = [E1mass;E1tempK;E1pres;E1x]; 
      E2me = [E2mass;E2tempK;E2pres;E2x]; 

%     System output arrays
      E1uarray = E1me;
      E2uarray = E2me;

      endtxt = [10 'END OF UNIT (',Uname,') ENERGY BALANCE' 10];
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

%     Unit operation inforamtion for advanced transport  


%  ************************ ADVANCED TRANSPORT LEVEL ************************



      endtxt = [10 'END OF UNIT (',Uname,') ADVANCED TRANSPORT' 10]; 
      disp(endtxt); 

      end  % Advanced Transport Level 

%  *************** MICROBIAL & QUALITY VARIABLE ASSIGNMENTS *****************

%     Feed Microbial Death Kinetics
      Dvalm = F1uarray(ulocate('microDval'));   % D value (min)
      DtempCm =  F1uarray(ulocate('microDtemp'));% D value temperature basis (C)
      DtempKm = DtempCm+Tref; % D value temperature basis (K)
      Zvalm =  F1uarray(ulocate('microZval'));   % Z value (C)
      Eam =  F1uarray(ulocate('microEa'));       % activation energy (kcal/mol)
      Noin =  F1uarray(ulocate('microload'));     % microbial load
      F1mdk = [Dvalm;DtempCm;Zvalm;Eam;Noin];  % microbial array 

%     Feed Quality Kinetics
      Dvalq = F1uarray(ulocate('qualDval'));    % D value (min)
      DtempCq = F1uarray(ulocate('qualDtemp')); % D value temperature basis (C)
      DtempKq = DtempCq+Tref; % D value temperature basis (K)
      Zvalq = F1uarray(ulocate('qualZval'));     % Z value (C)
      Eaq = F1uarray(ulocate('qualEa'));       % activation energy (kcal/mole)
      Qoin =  F1uarray(ulocate('qualload'));     % quality load
      F1qual = [Dvalq;DtempCq;Zvalq;Eaq;Qoin];  % quality array 
 
%  ********************* MICROBIAL & QUALITY ASSESSMENT *********************

%     Pass microbial and quality assessment variables to exit streams
      E1mdk = Fmdk;  
      E1qual = Fqual;  

      if (mandq == 1) & (dlevel >= ?)   

%     Bulk temperature for microbial assessment (K)
      Tbulk = ?   

%     Minimum Residence Time (min)
      tmin = ?   
 
%     Final microbial load 
      Nout = mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,Tbulk,tmin)
 
%     Log reduction in microbial load
      DlogNred = log10(Noin/Nout);

%     Temperature for quality assessment (K)
      Tqual = ?;

%     Maximum processing time = Minimum processing time (min) 
      tqmax = ?;

%     Final quality load  
      Qout =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,tqmax); 

%     Log reduction in quality 
      DlogQred = log10(Qoin/Qout);

%     Exit stream 1 microbial and quality arrays
      E1mdk(size(F1mdk,1)) = Nout;
      Equal(size(F1qual,1)) = Qout;

%     Update system output array
      E1uarray = [E1me;E1mdk;E1qual];

      endtxt = [10 'END OF UNIT (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

      end  % microbial and quality assessment 

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(19:28); % economic analysis array 
      end  

%  ********************* SYSTEM OUTPUT ARRAY ASSIGNMENT *********************

%     Size stream info arrays to full system output u-arrays  
      E1uarray = ustd(E1uarray);
      E2uarray = ustd(E2uarray);

      sys =[E1uarray;E2uarray];    
 
%  ****************************** PRINTED OUTPUT ****************************

  if flag == 9    % print only after a solution is found 

%   Print results to specified file name; 
%    call filename manager function rfile.m 
    fid = rfile;

%   Name of unit operation; system input streams; system output streams  
    UNAME = Uname; 
    SINAMES = str2mat(F1name,F2name,Sname);
    SONAMES = str2mat(E1name,E2name);

    if (dlevel == 1)
%     System inputs array 
      SIARRAY = [F1mass F2mass ...]; 
%     System outputs array 
      SOARRAY = [E1mass E2mass ...];
    elseif (dlevel >= 2) 
%     System inputs array
      SIARRAY = [F1mass F2mass ... Smass];
%     System outputs array
      SOARRAY = [E1mass E2mass ...];
    end  

%   No intermediate streams
    SNNAMES = [];
    SNARRAY = []; 

%   Print mass balance results table
    msprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    if (dlevel >= 2) 
%     System inputs array  
      SIARRAY = [F1me F2me ... Sme]; 
%     System outputs array 
      SOARRAY = [E1me E2me ...];

%     Print energy balance results table
      eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    end   % energy balance tables 
    
%   Open the printable results file
    fid = fopen(fid,'a+');
    
    fprintf(fid,'UNIT MODEL INPUTS\n');     
    fprintf(fid,'     \n');
    fprintf(fid,'  First parameter                         %6.5f w/w \n',Ltsinit); 

    if (dlevel >= 2) 
      fprintf(fid,'  Second parameter                         %6.5f w/w \n',Ltsinit); 
    end

    if (dlevel >= 3) 
      fprintf(fid,'  Third parameter                      %5.3e kJ/h-m^2-K\n',U); 
    end  


    if (dlevel >= 4)   

    end  

    fprintf(fid,'     \n');
    fprintf(fid,'UNIT MODEL OUTPUTS\n');     
    fprintf(fid,'     \n');
    if (dlevel >=2) 
      fprintf(fid,'  Something                 %5.3e kJ/h\n',Q); 
    end  

    if (dlevel == 3)  
      fprintf(fid,'  Something                      %5.3e m^2\n',Area);
    end 

    if (dlevel >=4) 
    fprintf(fid,'     \n');
    
    endtxt = [10 'END OF (',Ename,') TRANSPORT RESULTS PRINT TO FILE' 10];  
    disp(endtxt);  

    end  

%   Microbial and Quality Assessment results printout 
    if (mandq == 1) & (dlevel >= 3)

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
      SIMQNAME = str2mat(F1name,F2name,...);
      SOMQNAME = str2mat(E1name,E2name,...); 

%     System input microbial array
      SIMARRAY = [F1mdk F2mdk ...];
%     System outputs microbial array
      SOMARRAY = [E1mdk;tmin;Tbulk-Tref;];

%     No intermediate streams
      SNMQNAME = [];
      SNMARRAY = [];
      SNQARRAY = [];

      mdkprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIMARRAY,SOMARRAY,SNMARRAY);

%     System input quality array
      SIQARRAY = [F1qual F2qual ...];
%     System outputs quality array
      SOQARRAY = [E1qual;tqmax;Tqual-Tref;];

      qprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIQARRAY,SOQARRAY,SNQARRAY);
    end 

    if (econ >= 1) 

      Uutil = [S W E NG FO CA]; 
      Usize = ?;
      Unum = 1; 
      Uecon = ECON;  
      uofile(Uutil,Usize,Unum,Uecon,Uname);

    end  


    fprintf(fid,'     \n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'     \n');

    % Close the printable results file
    fclose('all');

    sys = []; 
  end % flag 9 print 

else
   % Otherwise, no need to return anything since this is continuous system
   sys = [];

end
 
