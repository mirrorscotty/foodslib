      function [sys,x0] = evaps(t,x,u,flag,Ename,Fname,Sname,Vname,Lname,Cname)

%  EVAPS     Single effect evaporator model 
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for single effect evaporator steady-state 
%            model. Determines mass and energy balances, performs design at 
%            basic and advanced transport level. 
%   
%            Call: evaps(t,x,u,flag,Ename,Fname,Sname,Vname,Lname,Cname)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Ename - evaporator name
%             Fname - feed stream name 
%             Sname - steam stream name 
%             Vname - vapor stream generated by evaporator name  
%             Lname - Liquid concentrate stream name 
%             Cname - Steam condensate stream name   
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
%             Billet, R. 1989. "Evaporation Technology: Principles, 
%               Applications, Economics," Verlagsgesellschaft, mbH, Germany.


%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = 2;                            % Number of system input streams
SO = 3;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = 26;                     % Number of evaporator operating inputs 

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
      Suarray = u(usize+1:2*usize);  % full steam stream array 
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

%     Evaporator operation information 
      Lts = UOuarray(4);            % liquid concentrate total solids (w/w)  
      Ltsinit = Lts;

%  ***************************  MASS BALANCE LEVEL ************************** 

%     If feed stream mass flow rate is unknown, need to return
      if (F <= 0)  
       Vuarray(1:usize,1) = 0;   % no vapor stream  
       Luarray = Fuarray;           % send-out feed info; for system recycle 
       Cuarray(1:usize,1) = 0;   % no steam condensate stream  
       sys = [Vuarray;Luarray;Cuarray;];     
       return;
      end

%     Total solids content of feed (w/w); call tsolids.m function  
      Fts = tsolids(Fcomp,Fctype);  

%     User error message:  
      if (Lts <= Fts)  % Liquid concentrate total solids cannot be greater  
                       %  than feed total solids   
        errortxt = [10 'EVAPORATOR (',Ename,') ERROR: Specified total solids of' 10 ...
        'liquid concentrate is less than total solids of feed stream. ' 10];
        disp(errortxt);
        return;        % Stop execution of algorithm   
      end  
 
%     Liquid concentrate flow rate (kg/h)
      L = F*Fts/Lts;  

%     Overall Mass Balance - vapor flow rate (kg/h)
      V = F - L;  

%     Liquid Concentrate stream component types and composition (w/w) 
      Lctype = Fctype;            % All feed and liquid stream component types
                                  % are assumed to be the same    
      for N = 1:length(Fcomp)     % For each component, that is defined but
        if (Lctype(N) ~= comptype('water'))    %  is not water... 
          Lcomp(N,1) = Fcomp(N).*F./L;    % perform a component mass balance.  
        else Lcomp(N,1) = 0; 
        end  
        if (Lctype(N) == comptype('water')) xh2o = N; end % 
                                         %  If the component is water, note  
                                         %  index location in component array
      end 
      Lcomp(xh2o,1) = 1 - Lts;   % Set water component mass fraction       

%     Vapor stream stream component types and composition (w/w)
      Vctype = Fctype;            % All feed and liquid stream component types 
                                  % are assumed to be the same     
      Vcomp = abs((F.*Fcomp - L.*Lcomp)./V); % Component mass balances  

%     Exit stream mass arrays; used for print results table  
      Lmass = [L;Lcomp;Lctype];
      Vmass = [V;Vcomp;Vctype];

%     System output arrays 
      Luarray = Lmass;  
      Vuarray = Vmass;   
      Cuarray = 0; 

      endtxt = [10 'END OF EVAPORATOR (',Ename,') MASS BALANCE' 10];  
      disp(endtxt); 

%  **************** ENERGY BALANCE LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 2)  

%     Reference temperature (K)
      Tref = tref;

%     Feed stream information  
      FtempC = Fuarray(ulocate('temp'));% temperature (C) 
      FtempK = FtempC + Tref;            % temperature (K) 
      Fme = [Fmass;FtempK;Fuarray(ulocate('pressure'):ulocate('quality'))]; 
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
%     Evaporator operation information 
      Pin = UOuarray(5);                  % operating pressure (MPa)  
      Pininit = Pin;
      Cx = UOuarray(6);                   % steam condensate exit quality    
      Cxinit = Cx;
      BPRmodel = UOuarray(7);             % boiling point rise model 

%  *********************** ENERGY BALANCE LEVEL *****************************

%     Enthalpy of feed stream (kJ/kg(stream))
      Fh = enthalpy(Fme);   
 
%     Enthalpy of source steam (kJ/kg(steam)) 
      Sh = enthalpy(Sme);  

%     Assume boiling point of feed stream is at saturation temperature (K) 
%      associated with vapor space (operating) pressure (MPa)  
      Lpres = Pin;
      Tinsat = tpsteam(Pin); % Call to tpsteam.m function - 
                             %  temperature (K) given pressure at saturation 
%     Boiling Point Rise (K) 
      BPR = bpr(Lcomp,Lctype,Tinsat,BPRmodel); 

%     Liquid concentrate temperature 
      LtempK = Tinsat + BPR;  

      Lme = [Lmass;LtempK;Lpres;0];  % mass and energy information array 

%     Enthalpy of liquid condenstate stream (kJ/kg(stream))
      Lh = enthalpy(Lme);

%     Pressure (MPa) of vapor is at operating pressure
      Vpres = Pin;
%     Vapor temperature (K) is at saturation temp. for operating pressure 
      VtempK = Tinsat + BPR;      

%     Enthalpy of vapor leaving evaporator (kJ/kg(steam))
      Vx = 1;
      Vme = [Vmass;VtempK;Vpres;Vx]; % mass and energy information array  

%     Enthalpy of vapor stream (kJ/kg)  
      Vh = enthalpy(Vme);

%     If feed stream and vapor stream are not fully defined, need to return
      if (FtempK <= Tref) | (StempK <= Tref)
        Vuarray = ustd(Vme);    % vapor stream knowns 
        Luarray = ustd(Lme);    % liquid concentrate knowns  
        Cuarray(usize,1) = 0;              % no steam condensate stream
        sys = [Vuarray;Luarray;Cuarray;];
       return;   
      end  

%     Steam condensate component mass fractions (w/w) and types 
      Ccomp = Scomp;   
      Cctype = Sctype;    
      Cmass = Smass;               % mass information array

%     Steam condenstate outlet conditions where condensation of steam is a
%      constant pressure process (MPa)
      CtempK = tpsteam(Spres);     % Temperature (K)  
      Cpres = Spres;               % Pressure (MPa)   
      Cme = [Cmass;CtempK;Cpres;Cx]; % Mass and energy info array   

%     Enthalpy of outlet steam (kJ/kg(steam))
      Ch = enthalpy(Cme);

%     Energy Balances   
      if (S == 0)    % Case I - steam mass flow rate unknown 

        nosteam = 1;  % flag for not from steam source  

%       Steam requirement (kg/h)
        S = (L*Lh + V*Vh - F*Fh)/(Sh-Ch);     
        C = S; 

%       Heat transfer (kJ/h)
        Q = S*(Sh-Ch);
 
      else      % Case II - steam flow rate known 

        nosteam = 0;     % flag for steam not from source 

%       Steam condensate mass flow rate (kg/h)  
        C = S;

%       Maximum liquid concentrate total solids (w/w)
        Ltsmax = Lts;  
%       Minimum steam condensate exit quality 
        Cxmin = Cx;   
%       Maximum operating pressure (MPa)  
        Pinmax = Pin;  

%       Determination of liquid concentrate total solids (w/w)  
%       Bisection method 
        Ltslow = Fts;          % Lowest total solids is feed total solids 
        Ltshigh = Ltsmax;      % Highest total solids set by user  
        Qdiff = 1;  
        iter = 0;
        while (abs(Qdiff) > 1e-2)  
          iter = iter + 1;
          if iter == 50 break; end 
%         Liquid concentrate total solids (w/w)
          Lts = (Ltslow+Ltshigh)/2;
%         New liquid concentrate mass flow rate (kg/h)  
          L = F*Fts/Lts;  
%         New liquid Concentrate stream composition (w/w)
          for N = 1:length(Fcomp)     % For each component, that is defined but
            if (Lctype(N) > 1)        %  is not water...
              Lcomp(N,1) = Fcomp(N).*F./L;    % perform a component mass balance.
            else Lcomp(N,1) = 0;
            end
            if (Lctype(N) == 1) xh2o = N; end % If the component is water, note
                                          % index location in component array
          end
          Lcomp(xh2o,1) = 1 - Lts;     % Set water component mass fraction

%         Boiling Point Rise (K)
          BPR = bpr(Lcomp,Lctype,Tinsat,BPRmodel);  

%         Liquid concentrate temperature (K)  
          LtempK = Tinsat + BPR;

          Lmass = [L;Lcomp;Lctype];     % mass info array  
          Lme = [Lmass;LtempK;Lpres;0]; % mass and energy info array  
 
%         New enthalpy of exit stream (kJ/kg(stream))
          Lh = enthalpy(Lme);
 
%         Overall Mass Balance - vapor flow rate (kg/h)
          V = F - L;

%         Heat transfer (kJ/h)
          Q = S*(Sh-Ch);
          Qdiff = Q + F*Fh - L*Lh - V*Vh;

          if (Qdiff > 0)    % Still have heat available for evaporation 

            Ltslow = Lts;   % Raise lower limit of liq. conc. total solids 
%           If at max liq conc total solids; raise the steam cond. exit quality
            if (abs(Lts-Ltsmax) < 1e-6)    

              warntxt = [10 'EVAPORATOR (',Ename,') WARNING: Steam mass flow rate is too high' 10 ...
              'to achieve the maximum total solids and still condense steam to the specified quality.' 10];
              disp(warntxt);

              Cxhigh = 1;
              Cxlow = Cxmin;
              while (abs(Qdiff) > 1e-2)
%               Steam condensate quality
                Cx = (Cxhigh+Cxlow)/2;
%               Steam condensate enthalpy (kJ/kg)
                Ch = hsteam(0,Cpres,Cx);
 
%               Heat transfer (kJ/h)
                Q = S*(Sh-Ch);
                Qdiff = Q + F*Fh - L*Lh - V*Vh;
 
                if (Qdiff > 0) % Still have heat available for evaporation
                  Cxlow = Cx;  % Reset lower limit for steam cond. quality
                else           % Requires more heat to evaporate
                  Cxhigh = Cx; % Reset upper limit for steam cond quality
                end
              end % while Qdiff - Cx 
            end % if Lts-Ltsmax     
          else              % Requires more heat to evaporate  
            Ltshigh = Lts;  % Reset upper limit of liq. conc. total solids   
            if (abs(Lts-Fts) < 1e-6) 
              warntxt = [10 'EVAPORATOR (',Ename,') WARNING: Steam mass flow rate is insufficient' 10 ...
              'to raise feed stream temperature to the boiling point.' 10];
              disp(warntxt);

              Pinhigh = Pinmax;              % Highest operating pressure set by user (MPa)  
              Pinlow = ptsteam(FtempK);      % Lowest operating pressure (MPa)  
              while (abs(Qdiff) > 1e-3)  
%               Operating Pressure (MPa)    
                Pin = (Pinhigh+Pinlow)/2; 

%               Saturation temperature (K)  
                Tinsat = tpsteam(Pin); 
%               Boiling Point Rise (K)
                BPR = bpr(Lcomp,Lctype,Tinsat,BPRmodel);  
%               Liquid concentrate temperature
                LtempK = Tinsat + BPR;
%               Liquid mass and energy info array 
                Lme = [Lmass;LtempK;Lpres;0];    
%               Liquid concentrate enthalpy (kJ/kg) 
                Lh = enthalpy(Lme);

%               Vapor pressure (MPa) and temperature (K)  
                Vpres = Pin;
                VtempK = tpsteam(Pin);  
%               Vapor mass and energy info array 
                Vme = [Vmass;VtempK;Vpres;Vx];
%               Vapor enthalpy (kJ/kg)  
                Vh = enthalpy(Vme);   
 
%               Heat transfer (kJ/h)
                Q = S*(Sh-Ch);
                Qdiff = Q + F*Fh - L*Lh - V*Vh;  
 
                if (Qdiff > 0) % Still have heat available for evaporation
                  Pinlow = Pin;  % Reset lower limit for operating pressure 
                else           % Requires more heat to evaporate
                  Pinhigh = Pin; % Reset upper limit for operating pressure 
                end
              end % while Qdiff - Pin   
            end % if Lts-Fts  
          end % if Qdiff  

        end % while Qdiff  

      end  % no steam/steam mass flow rate calculations   

%     Update mass information arrays
      Smass = [S;Scomp;Sctype];
      Vmass = [V;Vcomp;Vctype];   
      Cmass = [C;Ccomp;Cctype];
%     Update mass and energy info arrays 
      Sme = [Smass;StempK;Spres;Sx]; 
      Vme = [Vmass;VtempK;Vpres;Vx]; 
      Cme = [Cmass;CtempK;Cpres;Cx]; 

%     System output arrays
      Luarray = Lme;
      Vuarray = Vme;
      Cuarray = Cme;

      endtxt = [10 'END OF EVAPORATOR (',Ename,') ENERGY BALANCE' 10];
      disp(endtxt);

      end % Energy Balance Level  

%  *************** BASIC TRANSPORT LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 3)

%     Evaporator operation information
      U = UOuarray(8);       % Overall heat transfer coefficient (kJ/h-m^2-K)  
      Din = UOuarray(9);     % inner pipe diameter (m)
      vel = UOuarray(10);    % average velocity per tube (m/s)
      fthick = UOuarray(11);   % film thickness (m)

      end % Basic Transport Level Assignments  

%  ************************** BASIC TRANSPORT LEVEL *************************

      if (dlevel >= 3) 

%     Heat Transfer Area (m^2)
      Area = Q/(U*(StempK - LtempK));

%     Flow cross sectional area per tube (m^2)
      Ax = pi*(Din/2)^2 - pi*((Din-fthick)/2)^2;  

%     Feed density (kg/m^3)
      Frho = rho(Fcomp,Fctype,FtempK); 

%     Volumetric flow rate (m^3/h) 
      Vq = F/Frho;

%     Number of tubes 
      Nt = Vq/(Ax*vel*3600); 

%     Length of tubes (m) 
      length = Area/(Nt*pi*Din); 

%     Time (min)
      tmin = Nt*Ax*length/Vq*60;  

      endtxt = [10 'END OF EVAPORATOR (',Ename,') BASIC TRANSPORT' 10];
      disp(endtxt);

      end  % Basic Transport Level 

%  ************* ADVANCED TRANSPORT LEVEL VARIABLE ASSIGNMENTS **************

      if (dlevel >= 4) 

%     Evaporator operation inforamtion 
      evaptype = UOuarray(12);  % evaporator model   
      vertical = UOuarray(13);  % tube orientation; 
                               %  vertical (1); horizontal (0)
      mumodel = UOuarray(14);   % viscosity model 
      wallthk = UOuarray(15);  % pipe thickness (m)
      Do = Din + 2*wallthk;    % outer pipe diameter (m)
      kpipe = UOuarray(16);    % pipe thermal conductivity (kJ/h-m-K)

%  ************************ ADVANCED TRANSPORT LEVEL ************************

%     Boiling temperature (K) 
      Tin = Tinsat + BPR;   

%     Set first guess at recirculation mix in evap mass and 
%      energy information array 
      Mme = Fme;
      Mmenew = Fme;
      Mmenew(1) = Mme(1)*1e5;

%     Looking for recirculation mix mass flow rate to be within a specified
%      range for two consectutive iterations 
      J = 0;                       % index for tracking number of iterations 
      while (abs(Mme(1) - Mmenew(1)) > 1e4) 
        J = J+1;  
        if J == 30 break; end      % break from loop after 30 iterations 

%       Reset recirculation mix mass flow rate for next iteration 
        if (J > 1) Mme = Mmenew; end 

        lengthdif = 2;
        while (lengthdif > 1)  

%       Trial & Error solution for inside and outside heat transfer 
%        coefficients
%       Looking for steam/vapor side wall temp guess to equal calculated value
%       Guess at Temperature at steam side tube wall (K)
        Twvguess = linspace(StempK-0.1,Tin+0.1,100);
        for N = 1:100   

%         Outside heat transfer coefficient (kJ/h-m^2-K)
          hoarray = hostcond(StempK,Spres,Twvguess(N),Do,length,vertical,Nt);
          houtN(N) = hoarray(1);
          NReoN(N) = hoarray(2);  
          NNuoN(N) = hoarray(3);
 
%         Temperature at product side tube wall (K)
          TwlN(N) = Twvguess(N) - houtN(N)*wallthk/kpipe*(StempK-Twvguess(N));

%         Inside heat transfer coefficient (kJ/h-m^2-K)
          evapdesc = [Din fthick vel length TwlN(N) evaptype mumodel]; 
          hiarray = hievap(Mme,Lme,evapdesc);
          hinN(N) = hiarray(1);  
          NReiN(N) = hiarray(2); 
          NPriN(N) = hiarray(3);
          NNuiN(N) = hiarray(4);  
 
%         Calculated temperature at steam side tube wall (K)
          Twvnew(N) = ((houtN(N)/hinN(N))*StempK+ ... 
                       houtN(N)*wallthk/kpipe*StempK+Tin)/ ...
                     (1+houtN(N)*wallthk/kpipe+houtN(N)/hinN(N));

          Twvdiff(N) = Twvguess(N) - Twvnew(N);
          if (Twvdiff(N) < 0) break; end  

        end % while loop for wall temperature 

%       Locate index for least difference between calculated value 
%        and guessed steam side wall temperature
        Twvdiffmin = min(abs(Twvdiff(N)),abs(Twvdiff(N-1)));
        index = findval1(abs(Twvdiff),Twvdiffmin);
%       Use index to find temperature profile information 
        Twv = Twvguess(index);
        Twl = TwlN(index);   
        hout = houtN(index);
        NReo = NReoN(index);  
        NNuo = NNuoN(index);  
        hin =  hinN(index);   
        NRei = NReiN(index);
        NPri = NPriN(index); 
        NNui = NNuiN(index);

%       Outside Area per tube (m^2)
        Aoutt = pi*length*Do;

%       Overall Heat Transfer Coefficient (kJ/h-m^2-K)
        U = 1/(1/hout + wallthk/kpipe + 1/hin);
 
%       Effective area required for heat transfer (m^2)
        Area = Q/(U*(StempK-Tin));    % Area for costing 

%     Volumetric flow rate (m^3/h)
      Vq = Mme(1)/Frho; 

%     Number of tubes
      Nt = Vq/(Ax*vel*3600);  
 
%     Length of tubes (m)
      lengthold = length;
      length = Area/(Nt*pi*Din); 
      lengthdif = abs(lengthold-length);
      
      end  


%       Initialize recycle ratio
        recyrat = 0;
        R = 0;

        if evaptype <= 20 
%         Recirculation mass flow to achieve minium velocity  
          if (Nt < 10) | (length > 20) Nt = Nt + 10; end  
          M = Nt*vel*Ax*Frho*3600;  
          R = M - F;              % recirculation flow rate
          recyrat = R/L;          % recirculation ratio
%         Component mass fractions (w/w) and types for recirculation 
          Mcomp = (F.*Fcomp + R.*Lcomp)./M;  
          Mctype = Fctype;  
%         Recirculation temperature (K)
          MtempK = mixtemp([Fme Lme]);  
%         Recirculation mass and energy information array  
          Mmenew = [M;Mcomp;Mctype;MtempK;Pin;0];  
        end    

      end % while loop for recirculation 

%     Vapor space requirement
%     Vapor load (m^3/m^3-h) at 1 bar
      vc = 7300;
 
%     Adjustment factor for vapor spaces at pressure other than 1 bar
      f = [1.35; 1.01; 0.93; 0.91; 0.885; 0.875; 0.865; 0.86; 0.855; 0.85;];
      P = [0.05; 0.1; 0.15; 0.2; 0.25; 0.3; 0.35; 0.4; 0.45; 0.5;];
      if (Pin < min(P))
        factor = max(f);  
      elseif (Pin > max(P))
        factor = min(f); 
      else 
        factor = interp1(P,f,Pin);
      end  
      vcPin = factor*vc;
 
%     Vapor velocity (m/s): assumes at 1 bar that velocity is 1 m/s
      rhov1 = rhosteam(0,0.1,1);
      rhovPin = rhosteam(0,Pin,Vx);
      velvPin = sqrt(rhov1/rhovPin);
 
%     Vapor space capacity (m^3)
      VB = V/(rhovPin*vcPin);
 
%     Diameter of vapor space (m)
      dB = sqrt(4*V/(3600*pi*velvPin*rhovPin));
 
%     Height of vapor space (m)
      Hv = 4*V/(pi*dB^2*rhovPin*vcPin);

      endtxt = [10 'END OF EVAPORATOR (',Ename,') ADVANCED TRANSPORT' 10]; 
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

%     Pass microbial and quality assessment variables to liq. conc stream
      Lmdk = Fmdk;
      Lqual = Fqual;

      if (mandq == 1) & (dlevel >= 3)   

%     Volume hold up in tubes (design is for 100% full)
%     if (dlevel == 3)
%       Vhold = Nt*pi*(Din/2)^2*length
%     else  
%       if (Din <= 2*fthick) 
%         Vhold = Nt*pi*(Din/2)^2*length
%       else  
          Vhold = Nt*Ax*length;
%       end 
%     end 
%     Bulk temperature on product side (K)
      Tbulk = (FtempK+LtempK)/2;
%     Bulk density of concentrated product (kg/m^3)
      rhobulk = rho(Fcomp,Fctype,Tbulk);
%     Volumetric liquid flow rate
      Lvol = F/rhobulk;
%     Minimum Residence Time (min)
      tmin = Vhold/Lvol*60;
 
%     Final microbial load 
      Nout = mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,Tbulk,tmin);  
 
%     Log reduction in microbial load
      DlogNred = log10(Noin/Nout);

      if (dlevel == 3)   
%       Temperature for quality assessment taken as liquid conc. exit temp (K)
        Tqual = LtempK;
%       Maximum processing time = Minimum processing time (min) 
        tqmax = tmin;
%       Final quality for single pass at liq. conc. boiling point  
        Qout =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,tqmax); 

      elseif (dlevel >= 4) & (Mme(ulocate('mass')) > Fme(ulocate('mass'))) 
%       Temperature for quality assessment taken as liquid conc. wall temp (K) 
        if (NRei < 2100) 
          Tqual = Twl;
        else 
          Tqual = Tbulk;
        end 
%       Use recycle ratio to find time for 50 per. and 90 per. leaving  
        Fnew = F;  
        X = 0;  
        first = 0;
        notrecy = 1/(recyrat+1);  
 
        while (Fnew > 0.1*F) 
          X = X+1;
          Fnew = Fnew - notrecy*Fnew; 
          if Fnew <= 0.5*F & first == 0  
            first = 1; 
            J = X;
          end 
        end 
 
        t50 = tmin*J;
        tqmax = tmin*X; 
 
%       Final quality for time when 50% of feed leaves the evaporator  
        Qout50 =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,t50); 
%       Final quality for time when 90% of feed leaves the evaporator  
        Qout =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,tqmax);  

      elseif (dlevel >= 4) & (Mme(ulocate('mass')) == Fme(ulocate('mass'))) 
%       Temperature for quality assessment taken as liquid conc. wall temp (K)  
        if (NRei < 2100)
          Tqual = Twl;
        else
          Tqual = Tbulk;
        end
%       Maximum processing time = Minimum processing time (min)
        tqmax = tmin;
%       Final quality for single pass evaporator  
        Qout =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,tqmax); 
      end 

%     Log reduction in quality 
      DlogQred = log10(Qoin/Qout);

%     Liquid concentrate microbial and quality arrays
      Lmdk = Fmdk;  
      Lmdk(size(Fmdk,1)) = Nout;
      Lqual = Fqual;
      Lqual(size(Fqual,1)) = Qout;

      endtxt = [10 'END OF EVAPORATOR (',Ename,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

      end  % microbial and quality assessment 

%     Update system output array
      Luarray = [Lme;Lmdk;Lqual];

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECONarray = UOuarray(17:26); % economic analysis array 
      end  

%  ********************* SYSTEM OUTPUT ARRAY ASSIGNMENT *********************

%     Size stream info arrays to full system output u-arrays  
      Vuarray = ustd(Vuarray);
      Luarray = ustd(Luarray);
      Cuarray = ustd(Cuarray);

      sys =[Vuarray;Luarray;Cuarray];    
 
%  ****************************** PRINTED OUTPUT ****************************

  if flag == 9    % print only after a solution is found 

%   Print results to specified file name; 
%    call filename manager function rfile.m 
    fid = rfile;
    
%   Name of unit operation; system input streams; system output streams  
    UNAME = Ename; 
    SINAMES = str2mat(Fname,Sname);
    SONAMES = str2mat(Lname,Vname,Cname);

    if (dlevel == 1)
%     System inputs array 
      SIARRAY = [Fmass]; 
%     System outputs array 
      SOARRAY = [Lmass Vmass];
    elseif (dlevel >= 2) 
%     System inputs array
      SIARRAY = [Fmass Smass];
%     System outputs array
      SOARRAY = [Lmass Vmass Cmass];
    end  

%   No intermediate streams
    SNNAMES = [];
    SNARRAY = []; 

%   Print mass balance results table
    msprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 
    
    if (dlevel >= 2) 
%     System inputs array  
      SIARRAY = [Fme Sme]; 
%     System outputs array 
      SOARRAY = [Lme Vme Cme];

%     Print energy balance results table
      eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 
      
    end  
    % THIS IS WHERE THE CHANGE WAS MADE 4/21/01
    fid = fopen(fid,'a+');
    % END CHANGE 4/21/01
    
    fprintf(fid,'EVAPORATOR MODEL INPUTS\n');     
    fprintf(fid,'     \n');
    fprintf(fid,'  Liquid Concentrate Total Solids      %6.5f w/w \n',Ltsinit); 

    if (dlevel >= 2) 
      fprintf(fid,'  Operating Pressure                     %5.3f MPa \n',Pininit);
      fprintf(fid,'  Steam Condensate Quality               %5.3f \n',Cxinit); 
      fprintf(fid,'  Boiling Point Rise Model               %5.0f \n',BPRmodel);     
    end

    if (dlevel == 3) 
      fprintf(fid,'  Overall Heat Transfer Coefficient  %5.3e kJ/h-m^2-K\n',U); 
    end  

    if (dlevel >=3) 
 %    fprintf(fid,'  Inner Tube Diameter                   %5.4f m \n',Din);
 %    fprintf(fid,'  Minimum Average Velocity               %5.2f m/s \n',vel);
 %    fprintf(fid,'  Film Thickness                        %5.4f m \n',fthick);
    end 
  
    if (dlevel >= 4)   
      fprintf(fid,'  Evaporator Model                       %5.0f \n',evaptype);
      fprintf(fid,'  Tube Orientation                       %5.0f \n',vertical);
      fprintf(fid,'  Viscosity Model                        %5.0f \n',mumodel);
      fprintf(fid,'  Wall Thickness                        %5.4f m \n',wallthk);
      fprintf(fid,'  Pipe Thermal Conductivity         %5.4e kJ/h-m-K \n',kpipe);
      fprintf(fid,'  Minimum Average Velocity               %5.2f m/s \n',vel);

    end  

    if (dlevel >=2) 
    fprintf(fid,'     \n');
    fprintf(fid,'EVAPORATOR MODEL OUTPUTS\n');     
    fprintf(fid,'     \n');
      fprintf(fid,'  Boiling Point Rise                     %5.3f K \n',BPR); 
      fprintf(fid,'  Steam Economy                          %5.3f \n',V/S); 
      fprintf(fid,'  Heat Transfer Rate                 %5.3e kJ/h\n',Q); 
    end  

    if (dlevel == 3)  
      fprintf(fid,'  Heat Transfer Area                 %5.3e m^2\n',Area);
 %    fprintf(fid,'  Tube Length                            %5.2f m \n',length);
 %    fprintf(fid,'  Number of Tubes                        %5.0f \n',Nt);
    end 

    if (dlevel >=4) 
    fprintf(fid,'     \n');
    fprintf(fid,'  Evaporator Detailed Design:\n');
    fprintf(fid,'    Tube Length                          %5.2f m \n',length);
    fprintf(fid,'    Number of Tubes                      %5.0f \n',Nt);
    fprintf(fid,'    Heat Transfer Area                   %5.1f m^2\n',Area);
    fprintf(fid,'     \n');
    if evaptype <=  20 
      fprintf(fid,'    Recirculation Ratio                  %5.2f \n',recyrat);
      fprintf(fid,'    Recirculation Flow Rate           %5.2e kg/h \n',R);
    end 
    fprintf(fid,'     \n');
    fprintf(fid,'   Inside tube conditions  \n');
    fprintf(fid,'    Heat Transfer Coefficient       %5.4e kJ/h-m^2-K\n',hin);
    fprintf(fid,'    Reynolds Number                 %5.4e \n',NRei);
    fprintf(fid,'    Prandtl Number                  %5.4e \n',NPri);
    fprintf(fid,'    Nusselt Number Or Dimensionless \n');
    fprintf(fid,'     Heat Transfer Coefficient      %5.4e \n',NNui);
    fprintf(fid,'     \n');
    fprintf(fid,'   Heat Conduction Through Wall     %5.4e kJ/h-m^2-K\n',kpipe/wallthk);
    fprintf(fid,'     \n');
    fprintf(fid,'   Outside tube conditions  \n');
    fprintf(fid,'    Heat Transfer Coefficient       %5.4e kJ/h-m^2-K\n',hout);
    if (vertical == 1)
      fprintf(fid,'    Reynolds Number                 %5.4e \n',NReo);
    end   
    fprintf(fid,'    Nusselt Number                  %5.4e \n',NNuo);
    fprintf(fid,'     \n');
    fprintf(fid,'   O.H.T. Coefficient               %5.4e kJ/h-m^2-K\n',U);
    fprintf(fid,'     \n');
    fprintf(fid,'   Inside Wall Temperature              %5.2f C \n',Twl-Tref);
    fprintf(fid,'   Outside Wall Temperature             %5.2f C \n',Twv-Tref);
    fprintf(fid,'     \n');
    fprintf(fid,'  Vapor Space:   \n');
    fprintf(fid,'   Capacity                              %5.2f m^3 \n',VB);
    fprintf(fid,'   Height                                %5.2f m \n',dB);
    fprintf(fid,'   Diameter                              %5.2f m \n',Hv);
    
    endtxt = [10 'END OF (',Ename,') TRANSPORT RESULTS PRINT TO FILE' 10];  
    disp(endtxt);  
    disp('PRINTED TRANSPORT RESULTS') 
    end  

%   Microbial and Quality Assessment results printout 
    if (mandq == 1) & (dlevel >= 3)

    fprintf(fid,'     \n');
    fprintf(fid,'  Microbial (Worst Case) Assessment:   \n');
    fprintf(fid,'   Single Pass (minimum) residence time %5.4f min \n',tmin);
    fprintf(fid,'   Minimum Processing Temperature        %5.2f C \n',Tbulk-Tref);
    fprintf(fid,'   Log Reduction in Microbial Load       %5.2e \n',DlogNred);
    fprintf(fid,'     \n'); 
    fprintf(fid,'  Quality (Worst Case) Assessment:   \n');
    if (dlevel == 3) 
    fprintf(fid,'   Single Pass (minimum) residence time %5.4f min \n',tqmax);
    fprintf(fid,'   Maximum Processing Temperature        %5.2f C \n',Tqual-Tref);  
    elseif (dlevel >= 4) 
      if (Mme(ulocate('mass')) > Fme(ulocate('mass')))
    fprintf(fid,'   90 Percent Residence Time            %5.4f min \n',tqmax);
      else 
    fprintf(fid,'   Single Pass Residence Time           %5.4f min \n',tqmax);
      end  
    fprintf(fid,'   Maximum Processing Temperature       %5.2f C \n',Tqual-Tref);
    end   
    fprintf(fid,'   Log Reduction in Quality Load        %5.2e \n',DlogQred);



%     System input and exit stream names for microbial and quality assessment 
      SIMQNAME = str2mat(Fname);
      SOMQNAME = str2mat(Lname); 

%     System input microbial array
      SIMARRAY = [Fmdk];
%     System outputs microbial array
      SOMARRAY = [Lmdk;tmin;Tbulk-Tref;];

%     No intermediate streams
      SNMQNAME = [];
      SNMARRAY = [];
      SNQARRAY = [];

      mdkprint(Ename,SIMQNAME,SOMQNAME,SNMQNAME,SIMARRAY,SOMARRAY,SNMARRAY);

%     System input quality array
      SIQARRAY = [Fqual];
%     System outputs quality array
      SOQARRAY = [Lqual;tqmax;Tqual-Tref;];

      qprint(Ename,SIMQNAME,SOMQNAME,SNMQNAME,SIQARRAY,SOQARRAY,SNQARRAY);
    end 

    if (econ >= 1) 

      if (nosteam == 1) & (dlevel >= 2)  
        Eutil = [S 0 0 0 0 0]; 
      else Eutil = [0 0 0 0 0 0]; end  
      if (dlevel >= 3) 
        Esize = Area;
        Enum = 1; 
        Eecon = ECONarray;  
      else  
        Esize = 0;   
        Enum = 0; 
        Eecon = 0;  
      end  
      uofile(Eutil,Esize,Enum,Eecon,Ename);

    end  


    fprintf(fid,'     \n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'     \n');
    
    % THIS IS WHERE THE CHANGE WAS MADE 4/21/01
    fclose('all');
    % END CHANGE 4/21/01
    
    sys = []; 
  end % flag 9 print 

else
   % Otherwise, no need to return anything since this is continuous system
   sys = [];

end
 
