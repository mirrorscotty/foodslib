      function [sys,x0] = spdrys(t,x,u,flag,Uname,Fname,AAname,Sname,Pname,AOname,Cname)

%  SPDRYS    Spray dryer model   
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for spray dryer steady-state model. 
%            Determines mass and energy balances, performs design at 
%            basic and advanced transport level. 
%   
%            Call: unitops(t,x,u,flag,Uname,Fname,AAname,Sname,Pname,AOname,Cname)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Uname - spray dryer name
%             Fname - feed stream name 
%             Aname - air stream name 
%             Sname - water/steam stream name 
%             Pname - dry product stream name 
%             Cname - water/steam condenstate stream name 
%  
%            Returned Arguments: 
%             sys - SIMULINK S-function variable: for flag = 0, 
%                   sys = sizes of parameters and initial conditions. 
%                   For flag = 3, sys = system outputs which are the 
%                   u-arrays for each stream exiting the unit operation 
%
%            References:  
%             Heldman, D.R. and Singh, R.P.  1981. "Food Process
%              Engineering," 2nd ed. Van Nostrand Reinhold, NY
%             Masters, K. 1985. "Spray Drying Handbook," 4th ed. 
%              John Wiley & Sons, NY   

%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = 3;                            % Number of system input streams
SO = 3;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = 30;                      % Number of unit operating inputs 

if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, (SO*usize), ((SI*usize)+UOinputs), 0, 1];
   x0 = ones(size(u));

elseif flag == 3 | flag == 9
   % if flag=3, return system outputs
  
%     Breakup of SIMULINK constructed u-array into stream arrays and 
%      operating information array  
 
      Fuarray = u(1:usize);            % full feed stream array 
      AAuarray = u(usize+1:2*usize);   % full ambient air stream array 
      Suarray = u(2*usize+1:3*usize);  % full water/steam stream array 
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

%     Spray dryer information for mass balances  
      Pts = UOuarray(4);            % product total solids (w/w)   

%  ***************************  MASS BALANCE LEVEL ************************** 

%     Check for recycle problems 
      if (F <= 0) 
        AOuarray = 0;
        AOuarray = ustd(AOuarray);
        Cuarray = 0; 
        Cuarray = ustd(Cuarray); 
        sys = [Fuarray;AOuarray;Suarray];
        return; 
      end  

%     Feed total solids (w/w) 
      Fts = tsolids(Fcomp,Fctype); 
      
%     Check that total solid of feed is less than total solid of exit
      if (Fts > Pts)
        errortxt = [10 'SPRAY DRYER (',Uname,') ERROR: Total solids of feed is greater than ' 10 ... 
        'desired total solids of product.' 10];  
        disp(errortxt);  % Stop execution of algorithm
        return;
      end

%     Dry product flow rate (kg/h)
      P = Fts*F/Pts;
 
%     Overall Mass Balance - Water removal from product (kg/h)
      H = F - P;

%     Dry Product stream component types and composition (w/w)
      Pctype = Fctype;          % All feed and product stream component types
                                % are assumed to be the same
      for N = 1:length(Fcomp)     % For each component, that is defined but
        if (Pctype(N) ~= comptype('water'))  %  is not water...
          Pcomp(N,1) = Fcomp(N).*F./P;    % perform a component mass balance.
        else Pcomp(N,1) = 0;
        end 
        %  If the component is water, note index location in component array
        if (Pctype(N) == comptype('water')) xh2o = N; end % 
      end
      Pcomp(xh2o,1) = 1 - Pts;   % Set water component mass fraction
         
%     Water removed from feed stream component types and composition (w/w)
      Hctype = Fctype;            % All feed and water stream component types
                                  % are assumed to be the same
      Hcomp = abs((F.*Fcomp - P.*Pcomp)./H); % Component mass balances

%     Exit stream mass arrays; used for print results table  
      Pmass = [P;Pcomp;Pctype];
      Hmass = [H;Hcomp;Hctype];  

%     System output arrays 
      Puarray = Pmass;  
      AOuarray = 0;   
      Cuarray = 0;

      endtxt = [10 'END OF SPRAY DRYER (',Uname,') MASS BALANCE' 10];  
      disp(endtxt); 

%  **************** ENERGY BALANCE LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 2)  

%     Reference temperature (K)
      Tref = tref;

%     Feed stream information  
      FtempC = Fuarray(ulocate('temp'));% temperature (C) 
      FtempK = FtempC + Tref;            % temperature (K) 
      Fpres = Fuarray(ulocate('pressure')); % pressure (MPa)  
      Fx = Fuarray(ulocate('quality')); 
      Fme = [Fmass;FtempK;Fpres;Fx]; % mass and energy information array  

%     Ambient Air stream information
      AA = AAuarray(ulocate('mass'));      % mass flow rate (kg/h)
      AAcomp = AAuarray(ulocate('comp1st'):ulocate('complast'));
                                         % component mass fractions (w/w)
      AActype = AAuarray(ulocate('ctype1st'):ulocate('ctypelast'));
                                         % component types
      AAmass = [AA;AAcomp;AActype];          % mass information array
      AAtempC = AAuarray(ulocate('temp'));% temperature (C)
      AAtempK = AAtempC + Tref;            % temperature (K)
      AApres = AAuarray(ulocate('pressure')); % pressure (MPa) 
      AAx = AAuarray(ulocate('quality'));  
      AAme = [AAmass;AAtempK;AApres;AAx];
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

%     Spray Dryer information for energy balances  
      PtempC = UOuarray(5);     % product temperature (C)
      PtempK = PtempC + tref;  
      AItempC = UOuarray(6);    % air inlet temperature (C)
      AItempK = AItempC + tref;  
      AOtempC = UOuarray(7);    % air outlet temperature (C)
      AOtempK = AOtempC + tref;  
      AIpres = UOuarray(8);     % operating pressure (MPa)
      Ql = UOuarray(9);   % heat loss rate from dryer (kJ/h)
      frecycle = UOuarray(10);   % air recycle fraction (w/w)

%  *********************** ENERGY BALANCE LEVEL *****************************

%     Check for recycle problems 
      if (FtempK <= Tref) | (StempK <= Tref) | (AAtempK <= Tref)
        sys = [Puarray;AAuarray;Suarray];
        return;
      end 
 
%     Simultaneous solution of water balances on heat exchanger and
%     spray drier and energy balance on spray drier. Cocurrent configuration.

%     Enthalpy of feed (kJ/kg)
      Fh = enthalpy(Fme);
%     Total energy of feed (kJ/h)
      Fhtot = F*Fh;

%     Dry product mass and energy information array 
      Ppres = AIpres;
      Px = Fx;
      Pme = [Pmass;PtempK;Ppres;Px];

%     Enthalpy of product (kJ/kg)
      Ph = enthalpy(Pme);
%     Total energy of product (kJ/h)
      Phtot = P*Ph;

%     Water and air component fractions (w/w) 
      AAxh2o = findcomp(AAcomp,AActype,comptype('water'));  
%     Absolute humidity of ambient air (kg water/kg dry air)
      AAw = AAxh2o/(1-AAxh2o);
%     Relative humidity of ambient air (percent)
      AARH = airwater(4,AAtempK,AApres,4,AAw);  

%     Enthalpy constants for dryer air inlet
      m1 = 2500.45 + 1.805*(AItempK-Tref);
%     m1 = 2501 + 1.805*(AItempK-Tref);
      m2 = 1.006*(AItempK-Tref);
%     Enthalpy constants for dryer air outlet
      n1 = 2500.45 + 1.805*(AOtempK-Tref);
%     n1 = 2501 + 1.805*(AOtempK-Tref);
      n2 = 1.006*(AOtempK-Tref);

%     Heat loss through chamber wall (kJ/h)   
%     Ql = 0;
 
%     Solution of simultaneous M&E balances
      Qtotal = -Fhtot + Phtot + Ql;
      Ganum = Qtotal + (1+frecycle)*n1*H - frecycle*m1*H;
      Gaden = (1+frecycle)*(AAw*(m1-n1) + m2 - n2);
%     Mass flow rate of air (kg dry air/h)
      Ga = Ganum/Gaden;
 
%     Absolute humidity of dryer inlet air (kg water/kg dry air)
      AIw = (Ga*AAw+frecycle*(Ga*AAw+H))/((1+frecycle)*Ga);
%     Enthalpy of dryer inlet air (kJ/kg-dry air)
      AIh = airwater(1,AItempK,AIpres,4,AIw);
%     Relative humidity of dryer inlet air
      AIRH = airwater(4,AItempK,AIpres,4,AIw);
%     Absolute humidity of outlet air (kg water/kg dry air)
      AOw = H/Ga + AAw;
%     Enthalpy of dryer outlet air (kJ/h)
      AOh = airwater(1,AOtempK,AIpres,4,AOw);
%     Relative humidity of dryer outlet air
      AORH = airwater(4,AOtempK,AIpres,4,AOw);

%     Total energy of air inlet to heater (kJ/h)
      AIhtot = (1+frecycle)*Ga*AIh;
%     Total energy of air outlet from heater (kJ/h)
      AOhtot = AIhtot - Qtotal;
 
%     Total air-water mass flow into heat exchanger (kg/h)
      AA = Ga + Ga*AAw;
 
%     Total air-water mass flow into spray dryer (kg/h)
      AI = AA + frecycle*Ga*(1+AOw);
%     Water fraction in total air (w/w)
      AIxh2o = AIw/(1+AIw);
%     Component fractions (w/w) and types 
      AIctype = AActype;
      AIcomp = zeros(size(AAcomp,1),1);
      xh2o = findval1(AIctype,comptype('water'))
      AIcomp(xh2o) = AIxh2o;
      xair = findval1(AIctype,comptype('air'))
      AIcomp(xair) = 1-AIxh2o;  
 
%     Total air-water mass flow out spray dryer (kg/h)
      AO = AA + H;
%     Water fraction in total air (w/w)
      AOxh2o = AOw/(1+AOw);  
%     Component fractions (w/w) and types  
      AOctype = AIctype;  
      AOcomp = zeros(size(AIcomp,1),1);  
      xh2o = findval1(AOctype,comptype('water')) 
      AOcomp(xh2o) = AOxh2o; 
      xair = findval1(AOctype,comptype('air')) 
      AOcomp(xair) = 1-AOxh2o;

%     Enthalpy of ambient air (kJ/kg dry air)
      AAh = airwater(1,AAtempK,AApres,4,AAw);
%     Enthalpy of air into heater (kJ/kg dry air)
      AAmixh = (AAh+frecycle*AOh)/(1+frecycle);
%     Approximate air temp (K) into heater
      AIHXtempK = (AA*AAtempK+frecycle*Ga*(1+AOw)*AOtempK)/ ...
                   (AA+frecycle*Ga*(1+AOw)); 
 
%     Energy Balance on indirect (steam) air heater
%     Total energy change on air side (kJ/h)
      Qairhx = (1+frecycle)*Ga*(AIh - AAmixh);
%     Qairhx = Ga*(AIh - AAh);
 
%     Enthalpy of entering steam (kJ/kg)
      Sh = hsteam(StempK,Spres,Sx);
%     Water/Steam Condensate exit conditions  
      Ccomp = Scomp;
      Cctype = Sctype;  
      CtempK = StempK;
      Cpres = Spres;  
      Cx = 0; 
%     Water/Steam condensate enthalpy (kJ/kg);
      Ch = hsteam(CtempK,Cpres,Cx);
 
%     Steam mass flow rate (kg/h)
      S = Qairhx/(Sh - Ch);
      C = S;

%     Update mass information arrays
      Smass = [S;Scomp;Sctype];
      Cmass = [C;Ccomp;Cctype];
      AAmass = [AA;AAcomp;AActype];   
      AOmass = [AO;AOcomp;AOctype];
%     Intermediate Streams
      ARmass = [frecycle*Ga*(1+AOw);AOcomp;AOctype]; 
      AIHXmass = [AI;AIcomp;AIctype];  
      AImass = [AI;AIcomp;AIctype];  
      AOSDmass = [AI+H;AOcomp;AOctype];  
       

%     Update mass and energy info arrays 
      Sme = [Smass;StempK;Spres;Sx]; 
      Cme = [Cmass;CtempK;Cpres;Cx]; 
      AAme = [AAmass;AAtempK;AApres;AAx]; 
      AOme = [AOmass;AOtempK;AIpres;AAx]; 
%     Intermediate Streams 
      ARme = [ARmass;AOtempK;AIpres;AAx]; 
      AIHXme = [AIHXmass;AIHXtempK;AIpres;AAx]; 
      AIme = [AImass;AItempK;AIpres;AAx]; 
      AOSDme = [AOSDmass;AOtempK;AIpres;AAx]; 

%     System output arrays
      Puarray = Pme;
      AOuarray = AOme;
      Cuarray = Cme; 
 
      endtxt = [10 'END OF SPRAY DRYER (',Uname,') ENERGY BALANCE' 10];
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

%     Spray Dryer inforamtion for advanced transport  

      diammean = UOuarray(11);   % mean droplet diameter (m)
      whRPM = UOuarray(12);      % wheel speed (m)
      whdiam = UOuarray(13);     % wheel diameter (m)
      whvanes = UOuarray(14);    % number of vanes
      vheight = UOuarray(15);    % vane height (m)
      FCcrit = UOuarray(16);     % critical moisture content (kg-water/kg-total)
      Prho = UOuarray(17);       % product density (w/w)
      htod = UOuarray(18);      % cyclindrical height to diameter ratio 
      phideg = UOuarray(19);     % cone angle (deg. from horizontal)
      phi = phideg*pi/(180);     % cone angle (rad. from horizontal)
      mumodel = UOuarray(20);    % viscosity model  

%  ************************ ADVANCED TRANSPORT LEVEL ************************

%     Determination of mean (50%) droplet diameter (m)
      if (diammean <= 0)
%       Atomizer wheel with vanes
        diammean = 1.4e4*(F)^0.24/ ...
                   ((whRPM*whdiam)^0.83*(whvanes*vheight)^0.12)/1e6;
      end
%     Determination of large (95%) droplet diameter (m)
      if (diammean < 60e-6)
        diam95 = 2*diammean;
      elseif (diammean >= 60e-6) & (diammean < 120e-6)
        diam95 = 2.5*diammean;
      else
        disp('SPRAY DRYER WARNING: Mean droplet diameter is outside the size range for a rotary vaned wheel atomizer. Standard maximum diameter is 120 microns.');
        diam95 = diammean; % made-up this criteria  
      end


%     Power to drive atomizer wheel (kW)
      power = 3.8e-10*F*whRPM^2*2*whdiam^2;

%     Feed density (kg/m^3)
      Frho = rho(Fcomp,Fctype,FtempK);

%     Feed moisture content (w/w) 
      Fxh2o = findcomp(Fcomp,Fctype,comptype('water'));   

%     Composition at critical moisture (end of constant rate period) (w/w)
      FC = F*Fts/(1-FCcrit);
      FCctype = Fctype;
      for I = 1:length(Fcomp)
        if (Fcomp(I) ~= comptype('water'))  
          FCcomp(I) = F/FC .*Fcomp(I);
        else 
          FCcomp(I) = FCcrit;
        end  
      end    


%     Droplet temperature (K) = wet bulb temperature of inlet air
      Tdroplet = airwater(5,AItempK,AIpres,4,AIw);
 
%     Droplet final diameter (m)
%     diamexit = diam95*(Frho/Prho*(1+FCcrit)/(1+Fxh2o))^(1/3) 
      diamexit = diammean*(Frho*Fts/(Prho*Pts))^(1/3);  
 
%     Average droplet diameter (m)
      diamavg = (diammean+diamexit)/2;
 
%     Determination of air temperature at critical point (K)
%     Water evaporated up to critical point (kg/h)
      h2ocrit = F - FC;
%     Mass balance solved for absolute humidity of air at critical point
      Waircrit = h2ocrit/Ga + AIw;
%     Heat capacity of product at critical point (kJ/kg-K)
      FCcp = cp(FCcomp,FCctype,Tdroplet);
%     Total energy of product at critical point (kJ/h)
      Qcrit = FC*FCcp*(Tdroplet - Tref);
%     Energy balance solved for enthalpy of air at critical point (kJ/kg dry air)
      Haircrit = (Fhtot + AIhtot - Qcrit)/Ga;
%     Critical air temperature (K)
      Taircrit = (Haircrit - 2501*Waircrit)/(1.006 + 1.805*Waircrit) +Tref;
%     Log mean temperature difference across dryer (K) - cocurrent
      delTlmc = ((AItempK-FtempK)-(Taircrit-Tdroplet))/log((AItempK-FtempK)/(Taircrit-Tdroplet));
%     Thermal conductivity of droplet (kJ/h-m-K)
      kfeed = thermc(Fcomp,Fctype,FtempK);
      kcrit = thermc(FCcomp,FCctype,Tdroplet);
      kavg = (kfeed+kcrit)/2;
 
%     Latent heat of vaporization at product temperature (kJ/kg)
%     latx1 = hsteam(PtempK,0,1);
      latx1 = hsteam(Tdroplet,0,1);
%     latx0 = hsteam(PtempK,0,0);
      latx0 = hsteam(Tdroplet,0,0);
      latheat = latx1-latx0;
%     Constant rate (kg/h)
      Rconst = 2*pi*kavg*diamavg*delTlmc/latheat;
 
%     Feed moisture content per droplet (kg water)
      Mo = (Frho*Fcomp(comptype('water'))*4*pi*(diammean/2)^3)/3;
%     Feed solids content per droplet (kg solids)
      Mos = (Frho*Fts*4*pi*(diammean/2)^3)/3;
%     Critical moisture content per droplet (kg water)
      Mcrit = Mos*FCcrit/(1-FCcrit);
%     Time to dry to critical moisture content (h-->s)
      tcrits = (Mo - Mcrit)/Rconst*3600;
 
%     Log mean temperature difference across dryer (K) - cocurrent
      delTlmf = ((Taircrit-Tdroplet)-(AOtempK-PtempK))/log((Taircrit-Tdroplet)/(AOtempK-PtempK));
 
%     Thermal conductivity of air (kJ/h-m-K)
      kairout = kair(AOtempK);
 
%     Density at critical moisture content (kg/m^3)
      rhocrit = (Mos+Mcrit)/(1/6*pi*diamexit^3);
%     Average density in falling rate period (kg/m^3)
      rhobulk = (rhocrit+Prho)/2;
 
%     Falling rate (kg water/kg total-h)
      Rfall = 12*kairout*delTlmf/(latheat*diamexit^2*rhobulk);
 
%     Evaporation rate in falling period (kg water/h)
      Erate = Rfall*Mos;
 
%     Exit moisture content per droplet (kg)
      Mexit = Mos*Pcomp(comptype('water'))/(1-Pcomp(comptype('water')));
%     Time to dry in falling rate period (h-->s)
      tfalls = (Mcrit-Mexit)/Erate*3600;
 
%     Total drying time (s)
      ttotals = tcrits + tfalls;
 
%     Determination of resultant velocity from rotary vaned wheel
%     Droplet viscosity (Pa-s --> cP) 
%     Power law values for viscosity at feed temperature
      powerlaw = viscos(Fcomp,Fctype,FtempK,mumodel);
      Kbulk = powerlaw(1);       % Consistency Coefficient (Pa-s^n)
      nbulk = powerlaw(2);       % Flow behavior index
 
%     Bulk viscosity (Pa-s); taken at max velocity in tube
      mudrop = Kbulk;
 
      mudrop = vissteam(FtempK,0,0)*1000;
%     Volumetric feed flow rate (m^3/h)
      Vfeed = F/Frho;
%     Radial velocity normal to wheel (m/s)
      Vr = 0.0024*(Frho*pi^2*whRPM^2*whdiam*Vfeed^2/ ...
            (mudrop*vheight^2*whvanes^2))^(1/3);
 
%     Tangential velocity tangent to wheel (m/s)
      Vt = pi*whdiam*whRPM/60;
 
%     Resultant velocity (m/s)
      Vres = (Vr^2 + Vt^2)^0.5;
 
%     Bulk viscosity of air (Pa-s)
      mairbulk = muair((AItempK+AOtempK)/2);
 
%     Bulk density of air (kg/m^3)
      AIrho = 1/airwater(3,AItempK,AIpres,4,AIw);
      AOrho = 1/airwater(3,AOtempK,AIpres,4,AOw);
      rairbulk = (AIrho+AOrho)/2;
 
%     Time increment (s)
      dt = 0.0005;
 
%     Total number of iterations plus 1 for time zero
      IT = ceil(ttotals/dt)+1;
 
%     Initial (t = 0) conditions
      time(1) = 0;
      Uh(1) = Vres;  % Horizontal and vertical velocities (m/s)
      Uv(1) = 0;
      U(1) = Uh;     % Resultant velocity
      Sh(1) = 0;     % Horizontal and vertical distance (m)
      Sv(1) = 0;
 
%     Reynold number for droplet in air
      NRedrop = rairbulk*diamavg*U(1)/mairbulk;
 
%     Drag coefficient on droplet (Perry's)
      CD = dragcoef(NRedrop);
 
%     Mass of droplet (kg)
      massdrop = Mo + Mos;
 
%     Perform 2-dimensional distance and velocity analysis until horizontal
%      distance remains nearly constant.
      for t = 2:IT
 
%       Terminal time increment conditions based on conditions at the
%        commensement of the time interval
 
%       Horizontal velocity differential (m/s)
        dUhc(t) = -3*rairbulk*CD*U(t-1)*Uh(t-1)*dt/(2*rhobulk*diamavg);
 
%       Vertical velocity differential (m/s)
        dUvc(t) = 9.81*(Frho-rairbulk)/Frho*dt-3*rairbulk*CD*U(t-1)*Uv(t-1)*dt/(2*diamavg*Frho);
 
%       Horizontal velocity at end of time incremement (m/s)
        Uh(t) = Uh(t-1) + dUhc(t);
 
%       Vertical velocity at end of time incremement (m/s)
        Uv(t) = Uv(t-1) + dUvc(t);
 
%       Resultant velocity (m/s)
        U(t) = (Uh(t)^2 + Uv(t)^2)^0.5;
 
%       Terminal time increment conditions based on conditions at the
%        end of the time interval
 
%       Reynold number for droplet in air
        NRedrop = rairbulk*diamavg*U(t)/mairbulk;
 
%       Drag coefficient on droplet (Perry's)
        CD = dragcoef(NRedrop);
 
%       Horizontal velocity differential (m/s)
        dUht(t) = -3*rairbulk*CD*U(t)*Uh(t)*dt/(2*rhobulk*diamavg);
 
%       Vertical velocity differential (m/s)
        dUvt(t) = 9.81*(Frho-rairbulk)/Frho*dt-3*rairbulk*CD*U(t)*Uv(t)*dt/(2*diamavg*Frho);
 
%       Average velocity changes based on commencement and end conditions
        dUhavg(t) = (dUhc(t)+dUht(t))/2;
        dUvavg(t) = (dUvc(t)+dUvt(t))/2;
 
%       Initialize new average velocites (m/s)
        dUhavgn(t) = dUhavg(t) + 1;
        dUvavgn(t) = dUvavg(t) + 1;
        iter = 0;
 
        while (abs(dUhavgn(t)-dUhavg(t))>1e-5) & (abs(dUvavgn(t)-dUvavg(t))>1e-5)
 
          if (iter > 0)           % Reassign average velocities to new values
            dUhavg(t) = dUhavgn(t);
            dUvavg(t) = dUvavgn(t);
          end
 
%         Horizontal velocity at end of time incremement (m/s)
          Uh(t) = Uh(t-1) + dUhavg(t);
 
%         Vertical velocity at end of time incremement (m/s)
          Uv(t) = Uv(t-1) + dUvavg(t);
 
%         Resultant velocity (m/s)
          U(t) = (Uh(t)^2 + Uv(t)^2)^0.5;
 
%         New Terminal time increment conditions based on average conditions
%          end of the time interval
 
%         Raynold number for droplet in air
          NRedrop = rairbulk*diamavg*U(t)/mairbulk;
 
%         Drag coefficient on droplet (Perry's)
          CD = dragcoef(NRedrop);
 
%         Horizontal velocity differential (m/s)
          dUhnew(t) = -3*rairbulk*CD*U(t)*Uh(t)*dt/(2*rhobulk*diamavg);
 
%         Vertical velocity differential (m/s)
          dUvnew(t) = 9.81*(Frho-rairbulk)/Frho*dt-3*rairbulk*CD*U(t)*Uv(t)*dt/(2*diamavg*Frho);
 
%         Average velocity changes based on commencement and end conditions
          dUhavgn(t) = (dUhnew(t)+dUhavg(t))/2;
          dUvavgn(t) = (dUvnew(t)+dUvavg(t))/2;
 
          iter = iter + 1;
        end
 
%       Determination of terminal mass of droplet (kg) - MAY NEED FOR
%         PREDICTING CHANGES IN DENSITY
%       Current time (s)
        time(t) = t*dt;
        if (time <= tcrits)
%         Moisture loss over time interval
          dM = Rconst*dt/3600;
        else
          dM = Erate*dt/3600;
        end
%       New water content (kg)
%       Mo = Mo - dM;
%       New mass of droplet (kg)
%       massdrop = Mo*Mos;
 
%       Distance (m)
        Sh(t) = Sh(t-1) + (Uh(t-1)+Uh(t))*dt/2;
        Sv(t) = Sv(t-1) + (Uv(t-1)+Uv(t))*dt/2;
 
%       Terminal velocity
%       Uterm = (4*(Frho-rairbulk)*9.81*diamavg/(3*CD*rairbulk))^0.5;
%       Uterm = 0.225*diamavg^1.14*9.81^0.71*(Frho-rairbulk)/(rairbulk^0.29*mairbulk^0.4);
        NReh = rairbulk*diamavg*Uh(t)/mairbulk;  
        if (NReh <= 0.2)
          subplot(2,2,1), plot(time,Uh);
          xlabel('Time (s)');
          ylabel('Horizontal Velocity (m/s)');
          title ('Spray Dryer');
          subplot(2,2,2), plot(time,Sh);
          xlabel('Time (s)');
          ylabel('Horizontal Distance (m)');
          title ('Spray Dryer');
          subplot(2,2,3), plot(time,Uv);
          xlabel('Time (s)');
          ylabel('Vertical Velocity (m/s)');
          title ('Spray Dryer');
          subplot(2,2,4), plot(time,Sv);
          xlabel('Time (s)');
          ylabel('Veritcal Distance (m)');
          title ('Spray Dryer');
          break;
        end
      end
 
%     Velocities and distances for remaining time
%      Time in laminar flow (s); assumed to be all remaining drying time
      trem = ttotals - time(t);
      t = t+1;
      time(t) = ttotals;
%     Constants
      a = 9.81*(1-rairbulk/Frho);
      b = 18*mairbulk/(diamavg^2*Frho);
%     Horizontal velocity (m/s)
      Uh(t) = Uh(t-1)*exp(-b*trem);
%     Horizontal distance (m)
      Sh(t) = Sh(t-1) + Uh(t)/a*(1-exp(-a*trem));
%     Vertical velocity (m/s)
      Uv(t) = b/a + (Uv(t-1) - b/a)*exp(-a*trem);
%     Vertical distance (m)
      Sv(t) = Sv(t-1) + b*trem/a+Uv(t-1)/a-b/a^2+(b/a^2-Uv(t-1)/a)*exp(-a*trem);
%     Total chamber height (m)
      htotal = Sv(t);
%     Chamber diameter & radius  (m)
      dcyl = htotal/htod;
      rcyl = dcyl/2;
%     Cone height (m)
      hcone = rcyl*sin(phi);
%     Cyclindrical height (m)
      hcyl = htotal - hcone;
 
%     Volume of cone (m^3)
      Vcone = pi/3*rcyl^3*sin(phi);
%     Volume of cyclinder (m^3)
      Vcyl = pi*rcyl^2*hcyl;
%     Total volume (m^3)
      Vtotal = Vcone + Vcyl;
 
%     Minimum Residence Time (s)
      tmin = Vtotal*rairbulk/((AI+AO)/2)*3600;
 
%     Surface area of cone (m^2)
      Acone = pi*rcyl*(rcyl^2 + hcone^2)^0.5;
%     Surface area of cyclinder (m^2)
      Acyl = 2*pi*rcyl*hcyl + pi*rcyl^2;
%     Total surface area (m^2)
      Atotal = Acone + Acyl;

      endtxt = [10 'END OF SPRAY DRYER (',Uname,') ADVANCED TRANSPORT' 10]; 
      disp(endtxt); 

      end  % Advanced Transport Level 

%  *************** MICROBIAL & QUALITY VARIABLE ASSIGNMENTS *****************

%     Feed Microbial Death Kinetics
      Dvalm = Fuarray(ulocate('microDval'));   % D value (min)
      DtempCm = Fuarray(ulocate('microDtemp'));% D value temperature basis (C)
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
      Qoin = Fuarray(ulocate('qualload'));     % quality load
      Fqual = [Dvalq;DtempCq;Zvalq;Eaq;Qoin];  % quality array 
 
%  ********************* MICROBIAL & QUALITY ASSESSMENT *********************

%     Pass microbial and quality assessment variables to exit streams
      Pmdk = Fmdk;  
      Pqual = Fqual;  

      if (mandq == 1) & (dlevel >= 3)   

%     Bulk temperature for microbial assessment (K)
      Tbulk = (FtempK+PtempK)/2;  

%     Minimum Residence Time (min)
      tminm = tmin/60;  
 
%     Final microbial load 
      Nout = mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,Tbulk,tminm);  
 
%     Log reduction in microbial load
      DlogNred = log10(Noin/Nout);

%     Temperature for quality assessment (K)
      Tqual = PtempK;

%     Maximum processing time = Minimum processing time (min) 
      tqmax = ttotals/60;

%     Final quality load  
      Qout =  qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tqual,tqmax); 

%     Log reduction in quality 
      DlogQred = log10(Qoin/Qout);

%     Product stream microbial and quality arrays
      Pmdk = Fmdk;  
      Pmdk(size(Fmdk,1)) = Nout;
      Pqual = Fqual;
      Pqual(size(Fqual,1)) = Qout;

      endtxt = [10 'END OF SPRAY DRYER (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
      disp(endtxt); 

%     Update system output array
      Puarray = [Pme;Pmdk;Pqual];

      end  % microbial and quality assessment 

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(21:30); % economic analysis array 
      end  

%  ********************* SYSTEM OUTPUT ARRAY ASSIGNMENT *********************

%     Size stream info arrays to full system output u-arrays  
      Puarray = ustd(Puarray);
      AOuarray = ustd(AOuarray);
      Cuarray = ustd(Cuarray);

      sys =[Puarray;AOuarray;Cuarray];    
 
%  ****************************** PRINTED OUTPUT ****************************

  if flag == 9    % print only after a solution is found 

%   Print results to specified file name; 
%    call filename manager function rfile.m 
    fid = rfile;

%   Name of unit operation; system input streams; system output streams  
    UNAME = Uname; 
    SINAMES = str2mat(Fname,AAname,Sname);
    if (dlevel == 1) 
      Hname = 'WaterRemoved'  
      SONAMES = str2mat(Pname,Hname);
    else 
      SONAMES = str2mat(Pname,AOname,Cname);
    end 

    if (dlevel == 1)
%     System inputs array 
      SIARRAY = [Fmass]; 
%     System outputs array 
      SOARRAY = [Pmass Hmass];
    elseif (dlevel >= 2) 
%     System inputs array
      SIARRAY = [Fmass AAmass Smass];
%     System outputs array
      SOARRAY = [Pmass AOmass Cmass];
    end  

%   No intermediate streams 
    SNNAMES = [];
    SNARRAY = []; 
    if (dlevel >= 2) 
      ARname = 'AirRecycle'; 
      AIHXname = 'AirInHX';
      AIname = 'AirOutHX';
      AOSDname = 'AirOutSD';
      if (frecycle <= 0) 
        SNNAMES = str2mat(AIname);
        SNARRAY = [AImass];
      else 
        SNNAMES = str2mat(AIHXname,AIname,AOSDname,ARname);
        SNARRAY = [AIHXmass AImass AOSDmass ARmass];
      end  
    end   

%   Print mass balance results table
    msprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    if (dlevel >= 2) 
%     System inputs array  
      SIARRAY = [Fme AAme Sme]; 
%     System outputs array 
      SOARRAY = [Pme AOme Cme];
%     System intermediate stream array 
      if (frecycle <= 0) 
        SNARRAY = [AIme];
      else 
        SNARRAY = [AIHXme,AIme,AOSDme,ARme];
      end  

%     Print energy balance results table
      eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    end  

    fprintf(fid,'UNIT MODEL INPUTS\n');     
    fprintf(fid,'     \n');
    fprintf(fid,'  Dry Product Total Solids                %6.5f w/w \n',Pts); 

    if (dlevel >= 2) 
      fprintf(fid,'  Ambient Air Relative Humidity             %5.2f percent \n',AARH); 
      fprintf(fid,'  Dry Product Temperature                   %5.2f C \n',PtempC); 
      fprintf(fid,'  Dryer Air Inlet Temperature              %5.2f C \n',AItempC); 
      fprintf(fid,'  Dryer Air Outlet Temperature              %5.2f C \n',AOtempC); 
      fprintf(fid,'  Operating Pressure                    %5.3e MPa \n',AIpres); 
      fprintf(fid,'  Dry Air Recycle Fraction                  %5.3f \n',frecycle); 
    end

    if (dlevel >= 4) 
      fprintf(fid,'  Critical Moisture Content               %6.5f w/w\n',FCcrit);
      fprintf(fid,'  Product Density                        %8.3f kg/m^3\n',Prho);
      fprintf(fid,'  Viscosity Model                        %8.0f \n',mumodel);

      fprintf(fid,'     \n');
      fprintf(fid,' Atomizer Design (vaned wheel):\n');
      fprintf(fid,'  Wheel Speed                            %8.0f RPM\n',whRPM);
      fprintf(fid,'  Wheel Diameter                         %8.4f m\n',whdiam);
      fprintf(fid,'  Number of Vanes                        %8.0f \n',whvanes);
      fprintf(fid,'  Vane Height                            %8.4f m\n',vheight);
      fprintf(fid,'  Mean Droplet Diameter                 %4.3e m\n',diammean);

      fprintf(fid,'     \n');
      fprintf(fid,' Chamber Dimensions:\n');
      fprintf(fid,'  Height:Diameter                        %8.2f \n',htod);
      fprintf(fid,'  Cone Angle                             %8.2f deg.\n',phideg);
    end  



    if (dlevel >=2) 
    fprintf(fid,'     \n');
    fprintf(fid,'UNIT MODEL OUTPUTS\n');     
    fprintf(fid,'     \n');
      fprintf(fid,'  Dry Air Mass Flow Rate               %10.2f kg air/h \n',Ga); 
      fprintf(fid,'  Dryer Air Inlet Relative Humidity         %5.2f percent \n',AIRH); 
      fprintf(fid,'  Dryer Air Outlet Relative Humidity        %5.2f percent \n',AORH); 
    end 

    if (dlevel >=4) 
    fprintf(fid,'     \n');
    fprintf(fid,' Atomizer Design (vaned wheel):\n');
    fprintf(fid,'  Droplet Diameter (95 percentile)      %4.3e m\n',diam95);
    fprintf(fid,'  Consistency Coefficient               %8.3e Pa-s^n\n',Kbulk);
    fprintf(fid,'  Flow Behavior Index (n)                %8.3f \n',nbulk);
    fprintf(fid,'  Droplet Viscosity                     %8.3e Pa-s \n',mudrop);
    fprintf(fid,'  Resultant Velocity                     %8.4f m/s\n',Vres);
    fprintf(fid,'  Power Consumption                      %8.2f kW\n',power);
    fprintf(fid,'     \n');
 
    fprintf(fid,' Constant Rate Period:\n');
    fprintf(fid,'  Feed Moisture Content (per drop)      %8.3e kg water\n',Mo);
    fprintf(fid,'  Critial Moisture Content (per drop)   %8.3e kg water\n',Mcrit);
    fprintf(fid,'  Average Droplet Diameter              %8.3e m\n',diamavg);
    fprintf(fid,'  Product Droplet Diameter              %8.3e m\n',diamexit);
    fprintf(fid,'  Droplet Temperature                    %8.3f C\n',Tdroplet-tref);
    fprintf(fid,'  Critical Air Temperature               %8.3f C\n',Taircrit-tref);
    fprintf(fid,'  Log Mean Temperature Difference        %8.3f C\n',delTlmc);
    fprintf(fid,'  Latent Heat of Vaporization            %8.3f kJ/kg\n',latheat);
    fprintf(fid,'  Feed Density                           %8.3f kg/m^3\n',Frho);
    fprintf(fid,'  Bulk Thermal Conductivity              %8.3f kJ/h-m-K\n',kavg);
    fprintf(fid,'  Constant Drying Rate                  %8.3e kg/h\n',Rconst);
    fprintf(fid,'  Time in Constant Rate Period           %8.3f s\n',tcrits);
    fprintf(fid,'     \n');
 
    fprintf(fid,' Falling Rate Period:\n');
    fprintf(fid,'  Product Moisture Content (per drop)   %8.3e kg water\n',Mexit);
    fprintf(fid,'  Log Mean Temperature Difference        %8.3f C\n',delTlmf);
    fprintf(fid,'  Exit Air Thermal Conductivity          %8.3f kJ/h-m-K\n',kairout);
    fprintf(fid,'  Bulk Density                           %8.3f kg/m^3\n',rhobulk);
    fprintf(fid,'  Average Rate of Evaporation            %8.3f kg water/kg product-h\n',Rfall);
    fprintf(fid,'  Evaporation Rate (per drop)           %4.3e kg/h\n',Erate);
    fprintf(fid,'  Time in Falling Rate Period            %8.3f s\n',tfalls);
    fprintf(fid,'     \n');
 
    fprintf(fid,' Total Drying Time                       %8.3f s\n',ttotals);
    fprintf(fid,'     \n');
 
    fprintf(fid,' Chamber Dimensions:\n');
    fprintf(fid,'  Total Chamber Height                   %8.2f m\n',htotal);
    fprintf(fid,'   Cylinder Height                       %8.2f m\n',hcyl);
    fprintf(fid,'   Cone Height                           %8.2f m\n',hcone);
    fprintf(fid,'  Chamber Diameter                       %8.2f m\n',dcyl);
    fprintf(fid,'  Maximum Horizontal Distance            %8.2f m\n',Sh(t));
    fprintf(fid,'  Chamber Volume                         %8.2f m^3\n',Vtotal);
    fprintf(fid,'  Chamber Surface Area                   %8.2f m^2\n',Atotal);
    fprintf(fid,'     \n');
   
    end  

    if (dlevel >=4) 
    fprintf(fid,'     \n');
    
    endtxt = [10 'END OF (',Uname,') TRANSPORT RESULTS PRINT TO FILE' 10];  
    disp(endtxt);  

    end  

%   Microbial and Quality Assessment results printout 
    if (mandq == 1) & (dlevel >= 3)

    fprintf(fid,'     \n');
    fprintf(fid,'  Microbial (Worst Case) Assessment:   \n');
    fprintf(fid,'   Minimum Residence Time                  %5.4f min \n',tminm);
    fprintf(fid,'   Minimum Processing Temperature           %5.2f C \n',Tbulk-Tref);
    fprintf(fid,'   Log Reduction in Microbial Load          %5.2e \n',DlogNred);
    fprintf(fid,'     \n'); 
    fprintf(fid,'  Quality (Worst Case) Assessment:   \n');
    fprintf(fid,'   Maximum Residence Time                  %5.4f min \n',tqmax);
    fprintf(fid,'   Maximum Processing Temperature           %5.2f C \n',Tqual-Tref);  
    fprintf(fid,'   Log Reduction in Quality Load            %5.2e \n',DlogQred);

%     System input and exit stream names for microbial and quality assessment 
      SIMQNAME = str2mat(Fname);
      SOMQNAME = str2mat(Pname); 

%     System input microbial array
      SIMARRAY = [Fmdk];
%     System outputs microbial array
      SOMARRAY = [Pmdk;tminm;Tbulk-Tref;];

%     No intermediate streams
      SNMQNAME = [];
      SNMARRAY = [];
      SNQARRAY = [];

      mdkprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIMARRAY,SOMARRAY,SNMARRAY);

%     System input quality array
      SIQARRAY = [Fqual];
%     System outputs quality array
      SOQARRAY = [Pqual;tqmax;Tqual-Tref;];

      qprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIQARRAY,SOQARRAY,SNQARRAY);
    end 

    if (econ >= 1) 

      Uutil = [S 0 0 0 0 0]; 
      Usize = H;
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
 
