% **********************************************************************
%   Program to Model Two Dimensional Heat Flow in Can in Retort 9/8/92 
% **********************************************************************
%   Modified October 17, 1992 (Karen Grenus)  
%   Converted to Matlab Language August 29, 1996 (Heidi A. Diefes)  

    function [out] = can(u)    

%   Reference:  Geankoplis, C.J. 1983.  "Transport Processes and Unit
%     Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.

%   Assignment of variables supplied to function
%   Steam supply information  

    Psteam = u(1);       % pressure (MPa) - NOT USED IN ALGORITHM  
    xsteam = u(2);       % quality - NOT USED IN ALGORITHM  
    retT = u(3);         % retort temperature (K)

%   Food in can information
    initT = u(4);        % initial food material temperature (K)
    xpro = u(5);         % protein content (w/w)
    xcarb = u(6);        % carbohydrate content (w/w)
    xfat = u(7);         % fat content (w/w)
    xash = u(8);         % ash component (w/w)
    xh20 = 1 - (xpro+xcarb+xfat+xash);  % water component (w/w)  
%   Food microbial death kinetics 
    D = u(9);            % D value (min) for selected microbe    
    Dtemp = u(10);       % Temperature (C) basis for D value  
    Zo = u(11);          % Z value (C) for selected microbe 
%   Food quality kinetics
    Dq = u(12);          % D value (min) for selected microbe    
    Dqtemp = u(13);      % Temperature (C) basis for D value  
    Zq = u(14);          % Z value (C) for selected microbe  

%   Can and finite-difference information
    height = u(15);      % Can height (cm)  
    diam = u(16);        % Can diameter (cm)  
    Ztotal = u(17);      % Number of axial nodes; must be odd 
    Rtotal = u(18);      % Number of radial nodes; must be odd 
    heattime = u(19);    % Total heating time (min)     
    Dlogred = u(20);     % Desired log reduction    
    h = u(21);           % Heat transfer coefficent to can (W/cm^2 K)  

%   Increments in radial and vertical directions (cm)
    delR = diam/2/Rtotal;   
    delZ = height/Ztotal;  

%   Initialize Fo and Fq  
    for R = 1:Rtotal
      for Z = 1:Ztotal  
        Fo(R,Z) = 0.0;  
        Fq(R,Z) = 0.0;  
      end
    end  

%   Thermal properties of food material  
%   Density of the product (gm/cm^3): function units (kg/m^3)     
    den = rho(xpro,xcarb,xfat,xash,xh20,(initT+273.15))/1000;  
%   Heat capacity of the product (J/g-K)
    cpfood = cp(xpro,xcarb,xfat,xash,xh20,(initT+273.15));     
%   Thermal conductivity fo the product (W/cm-K): function units (kJ/h-m-K)  
    cond = thermc(xpro,xcarb,xfat,xash,xh20,(initT+273.15))/360;   
%   Thermal diffusivity, cm^2/s
    alpha=cond/(cpfood*den);   

%   Changing the value of h will radically change the number of times 
%    that the temperature is calculated it will also impede or enhance 
%    how well the heat is transferred to the product.  
%    If h is low, the step size is big and the time to heat up is low. 
%    If h is high, then the step size is small and the time to heat up 
%    is high.    

%   Calculating the time increment.  Depends which direction
%    needs the smallest time step in order to maintain stability.  

%   Calculate the time step in the r direction
	
    Rn = h*delR/cond;  
    Rm = 2*Rn+2;  
    Rdelt = delR^2/(alpha*Rm);  

%   Calculate the time step in the z direction

    Zn = h*delZ/cond;      
    Zm = 2*Zn+2;  
    Zdelt = delZ^2/(alpha*Zm);  

%   Selecting smaller time step (sec)
    if (Zdelt <= Rdelt) 
      delt = Zdelt;  
    else 
      delt = Rdelt; 
    end
%   Time step (min)  
    deltmin = delt/60;    
%   Total number of iterations  
    ntotal = round(heattime/deltmin);         

%   Print to results to specified file name  
    filem = input('Enter filename for results:   ');   
    fid = fopen(filem,'w+');   
    fprintf(fid,'     \n');  
    fprintf(fid,'CAN STERILIZATION RESULTS');
    fprintf(fid,'     \n');  
    fprintf(fid,'INPUT SUMMARY\n'); 
    fprintf(fid,'     \n');  
    fprintf(fid,' Retort Conditions: \n');  
    fprintf(fid,'  Steam Temperature       = %5.1f C \n',retT);
    fprintf(fid,'     \n');
    fprintf(fid,' Food Material Description: \n');  
    fprintf(fid,'  Initial Temperature     = %5.1f C \n',initT);
    fprintf(fid,'  Protein Content         = %5.3f w/w \n',xpro);
    fprintf(fid,'  Carbohydrate Content    = %5.3f w/w \n',xcarb);
    fprintf(fid,'  Fat Content             = %5.3f w/w \n',xfat);
    fprintf(fid,'  Ash Content             = %5.3f w/w \n',xash);
    fprintf(fid,'  Makeup Water Content    = %5.3f w/w \n',xh20);
    fprintf(fid,'     \n');
    fprintf(fid,' Microbial Death Kinetics: \n');  
    fprintf(fid,'  D Value                 = %5.1f min \n',D);
    fprintf(fid,'  D Temperature Basis     = %5.1f C \n',Dtemp);
    fprintf(fid,'  Z Value                 = %5.1f C \n',Zo);
    fprintf(fid,'     \n');
    fprintf(fid,' Quality Kinetics: \n');  
    fprintf(fid,'     \n');
    fprintf(fid,'  D Value                 = %5.1f min \n',Dq);
    fprintf(fid,'  D Temperature Basis     = %5.1f C \n',Dqtemp);
    fprintf(fid,'  Z Value                 = %5.1f C \n',Zq);

    fprintf(fid,'     \n');
    fprintf(fid,' Can Dimension and Model Information: \n');  
    fprintf(fid,'  Height                  = %5.1f cm \n',height);
    fprintf(fid,'  No. Axial Nodes         = %5.1f \n',Ztotal);
    fprintf(fid,'  Diameter                = %5.1f cm \n',diam);
    fprintf(fid,'  No. Radial Nodes        = %5.1f  \n',Rtotal);
    fprintf(fid,'  Total Heating Time      = %5.1f min \n',heattime);
    fprintf(fid,'  No. Log Reductions in D = %5.1f \n',Dlogred);
    fprintf(fid,'  Heat Transfer Coeff.    = %5.3f \n',h);
    fprintf(fid,'     \n');
    
    
    fprintf(fid,'     \n');  
    fprintf(fid,'MODEL EXECUTION SUMMARY\n');  
    fprintf(fid,'     \n');  
    fprintf(fid,' Time Step                = %8.5f (min)\n',deltmin);   
    fprintf(fid,' Axial Increment          = %8.3f (cm)\n',delZ);   
    fprintf(fid,' Radial Increment         = %8.3f (cm)\n',delR);   
    fprintf(fid,'     \n');  

%   Initial food material temperature (C)   
    for R = 1:Rtotal  
      for Z = 1:Ztotal    
        N = R*Ztotal - (Ztotal - Z);  % node number   
	T(N,1) = initT;               % temperature at node   
      end  
    end 
 
%   Boundary temperature (C) = Retort temperature  
    for IT = 1:ntotal  
      for R = 1:Rtotal  
        Z = 1;   
        N = R*Ztotal - (Ztotal - Z);  % nodes at bottom of can   
        T(N,IT) = retT;   
        Z = Ztotal;    
        N = R*Ztotal - (Ztotal - Z);  % nodes at top of can  
        T(N,IT) = retT;  
      end    
      for Z = 1:Ztotal  
        R = Rtotal;   
        N = R*Ztotal - (Ztotal - Z);  % nodes at side of can   
        T(N,IT) = retT;    
      end  
    end  

%   Initially not meeting sterilization criteria  
    steril = 0.0;  

%   Print header for time table  
    fprintf(fid,'Table 1. Conditions at center point of can over heating time.\n');  
    fprintf(fid,'     \n');  
    fprintf(fid,'  Time      Temp       Fo          Fq         Percent       Percent\n'); 
    fprintf(fid,'  (min)     (C)       (min)       (min)      Microbial      Quality\n');   
    fprintf(fid,'                                             Lethality    Degradation\n');   
    fprintf(fid,'-----------------------------------------------------------------------\n');   
    fprintf(fid,'%8.2f   %5.1f    %5.3e    %5.3e    %7.3f       %7.3f  \n',0,initT,0,0,0,0);

%   Temperature profiles over time (marked by iterations)
    for IT = 1:ntotal
      for IR = (Rtotal-1):-1:1
        for IZ = 2:(Ztotal-1)
          N = IR*Ztotal - (Ztotal - IZ);  
          NRP1 = (IR+1)*Ztotal - (Ztotal - IZ);  
          NRP2 = (IR+2)*Ztotal - (Ztotal - IZ); 
          NRM1 = (IR-1)*Ztotal - (Ztotal - IZ);  
          NZM1 = IR*Ztotal - (Ztotal - (IZ-1));  
          NZP1 = IR*Ztotal - (Ztotal - (IZ+1));  
          if (IR == 1)
%           Heat Diffusion Equation
            T(N,IT+1)=delt*alpha*(1/delR^2*(T(NRP1,IT)-2*T(N,IT)+T(NRP1,IT))+1/(delR^2*IR)*(T(NRP1,IT)-T(N,IT))+1/delZ^2*(T(NZP1,IT)-2*T(N,IT)+T(NZM1,IT)))+T(N,IT); 
          else
            T(N,IT+1)=delt*alpha*(1/delR^2*(T(NRP1,IT)-2*T(N,IT)+T(NRM1,IT))+1/(delR^2*IR)*(T(NRP1,IT)-T(N,IT))+1/delZ^2*(T(NZP1,IT)-2*T(N,IT)+T(NZM1,IT)))+T(N,IT);  
          end
        end
      end

%     Node at center of can
      Zcenter = round(Ztotal/2);  
      Ncenter = Zcenter; 
%     Node at axial center and half radius   
      Rcenter = round(Rtotal/2); 
      NRcenter = Rcenter*Ztotal - (Ztotal - Zcenter);  
%     Node at axial center and one radial node from edge
      Ecenter = (Rtotal-1)*Ztotal - (Ztotal - Zcenter);
 
%     Comparing; the kinetics for the center of the can to the
%     overall death for the microbe

%     Time (min) at can center  
      timeit = IT*delt/60;    
%     Temperature (C), Fo (min), Fq (min) at can center  
      Tcenter = T(Ncenter,IT);  
      Fo(1,Zcenter) = deltmin*10^((T(Ncenter,IT)-Dtemp)/Zo)+Fo(1,Zcenter);
      Fq(1,Zcenter) = deltmin*10^((T(Ncenter,IT)-Dqtemp)/Zq)+Fq(1,Zcenter);
%     Temperature (C), Fo (min), Fq (min) at can axial and radial center
      TRZcentr = T(NRcenter,IT);
      Fo(Rcenter,Zcenter) = deltmin*10^((T(NRcenter,IT)-Dtemp)/Zo)+Fo(Rcenter,Zcenter);
      Fq(Rcenter,Zcenter) = deltmin*10^((T(NRcenter,IT)-Dqtemp)/Zq)+Fq(Rcenter,Zcenter);
%     Temperature (C), Fo (min), Fq (min) at can outer edge   
      TEcenter = T(Ecenter,IT);
      Fo((Rtotal-1),Zcenter) = deltmin*10^((T(Ecenter,IT)-Dtemp)/Zo)+Fo((Rtotal-1),Zcenter);
      Fq((Rtotal-1),Zcenter) = deltmin*10^((T(Ecenter,IT)-Dqtemp)/Zq)+Fq((Rtotal-1),Zcenter);
 
%     Output variables for time plots    
      timeplot(IT) = timeit;    
      Tempplot(IT) = Tcenter;   
      TempplRZ(IT) = TRZcentr;   
      TempplE(IT) = TEcenter;   
      Foplot(IT) = Fo(1,Zcenter); 
      Fqplot(IT) = Fq(1,Zcenter);
      FoplotRZ(IT) = Fo(Rcenter,Zcenter);
      FqplotRZ(IT) = Fq(Rcenter,Zcenter);
      FoplotE(IT) = Fo((Rtotal-1),Zcenter);
      FqplotE(IT) = Fq((Rtotal-1),Zcenter);

%     Percent of Microbes Destroyed at can center   
      PML(1,Zcenter)=(10^(Fo(1,Zcenter)/D)-1)/10^(Fo(1,Zcenter)/D)*100;
%     Percent of Quality Destroyed at can center    
      PQD(1,Zcenter)=(10^(Fq(1,Zcenter)/Dq)-1)/10^(Fq(1,Zcenter)/Dq)*100;

%     Print contents of table for every 30 iterations  
      pcrit = rem(IT,30); 
      if (pcrit == 0) 
        IT  
        Tcenter
        timeit 
        fprintf(fid,'%8.2f   %5.1f    %5.3e    %5.3e    %7.3f       %7.3f  \n',timeit,Tcenter,Fo(1,Zcenter),Fq(1,Zcenter),PML(1,Zcenter),PQD(1,Zcenter));  
      end 
  
      if (Fo(1,Zcenter) >= D*Dlogred)
        istep = IT  
        steril = 1; 
        break;  
      end
    end
%   Print last iteration   
    fprintf(fid,'%8.2f   %5.1f    %5.3e    %5.3e    %7.3f       %7.3f  \n',timeit,Tcenter,Fo(1,Zcenter),Fq(1,Zcenter),PML(1,Zcenter),PQD(1,Zcenter));
    fprintf(fid,'-----------------------------------------------------------------------\n');   
    fprintf(fid,'     \n');  

    if (steril == 0)
      disp('PRODUCT DID NOT REACH DESIRED STERILITY!!! ');   
      fprintf(fid,'         PRODUCT DID NOT REACH DESIRED STERILITY!!! \n');   
      disp('NEED TO INCREASE TOTAL HEATING TIME!!! ');   
      fprintf(fid,'           NEED TO INCREASE TOTAL HEATING TIME!!! \n');   
      istep = ntotal;  
    end

    fprintf(fid,'     \n');  

 
%   Reset the lethality at each node in the can to zero
    for R = 1:Rtotal
      for Z = 1:Ztotal
        Fo(R,Z) = 0.0;   
        Fq(R,Z) = 0.0;   
      end
    end
 
     for IT = 2:istep   
      for R = 1:Rtotal
        for Z = 1:Ztotal
%         Time (min) at D associated temperature that will give same 
%          degree of sterility as actual processing temperature  
          N = R*Ztotal - (Ztotal - Z);  
          Fo(R,Z)=deltmin*10^((T(N,IT)-Dtemp)/Zo)+Fo(R,Z);
          Fq(R,Z)=deltmin*10^((T(N,IT)-Dqtemp)/Zq)+Fq(R,Z); 
%         Percent of Desired Lethality   
          Fop(R,Z)=(D-Fo(R,Z))/D*100;  
          Foq(R,Z)=(Dq-Fq(R,Z))/Dq*100;  
        end
      end
 
 
%    Print final lethality, quality, and temperature profile of can   

      if (IT == istep)
        timeit = IT*delt/60; 
        fprintf(fid,'     \n');
        fprintf(fid,'Table 2. Final can conditions at time t = %5.2f (min).\n',timeit);                    
        fprintf(fid,'     \n');
        fprintf(fid,'Axial         Radial Increment    \n');
        fprintf(fid,'Increments    Center ---> Edge    \n');
        fprintf(fid,' Top     \n');
        fprintf(fid,'  ^   \n');
        fprintf(fid,'  |   \n');
        fprintf(fid,' Bottom    \n');
        fprintf(fid,'     \n');
%       Print Final Temperature Profile  
        fprintf(fid,'     \n');
        fprintf(fid,' Temperature profile (C).\n');
        fprintf(fid,'     \n');
        for Z = 1:Ztotal
          for R = 1:Rtotal
            N = R*Ztotal - (Ztotal - Z);  
            if (R < Rtotal)
              fprintf(fid,'  %7.3f ',T(N,IT));
            else
              fprintf(fid,'  %7.3f \n',T(N,IT));
            end
          end
        end
%       Take transpose of profiles to orient nodes for upright can          
        Fot=Fo';  
        Fqt=Fq'; 
%       Print Final Fo Profile   
        fprintf(fid,'     \n');  
        fprintf(fid,'     \n');  
        fprintf(fid,' Fo profile (min).\n'); 
        fprintf(fid,'     \n');  
        for Z = 1:Ztotal 
          for R = 1:Rtotal
            if (R < Rtotal) 
              fprintf(fid,'  %7.3f ',Fot(Z,R));
            else    
              fprintf(fid,'  %7.3f \n',Fot(Z,R));  
            end 
          end  
        end  
%       Print Final Fq Profile   
        fprintf(fid,'     \n');  
        fprintf(fid,'     \n');  
        fprintf(fid,' Fq profile (min).\n');
        fprintf(fid,'     \n');  
        for Z = 1:Ztotal 
          for R = 1:Rtotal
            if (R < Rtotal)  
              fprintf(fid,'  %7.3f ',Fqt(Z,R));
            else    
              fprintf(fid,'  %7.3f \n',Fqt(Z,R));  
            end  
          end  
        end   
      end
 
    end

%   Plot results  
    subplot(2,2,1), plot(timeplot,Tempplot,timeplot,TempplRZ,timeplot,TempplE);
    xlabel('Time (min)');
    ylabel('Temperature (C)');
    title ('Temp vs Time Along Can Axial Center');
    subplot(2,2,2), plot(timeplot,Foplot,timeplot,FoplotRZ,timeplot,FoplotE); 
    xlabel('Time (min)');
    ylabel ('Fo (min)');
    title ('Fo vs Time Along Can Axial Center');
    subplot(2,2,3), plot(timeplot,Fqplot,timeplot,FqplotRZ,timeplot,FqplotE);
    xlabel('Time (min)');
    ylabel ('Fq (min)');
    title ('Fq vs Time Along Can Axial Center');

    out = [];  
