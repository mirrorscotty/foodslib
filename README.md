FOODSLIB - Food Operations Oriented Design System Block Library 
===============================================================
Version 1.0c  1997 

Steady-State Library and Demos S-Functions
------------------------------------------
    foodslib     - Block library, similar to SIMULINK  
    evapsdemo    - Single effect evaporator block diagram - steam rate unknown 
    evapsdemo2   - Single effect evaporator block diagram - steam rate known 
    evapffdemo   - Forward feed multiple effect block diagram  
    evaprfdemo   - Reverse feed multiple effect block diagram  
    econanal     - Economic analysis block diagram
    mix2demo     - Point mix for 2 streams block diagram   
    mix3demo     - Point mix for 3 streams block diagram   

System Input and Output S-Functions 
-----------------------------------
    feeds        - Feed variables u-array re-assignment
    steams       - Steam or water source variables u-array re-assignment
    products     - Product variables u-array re-assignment
    wastes       - Waste variables u-array re-assignment  

Steady-State Unit Operation Models S-Functions
----------------------------------------------
    evaps        - Single effect evaporator  
    evapffs      - Multiple-effect feed-forward evaporator  
    evaprfs      - Multiple-effect reverse-feed evaporator 
    fittings     - Pipe fittings  
    holdtbs      - Holding tube
    hxs          - Heat exchanger for hot/cold water/steam and food streams
    mix2s        - Point mix for 2 streams 
    mix3s        - Point mix for 3 streams 
    pipings      - Straight pipe   
    pumps        - Pump   
    spdrys       - Spray dryer 
    splitms      - Mass split 
    sxpcons      - Sudden expansion and contraction in pipe    

Stream Management
-----------------
    ustd         - Standardize length of stream u-array
    usizemax     - Current size of u-array for all streams
    ulocate      - Locate index in u-array for specified information
    ucompmax     - Current number of components compositions in stream u-array

Stream Property Functions
-------------------------
    cp           - Heat capacity of multi-component stream  
    cpair        - Heat capacity of air 
    enthalpy     - Enthalpy of any stream type  
    density      - Density of any stream type 
    thcond       - Thermal conductivity of any stream type
    hsteam       - Enthalpy of steam 
    ptsteam      - Saturation pressure given steam temperature 
    tpsteam      - Saturation temperature given steam pressure  
    unfroz       - Frozen & unfrozen water fractions of multi-component stream
    tsolids      - Total solids of stream based on composition
    rho          - Density of multi-component stream 
    rhoair       - Density of air 
    thermc       - Thermal conductivity of multi-component stream 
    vissteam     - Dynamic viscosity of steam 
    muair        - Viscosity of air 
    kair         - Thermal conductivity of air  
    ksteam       - Thermal conductivity of steam 
    rhosteam     - Density of steam 
    molwts       - Molecular weights of stream components
    bpr          - Boiling point rise of multi-component stream
    mixtemp      - Temperature of mixed streams
    viscos       - Viscosity models  

Dimensionless Numbers
---------------------
    reynolds     - Reynolds number 
    prandtl      - Prandtl number

Constants
---------
    gaslawc      - Gas law constant 
    gravity      - Standard acceleration of gravity  
    tref         - Reference temperature 

Mass, Heat, and Momentum Transfer Functions
-------------------------------------------
    dragcoef     - Drag coefficient on a sphere  
    hostcond     - Steam side heat transfer coefficient for film-condensation 
    hievap       - Product side heat transfer coefficient models  
    shrwall      - Shear rate at the pipe wall 

Microbial and Quality Functions
-------------------------------
    mdeathk      - Thermal death-rate kinetics of microorganisms
    qualdegr     - Thermal degradation of quality

File Management and Printing Functions
--------------------------------------
    rfile        - Creates user-specified results file for design simulation  
    scfile       - Creates user-specified stream cost/value data file
    uofile       - Creates user-specified unit operation cost data file 
    eprint       - Prints energy balance results to results file  
    mprint       - Prints mass balance results to results file  
    mdkprint     - Thermal death-rate kinetics of microbes results printout
    qprint       - Quality kinetics results printout

Economic Analysis Functions
---------------------------
    econs        - Economic analysis model  
    inflate      - Inflates year 1 value to year N value  

Miscellaneous Functions
-----------------------
     findval1    - Find index of all elements in 1D array with value specified
     findcomp    - Find mass fraction of specified component in component array
     namesize    - Standard name length for printing tables    
     comptype    - Index for component type

