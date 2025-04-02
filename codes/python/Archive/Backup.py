def create_parameters():
    return     {
        "General simulation parameters" : { 
                "Diameter":    [500," The diameter of the asteroid"],                 
                "Density":     [1250," The density of the asteroid"],
                "Angular velocity":     [0.7," The non-dimensional angular velocity of the asteroid"],
                "Simulation period":    [1e+5,"Total simulation time of the asteroid in main-belt "],
                "Collision":   [True," Do you want to include the effect of collision on spin rate?"],
                "Landslide" :  [True,"Do you want to include the effect of lanslide on spin rate and shape? "],
                "YORP":        [True," Do you want to include the effect of YORP on spin rate?"]     
        },
        
        "Landslide" : {
                    "Gamma":       [0.1,"It is the non-dimensional scale of the topography. Keep it smaller than 1."] , 
                    "Friction angle":       [15, "Friction angle of the regolith"],
                    "dump":        [200, """The timesteps after which the output of the landslide to be stored in the file. The location of output will be avialable in the parameter file under option verbose_dir. The option is valid only when verbose is selected"""],
                    "epsilon":     [0.002, """The non-dimensional height of the failed layer of the regolith. 
                                      It should be smaller than Gamma and <<1"""],
                    "executable":  ["asteroid","This is the executable file name which is run for the landslide simulation."],
                    "offset":      [0.01," To circumvent the pole problem., we neglect the small region near the pole. The value is in radians"], 
                    "Resolution":         [500," The grid size for the landslide simulations"],
                    "theta":       [1.0," Top be used in the min-mod limiter"],                    
                    "weight":      [0.5," To be used for the Runge-Kutta scheme"],
                    "verbose":     [True," If you want to see the detailed simulation of landslides "],                    
                    "Profile":     ["Uniform","The profile of the failed regolith layer",["Uniform", "Gaussian"]],
                    "Landslide simulation period":        [20," The simulation period after which the landslide will stop. This is the highest time for which landslide will be simulated"],
                    "Friction type":   ["Constant","If you want a constant friction angle through out the lanslide, select constant otherwise choose variable",["Constant","variable"]]
                        },
        "Collision": {
                    "Cummulative distribution":    ["KAH-SFD","Select a population density for the impactor",["KAH-SFD","Bottke2005","SFDtruncated"]]
                    },
                  
        "YORP": {
                    "K":           [0.01,"Thermal inertia "],
                    "atype":       ["C-Type","Typeof asteroid, C-Type or S-Type",["C-Type","S-Type"]],
                    "c_YORP":      [ 0.7," For calculating YORP torques"],                   
                    "c_reor":      [0.9,"To calulate the probability of collisional reoriention after at very slow spin rate "],
                    "Obliquity":       [135.0,"The obliquity of the asteroid "],
                    "Semi major axis":         [2.3," The semi major axis of the asteroid"],
                    "h_yorp":      [50," The time step for the YORP calculation in years"],            
                    "stoc_yorp":   [False, "Do you want yorp to be stochastic"]
                }           
    }
                                    
                                    