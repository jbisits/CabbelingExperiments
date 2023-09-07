using Rasters, Glob

saved_output = glob("*.nc", SIMULATION_PATH)
isothermal_lp, stable_lp, verystable_lp = saved_output
