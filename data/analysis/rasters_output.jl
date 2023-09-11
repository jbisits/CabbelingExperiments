using Rasters, Glob

saved_output = glob("*.nc", SIMULATION_PATH)
cab_noise, isothermal_lp, stable_noise, stable_lp, verystable_lp = saved_output
