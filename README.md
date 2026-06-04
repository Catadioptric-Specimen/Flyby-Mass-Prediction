# Flyby-Mass-Prediction
Satellite orbit solver with dynamic timesteps to constrain the error in energy within a specified tolerance.
Also predicts the mass of the asteroid based on the probe velocity data by fitting a circle.

# How to use:
Fill the Asteroid_Mass and Asteroid_Location arrays with the masses and locations of spherical immobile asteroids
you want to flyby (SI units).

Fill the entries of Probe_Coords with the initial position and velocity of the probe (SI units).

Set tf (end time in seconds) and dt (timestep size in seconds).

If you wish, you can also adjust the error tolerances (note: often it is reccomended to set p_dEnergy to 0)

# Outputs:
The console returns:
1. The fitted hodograph circle distance from origin and radius
2. Estimated asteroid mass
3. Error in mass estimate
4. Relative fluctuation in probe energy throughout the simulation

The plots returned are:
1. Asteroid (black) and flyby trajectory of the probe (blue)
2. Plot of velocity vectors (blue) with fitted circle (red)
3. Plot of probe kinetic and potential energy
4. Plot of total probe mechanical energy
