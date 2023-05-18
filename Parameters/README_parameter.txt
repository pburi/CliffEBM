#===============================================================================
# README PARAMETERS
#===============================================================================
HEADERS:

Run	ad	ai	ed	ei	ft	lw	rh	sw	ta	ts	ws	z0



VARIABLES:

<< Run >>
ID-number for sensitivity run

<< ad >> 
ice albedo [-]

<< ai >> 
debris albedo [-]

<< ed >> 
debris emissivity [-]

<< ei >> 
ice emissivity [-]

<< ft >> 
value for filtering (smoothing) of dem1 in terms of aspect:
"filter_val x filter_val" mean filter around pixel (e.g. '3')

<< lw >> 
longwave radiation deviation [W/m2]

<< rh >> 
relative humidity deviation [%]

<< sw >> 
shortwave radiation deviation [W/m2]

<< ta >> 
T air deviation [K]

<< ts >> 
T surface deviation [K]

<< ws >> 
wind speed multiplier [-]

<< z0 >> 
surface roughness [m]



EXAMPLE STRUCTURE:

Run	ad	ai	ed	ei	ft	lw	rh	sw	ta	ts	ws	z0
1	0.15	0.2	0.95	0.97	9	0	0	0	0	0	1	0.003