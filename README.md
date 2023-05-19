# CliffEBM
A Gridded Ice Cliff Energy Balance Model

This repository provides a model code that calculates the distributed surface energy balance for ice cliffs, i.e. steep ice surfaces with complex, heterogeneous topographies. It was originally designed for supraglacial ice cliffs on debris-covered glaciers.
The model has been applied on debris-covered glaciers in the Nepalese Himalayas (Buri et al., 2016a) and in Southeast Tibet (Kneib et al., 2022).



================ Model example run ================

Here we provide example input data (digital elevation models, shapefiles, meteodata) to run the CliffEBM on one supraglacial cliff on the debris-covered Lirung Glacier (Nepal).
To run the model you need to download the entire branch on your machine, and adjust the paths in the model code (CliffEBM.R, section "primary definitions") according to the paths on your machine.

The model was run on the newest R versions (4.3.0), but should also run on older versions, too.



================ Versions ================

R: R version 4.3.0 (2023-04-21 ucrt) -- "Already Tomorrow"


RStudio: RStudio 2023.03.0+386 "Cherry Blossom" Release (3c53477afb13ab959aeb5b34df1f10c237b256c3, 2023-03-09) for Windows



================ References ================

Buri P, Pellicciotti F, Steiner JF, Miles ES and Immerzeel WW (2016a) A grid-based model of backwasting of supraglacial ice cliffs on debris-covered glaciers. Annals of Glaciology, 57(71), 199–211 (doi: 10.3189/2016AoG71A059)

Kneib M, Miles ES, Buri P, Fugger S, McCarthy MJ, Shaw TE and Pellicciotti F (2022) Sub-seasonal variability of supraglacial ice cliff melt rates and associated processes from time-lapse photogrammetry. The Cryosphere, 16, 4701–4725 (doi: 10.5194/tc-16-4701-2022)
