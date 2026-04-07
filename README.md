R code and data used for the manuscript:

Tablada, J., Bennion, M., Geange, S.W., Duffy, C.A., Stephenson, F., Hartill, B.W., Hiddink, J.G. and Stephenson, F., 2025. Spatial assessment of recreational fishing displacement from marine protected areas. *Ecosphere* 


All datasets used in this study are publicly available at the [Department of Conservation Marine Data Portal](https://doc-marine-data-deptconservation.hub.arcgis.com/) and the [Ministry for Primary Industries Open Data Portal](https://data-mpi.opendata.arcgis.com/). 

The analysis is divided into 3 scripts:

`01_main.R`

- Reads all relevant data and generates all figures.
- Calculates distance rasters for each of the boat ramps.
- Uses distance rasters to create masks for each of the 847 cells.
- Creates a Zonation setup for each one of the 847 cells.

`02_cells_zonation.R`

- Executes all 847 Zonation scenarios in sequence.

`03_combine_final.R`

- Reads all 847 Zonation prioritisation rasters (outputs from previous step), combines them (func=sum) and 
filters the top 847 cells, which is then saved as a raster file.

