# hgmp-rec-fish

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

