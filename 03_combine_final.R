config <- yaml.load_file("./config.yaml", eval.expr=TRUE)

rec_fish <- raster(paste(config$data, config$rec.fishing, sep="/"))
hpas <- read_sf(paste(config$data, config$hpas, sep="/"))

files <- list.files(config$zonation.cells, pattern="*rankmap.tif$", full.names=TRUE, recursive = TRUE)
rs <- stack(files)
rs.sum <- calc(rs, fun=sum, na.rm = TRUE)
writeRaster(rs.sum, paste(config$output, "sc_cells_combined.tif", sep=""), datatype="FLT4S")

# Top 847 cells
t <- as_tibble(rs.sum[rs.sum>0,]) %>% arrange(-value)
threshold <- min(t[1:847,]$value)

# Use min value of top 847 cells to filter final added raster rs.sum
tmpfilter <- rs.sum < threshold
filtered.rs.sum <- mask(rs.sum, tmpfilter, maskvalue=TRUE)
writeRaster(filtered.rs.sum, paste(config$output, "sc_cells_combined_top_847.tif", sep=""), datatype="FLT4S")


# Find max for each raster and then combine them

# first initialise raster using the HG mask
hg.mask <- raster(paste(config$data, config$mask, sep="/"))
hg.mask <- setValues(hg.mask, 0)

# Get cell with value == 1 for each raster
idxs <- sapply(unstack(rs), function(r){which(values(r)==1)})
# find cells' indexes
uidxs <- unique(idxs)
for (i in uidxs) {
    tryCatch(hg.mask[i] <- 2 , error= function(e) {return(0)})
}
writeRaster(hg.mask, paste(config$output, "sc_cells_highest.tif", sep=""), datatype="FLT4S")


hg.mask <- raster(paste(config$data, config$mask, sep="/"))
hg.mask <- setValues(hg.mask, 0)

# For each cell, find cell with highest value in its prioritisation result
# (Zonation run) and then fetch value from Rec Fishing layer. Accumulate values
# for all looped cells and save results in output raster.

for (hpa.name in unique(hpas$NAME)) {
    cells.hpa <-
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name,])$cell_
    
    for (cell in cells.hpa) {
        output_dir <- paste(config$zonation.cells, as.character(cell), sep = "/")
        rankmap <- paste(output_dir, "output", "rankmap.tif", sep = "/")
        cell.rec.fish.val <- raster::extract(rec_fish, cell)
        
        cell.raster <- raster(rankmap)
        i <- which(values(cell.raster)==1)
        tryCatch(hg.mask[i] <- hg.mask[i] + rec_fish[i] , error= function(e) {return(0)})
    }
}

writeRaster(hg.mask, paste(config$output, "sc_cells_highest_displacement_calculated.tif", sep=""), datatype="FLT4S")
