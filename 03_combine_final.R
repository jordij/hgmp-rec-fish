##############################################################
# Combines all cells Zonation outputs and generates a few
# useful rasters to filter final results to cells with the
# highest values overall.
# Grab another hot drink if you run this..
##############################################################

config <- yaml.load_file("./config.yaml", eval.expr = TRUE)

rec_fish <- raster(paste(config$data, config$rec.fishing, sep = "/"))
hpas <- read_sf(paste(config$data, config$hpas, sep = "/"))

files <-
    list.files(
        config$zonation.cells,
        pattern = "*rankmap.tif$",
        full.names = TRUE,
        recursive = TRUE
    )
rs <- stack(files)
rs.sum <- calc(rs, fun = sum, na.rm = TRUE)
writeRaster(
    rs.sum,
    paste(config$output, "sc_cells_combined.tif", sep = ""),
    overwrite = TRUE,
    datatype = "FLT4S"
)

# Top 847 cells
t <- as_tibble(rs.sum[rs.sum > 0, ]) %>% arrange(-value)
threshold <- min(t[1:847, ]$value)

# Use min value of top 847 cells to filter final added raster rs.sum
tmpfilter <- rs.sum < threshold
filtered.rs.sum <- mask(rs.sum, tmpfilter, maskvalue = TRUE)
writeRaster(
    filtered.rs.sum,
    paste(config$output, "sc_cells_combined_top_847.tif", sep = ""),
    overwrite = TRUE,
    datatype = "FLT4S"
)

# do the same but use a % (let's say keep the top 30%)
threshold <- min(t[1:(nrow(t) * 0.3), ]$value)
tmpfilter <- rs.sum < threshold
filtered.rs.sum <- mask(rs.sum, tmpfilter, maskvalue = TRUE)
writeRaster(
    filtered.rs.sum,
    paste(config$output, "sc_cells_combined_top_30percent.tif", sep = ""),
    overwrite = TRUE,
    datatype = "FLT4S"
)

# Generate 10% top for each HPA but raster values are normalised 0-1
for (hpa.name in unique(hpas$NAME)) {
    # Get cell ID
    cells.hpa <-
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name, ])$cell_
    # Get all zonation cells rasters
    hpa.files <- lapply(cells.hpa, function(cell) {
        cell.str <- paste("/", cell, "/", sep = "")
        files[grepl(cell.str, files)]
    })
    # Merge and filter top 10%
    rs.hpa <- stack(hpa.files)
    rs.hpa.sum <- calc(rs.hpa, fun = sum, na.rm = TRUE)
    t <- as_tibble(rs.hpa.sum[rs.hpa.sum > 0, ]) %>% arrange(-value)
    threshold <- min(t[1:(nrow(t) * 0.1), ]$value)
    tmpfilter <- rs.hpa.sum < threshold
    filtered.rs.sum <- mask(rs.hpa.sum, tmpfilter, maskvalue = TRUE)
    # Normalise to 0-1 using `val = (val - min / max - min)`
    minr <- cellStats(filtered.rs.sum, stat = min)
    maxr <- cellStats(filtered.rs.sum, stat = max)
    filtered.rs.sum <- (filtered.rs.sum - minr) / (maxr - minr)
    #
    # writeRaster(
    #     filtered.rs.sum,
    #     paste(
    #         config$output,
    #         paste("top10percent_", hpa.name, ".tif", sep = ""),
    #         sep = ""
    #     ),
    #     overwrite = TRUE,
    #     datatype = "FLT4S"
    # )
}

rm(hpa.files, rs.hpa, rs.hpa.sum)
# Find max for each raster and then combine them

# first initialise raster using the HG mask
hg.mask <- raster(paste(config$data, config$mask, sep = "/"))
hg.mask <- setValues(hg.mask, 0)

# Get cell with value == 1 for each raster
idxs <- sapply(unstack(rs), function(r) {
    which(values(r) == 1)
})
# find cells' indexes
uidxs <- unique(idxs)
for (i in uidxs) {
    tryCatch(
        hg.mask[i] <- 2 ,
        error = function(e) {
            return(0)
        }
    )
}
writeRaster(
    hg.mask,
    paste(config$output, "sc_cells_highest.tif", sep = ""),
    overwrite = TRUE,
    datatype = "FLT4S"
)


hg.mask <- raster(paste(config$data, config$mask, sep = "/"))
hg.mask <- setValues(hg.mask, 0)

# For each cell, find cell with highest value in its prioritisation result
# (Zonation run) and then fetch value from Rec Fishing layer. Accumulate values
# for all looped cells and save results in output raster.

for (hpa.name in unique(hpas$NAME)) {
    cells.hpa <-
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name, ])$cell_
    
    for (cell in cells.hpa) {
        output_dir <-
            paste(config$zonation.cells, as.character(cell), sep = "/")
        rankmap <-
            paste(output_dir, "output", "rankmap.tif", sep = "/")
        cell.rec.fish.val <- raster::extract(rec_fish, cell)
        
        cell.raster <- raster(rankmap)
        i <- which(values(cell.raster) == 1)
        tryCatch(
            hg.mask[i] <-
                hg.mask[i] + rec_fish[i] ,
            error = function(e) {
                return(0)
            }
        )
    }
}

writeRaster(
    hg.mask,
    paste(
        config$output,
        "sc_cells_highest_displacement_calculated.tif",
        sep = ""
    ),
    overwrite = TRUE,
    datatype = "FLT4S"
)
