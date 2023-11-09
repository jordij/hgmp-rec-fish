config <- yaml.load_file("./config.yaml", eval.expr=TRUE)

rec_fish <- raster(paste(config$data, config$rec.fishing, sep="/"))

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
