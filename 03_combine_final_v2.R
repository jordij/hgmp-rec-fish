##############################################################
# Combines all cells Zonation outputs and generates a few
# useful rasters to filter final results to cells with the
# highest values overall.
# Grab another hot drink if you run this..
##############################################################

config <- yaml.load_file("./config.yaml", eval.expr = TRUE)

rec_fish <-
    raster(paste(config$data, config$rec.fishing, sep = "/"))
hpas <- read_sf(paste(config$data, config$hpas, sep = "/"))
hg.mask <- raster(paste(config$data, config$mask, sep = "/"))


files <-
    list.files(
        config$zonation.cells,
        pattern = "*rankmap.tif$",
        full.names = TRUE,
        recursive = TRUE
    )
rs <- stack(files)


# Summed rasters
rs.sum <- calc(rs, fun = sum, na.rm = TRUE)
writeRaster(
    rs.sum,
    paste(config$output, "sc_cells_combined.tif", sep = ""),
    overwrite = TRUE,
    datatype = "FLT4S"
)

# Top 847 cells
t <- as_tibble(rs.sum[rs.sum > 0,]) %>% arrange(-value)
threshold <- min(t[1:847,]$value)

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
threshold <- min(t[1:(nrow(t) * 0.3),]$value)
tmpfilter <- rs.sum < threshold
filtered.rs.sum <- mask(rs.sum, tmpfilter, maskvalue = TRUE)
writeRaster(
    filtered.rs.sum,
    paste(config$output, "sc_cells_combined_top_30percent.tif", sep = ""),
    overwrite = TRUE,
    datatype = "FLT4S"
)

# MPA potential displacement (top 10%)
mpa.displacement.df <- data.frame()
# Generate 10% top for each HPA but with raster values normalised [0-1]
for (hpa.name in unique(hpas$NAME)) {
    # Get cell ID
    cells.hpa <-
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name[1],])$cell_
    # Get all zonation cells rasters
    hpa.files <- lapply(cells.hpa, function(cell) {
        cell.str <- paste("/", cell, "/", sep = "")
        files[grepl(cell.str, files)]
    })
    # Merge and filter top 10%
    rs.hpa <- stack(hpa.files)
    rs.hpa.sum <- calc(rs.hpa, fun = sum, na.rm = TRUE)
    t <- as_tibble(rs.hpa.sum[rs.hpa.sum > 0,]) %>% arrange(-value)
    threshold <- min(t[1:(nrow(t) * 0.1),]$value)
    tmpfilter <- rs.hpa.sum < threshold
    filtered.rs.sum <- mask(rs.hpa.sum, tmpfilter, maskvalue = TRUE)
    # Normalise to 0-1 using `val = (val - min / max - min)`
    minr <- cellStats(filtered.rs.sum, stat = min)
    maxr <- cellStats(filtered.rs.sum, stat = max)
    filtered.rs.sum <- (filtered.rs.sum - minr) / (maxr - minr)
    # Write to file
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
    #
    # To dataframe
    df <- as.data.frame(filtered.rs.sum, xy = TRUE, na.rm = TRUE)
    df <- rename(df, displacement = layer) %>%
        add_column(mpa = hpa.name) %>%
        add_column(distance = 0)
    # calculate distances to mpa cells
    hg.mask[cells.hpa] <- 2
    dist.grid <- gridDistance(hg.mask, 2, omit = NA)
    cells <-
        cellnumbers(filtered.rs.sum, as.matrix(df[, c("x", "y")]))$cell_
    
    df[as.character(cells), ]$distance <- dist.grid[cells]
    mpa.displacement.df <- rbind(mpa.displacement.df, df)
    # re-initialise original mask
    hg.mask <- raster(paste(config$data, config$mask, sep = "/"))
}

mpa.displacement.df$distance <- mpa.displacement.df$distance / 1000

mpa.displacement.df.sum <- mpa.displacement.df %>%
    group_by(mpa) %>%
    summarise(
        distance = mean(distance, na.rm = TRUE),
        displacement = mean(displacement, na.rm = TRUE),
        na.rm = TRUE
    )
#%>%
#    rename(MPA = mpa)

mpa.labels <- c(
    "Cape Rodney Okakari Point Marine Reserve extension" = "Cape Rodney-Okakari Point",
    "Rotoroa Island" = "Rotoroa",                                
    "Whanganui A Hei Marine Reserve extension" = "Whanganui A Hei",          
    "Motukawao" = "Motukawao",                                         
    "Colville Type 1" = "Cape Colville",                                   
    "Tiritiri Type 1" = "Tiritiri Matangi",                                   
    "Slipper" = "Slipper Island",                                           
    "Alderman South" = "Aldermen South",                                     
    "Aldermen North" = "Aldermen North",                                    
    "Little Barrier Island" = "Little Barrier - Hauturu",                             
    "Kawau Bay Type 1" = "Kawau Bay",                               
    "Rangitoto Motutapu" = "Rangitoto - Motutapu",                                
    "Mokohinau Type 1" = "Mokohinau",
    "The Noises" = "Noises"
)

pc <- ggplot(mpa.displacement.df, aes(distance, displacement)) +
    geom_point(alpha = 3 / 10, color = "black") +
    geom_point(data = mpa.displacement.df.sum,
                aes(colour = mpa),
                alpha = 0.9,
                size = 4) +
    facet_wrap( ~ mpa, ncol = 5 ,scales = "fixed", labeller = as_labeller(mpa.labels)) + theme_classic() +
    xlab("Distance from MPA (km)") + ylab("MPA Displacement Potential") +
    theme(legend.position = "none")

ggsave(
    paste(config$figs, "pointcloud_faceted.png", sep=""),
    pc,
    width = 10,
    height = 6,
    dpi = "retina"
)

pc <- ggplot(mpa.displacement.df, aes(distance, displacement)) +
    geom_point(alpha = 3 / 10, aes(colour = mpa)) +
    geom_point(data = mpa.displacement.df.sum,
               aes(colour = MPA),
               size = 4) +
    #facet_wrap(~mpa, scales="fixed") + theme_classic() +
    xlab("Distance from MPA (km)") + ylab("MPA Displacement Potential") +
    theme_classic()

ggsave(
    paste(config$figs, "pointcloud_faceted_all.png", sep=""),
    pc,
    width = 10,
    height = 6,
    dpi = "retina"
)

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
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name,])$cell_
    
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
