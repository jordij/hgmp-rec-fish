library(ggpubr)
library(glue)
library(raster)
library(sf)
library(tabularaster)
library(terra)
library(tidyverse)
library(yaml)

config <- yaml.load_file("./config.yaml", eval.expr=TRUE)

# Read all necessary layers
rec_fish <- raster(paste(config$data, config$rec.fishing, sep="/"))
hpas <- read_sf(paste(config$data, config$hpas, sep="/"))
mrs <- read_sf(paste(config$data, config$mrs, sep="/"))
boat.ramps <- shapefile(paste(config$data, config$boat.ramps, sep="/"))
hg.mask <- raster(paste(config$data, config$mask, sep="/"))
habitats <- raster(paste(config$data, config$habitats, sep="/"))

# Load zonation outputs
zon.sc.bio <- raster(paste(config$zonation.bio, "output", "rankmap.tif", sep="/"))
zon.rec.fish <- raster(paste(config$zonation.recfish, "output", "merged_rankmap.tif", sep="/"))

# Pearson for biodiversity (SDMs) against rec fishing zonation output
rstack <- raster::stack(zon.sc.bio, zon.rec.fish)
stats.pearson <- layerStats(rstack, stat = "pearson", na.rm = TRUE)
print(stats.pearson)
rm(stats.pearson, rstack)

# Plot Pearson's correlation raster (R and P-values)
cor.ras <- corLocal(zon.sc.bio, zon.rec.fish, test = T)
png(paste(config$figs, "correlation_raster.png", sep=""), width = 1200, height = 800)
plot(cor.ras)
dev.off()
rm(cor.ras)


# Zonation rec fishing and biodiversity to dataframes
df <- cbind.data.frame(values(zon.rec.fish), values(zon.sc.bio))
names(df) <- c("Biodiversity", "Rec fishing")
df.clean <- df[!is.na(df$Biodiversity) & !is.na(df$`Rec fishing`), ]

# random 500 points, 100 iterations
plot.df <- data.frame(Method=character(0), Value=numeric(0))

for (i in 1:100) {
    random.idx <- sample(nrow(df.clean), 500)
    sdf <- df.clean[random.idx,]
    p <- cor(sdf$Biodiversity, sdf$`Rec fishing`, method = "pearson")
    k <- cor(sdf$Biodiversity, sdf$`Rec fishing`, method = "kendall")
    s <- cor(sdf$Biodiversity, sdf$`Rec fishing`, method = "spearman")
    plot.df <- plot.df %>% 
        add_row(Method = "pearson", Value=p) %>%
        add_row(Method = "kendall", Value=k) %>%
        add_row(Method = "spearman", Value=s)
}

plot.df %>%
    group_by(Method) %>%
    summarise(total = n(), min=min(Value), max=max(Value)) 

rm(sdf, df.clean)
rm(p,k,s)

# Plot correlations using cells where the original rec fishing value was greater
# than 1 (which is 0-0.01 vessels/km2)

srf <- rec_fish * (rec_fish$VesselsPer > 1)
srf[srf==0] <- NA
sub.zon.rec.fish <- mask(zon.rec.fish, srf)
sub.zon.sc.bio <- mask(zon.sc.bio, srf)
rstack <- raster::stack(sub.zon.sc.bio, sub.zon.rec.fish)
stats.pearson <- layerStats(rstack, stat = "pearson", na.rm = TRUE)
print(stats.pearson)

df <- cbind.data.frame(values(sub.zon.rec.fish), values(sub.zon.sc.bio))
names(df) <- c("Biodiversity", "Rec fishing")
df.clean <- df[!is.na(df$Biodiversity) & !is.na(df$`Rec fishing`), ]

# random 500 points, 100 iterations
plot.df <- data.frame(Method=character(0), Value=numeric(0))

for (i in 1:100) {
    random.idx <- sample(nrow(df.clean), 500)
    sdf <- df.clean[random.idx,]
    p <- cor(sdf$Biodiversity, sdf$`Rec fishing`, method = "pearson")
    s <- cor(sdf$Biodiversity, sdf$`Rec fishing`, method = "spearman")
    plot.df <- plot.df %>% 
        add_row(Method = "pearson", Value=p) %>%
        add_row(Method = "spearman", Value=s)
}

plot.df %>%
    group_by(Method) %>%
    summarise(total = n(), min=min(Value), max=max(Value)) 

rm(sdf, df.clean)
rm(p,k,s)

# Create and save boxplot with correlation values

corr.plot <- ggplot(plot.df, aes(x=Method, y=Value)) + 
    geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, width=0.2)
ggsave(
    paste(config$figs, "correlation_boxplot_filtered.png", sep=""),
    corr.plot,
    width=4, height=5,
    dpi="retina"
)

maxlim <- max(plot.df$Value)
minlim <- min(plot.df$Value)

p <- ggplot(plot.df, aes(x=Method, y=Value, fill = Method)) +
    geom_violin(color="black") + # OR geom_boxplot() 
    #geom_point(position=position_jitter(0.06), alpha = 3/10, color="black") + 
    geom_jitter(color="black", size=0.4, width=0.2) +
    theme_bw() +
    ylim(minlim, maxlim) +
    theme(axis.text = element_text(size = 12)) + 
    theme(axis.title = element_text(size = 14)) + 
    theme(legend.title = element_text(size = 16)) + 
    theme(legend.text = element_text(size = 12)) +
    scale_fill_manual(values=c("#66c2a5","#fc8d62")) +
    ylab("Correlation coefficient") +
    xlab("Method") +
    theme(legend.position="none")

ggsave(
    paste(config$figs, "correlation_boxplot_filtered_v2.png", sep=""),
    p,
    width=4, height=5,
    dpi="retina"
)

rm(corr.plot)
rm(plot.df)

# Get cells that fall within the HPAs

cells.intersect <- cellnumbers(rec_fish, hpas)$cell_
df.intersect <- df[cells.intersect,]
ggscatter(
    df.intersect,
    x = "Biodiversity",
    y = "Rec fishing",
    add = "reg.line",
    conf.int = TRUE
) + stat_cor(method = "pearson")

# Plots for each one of the HPAs
for (hpa.name in unique(hpas$NAME)) {
    cells.hpa <- cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name,])$cell_
    df.intersect <- df[cells.hpa,]
    hpa.plot <- ggscatter(df.intersect,
                          x = "Biodiversity",
                          y = "Rec fishing",
                          title = hpa.name) +
        stat_cor(method = "pearson") 
    
    ggsave(
        paste(config$figs, hpa.name, ".png", sep=""),
        hpa.plot,
        width=6, height=4,
        dpi="retina"
    )
}
rm(hpa.plot)
rm(df.intersect, df)
rm(cells.hpa)
rm(cells.intersect)

# For each boat ramp, create a raster with the grid distance
# from all cells to the given boat ramp location
# Save each raster file to the output folder.
# Go get a cuppa..

df.boat.ramps <- data.frame(boat.ramps)
df.boat.ramps$file <- NA

for (i in 1:nrow(df.boat.ramps)) {
    cell.idx <- cellFromXY(hg.mask, boat.ramps[i, ])
    dest.file <- paste(
        config$output,
        config$distance.ramps,
        "dist_grid_",
        as.character(i),
        gsub("[^[:alnum:] ]", "", df.boat.ramps[i, ]$title),
        ".tif",
        sep = ""
    )
    df.boat.ramps[i, c("file")] <- dest.file
    val <- hg.mask[cell.idx]
    hg.mask[cell.idx] <- 2
    dist.grid <- gridDistance(hg.mask, 2, omit = NA)
    writeRaster(dist.grid, dest.file, overwrite = TRUE)
    # re-initialise original mask
    hg.mask <- raster(paste(config$data, config$mask, sep="/"))
}

rm(dest.file, dist.grid)

# Stack individual boat ramp distance rasters
grid.distances.files <-
    list.files(
        path = paste(config$output, config$distance.ramps, sep=""),
        pattern = "dist_grid_",
        all.files = FALSE,
        full.names = TRUE,
        recursive = TRUE
    )
stacked.grid.distances <- stack(grid.distances.files)
names(stacked.grid.distances) <- grid.distances.files


getCellsWithinDistance <- function(i, files, distances, cell) {
    ramp.raster <- raster(files[distances$ix[i]])
    return (ramp.raster <= ramp.raster[cell])
}

for (hpa.name in unique(hpas$NAME)) {
    cells.hpa <-
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name,])$cell_
    for (cell in cells.hpa) {
        # Get closest N (config$n.ramps) boat ramps to current cell and load the
        # boat ramps distance rasters. Then mask out all cells with distance
        # greater than current cell. Merge result rasters and save to file.
        dist.br <- raster::extract(stacked.grid.distances, cell)
        dist.br <- sort(dist.br, index.return = T, na.last = T)
        masks <-
            lapply(
                1:config$n.ramps,
                FUN = getCellsWithinDistance,
                files = grid.distances.files,
                distances = dist.br,
                cell=cell
            )

        masked.cell.area <- merge(stack(masks)) %>%
            mask(hpas, inverse = TRUE) %>%
            mask(mrs, inverse = TRUE)
        dest.file <- paste(
            config$output,
            config$masked.cells,
            "masked_cell_",
            as.character(cell),
            ".tif",
            sep = ""
        )
        # write raster (int)
        writeRaster(masked.cell.area,
                    dest.file,
                    overwrite = TRUE,
                    datatype = 'INT2S')
    }
}

rm(masked.cell.area, dest.file, dist.br, masks)

##################################################
# Setup Zonation-style folder folder for each cell
# 
# hgmp-rec-fish/sc_cells/cell_{ID}/
# |---> {ID}.cmd
# |---> features.txt
# |---> settings.z5
# |---> output/
#
##################################################

lvls <- levels(habitats)[[1]]

for (hpa.name in unique(hpas$NAME)) {
    cells.hpa <-
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name,])$cell_

    for (cell in cells.hpa) {
        habitatID <- raster::extract(habitats, cell)
        habitat <- gsub(" ", "_", lvls[lvls$ID == habitatID,]$Habitat,)
        habitat.path <- paste(
            "../.",
            config$data,
            "/",
            config$habitats,
            "/",
            habitatID, sep = "")
        mask.path <- paste(
            "../.",
            config$output,
            config$masked.cells,
            "masked_cell_",
            as.character(cell),
            ".tif",
            sep = ""
        )
        output_dir <- paste(config$zonation.cells, as.character(cell), sep = "/")
        # Create cell folder with nested settings and .cmd file
        if (!dir.exists(output_dir)) {
            dir.create(output_dir)
            features.file <- paste(output_dir, "/features.txt", sep = "")
            settings.file <- paste(output_dir, "/settings.z5", sep = "")
            cmd.file <- paste(output_dir, "/", as.character(cell) ,".cmd", sep = "")
            output_dir <- paste(output_dir, "output", sep = "/")
            dir.create(output_dir)
            # features file
            file.copy(from="./features_template.txt", to=features.file)
            # add a new line for physical habitat
            cat(paste("\n1.0", habitat.path, sep=" "), file = features.file, append = TRUE)
            # settings file
            file.copy(from="./settings_template.z5", to=settings.file)
            cat(paste("\nanalysis area mask layer = ", mask.path, sep=""), file = settings.file, append = TRUE)
            # bat file
            bat.content <- paste("start \"\" ", config$zonation.exe, '--mode=CAZ2 -wa')
            bat.content <- paste(bat.content, normalizePath(settings.file), normalizePath(output_dir), sep=" ")
            write(x=bat.content, file = cmd.file)
            
        }
    }
}
