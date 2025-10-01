##############################################################
# Runs the Zonation prioritisation process for each one of the 
# 847 cells that fall within an HPA.
# Grab a hot drink if you run this..
##############################################################


config <- yaml.load_file("./config.yaml", eval.expr=TRUE)
habitats <- raster(paste(config$data, config$habitats, sep="/"))


for (hpa.name in unique(hpas$NAME)) {
    cells.hpa <-
        cellnumbers(rec_fish, hpas[hpas$NAME == hpa.name,])$cell_
    for (cell in cells.hpa) {
        output_dir <- paste(config$zonation.cells, as.character(cell), sep = "/")
        cmd.file <- paste(output_dir, "/", as.character(cell) ,".cmd", sep = "")
        shell(normalizePath(cmd.file), wait = TRUE)
    }
}
