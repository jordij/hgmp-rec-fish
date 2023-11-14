dist.df <- as.data.frame(dist.raster, xy=TRUE)
rec.fish.df <- as.data.frame(rec_fish, xy=TRUE)
dfd <- merge(rec.fish.df, dist.df)

dfd$Density <- factor(dfd$Density , levels=c("0 - 0.01",
                                             "0.01 - 1 1.01",
                                             "1.01 - 2",
                                             "2.01 - 5",
                                             "5.01 - 10",
                                             "10.01 - 15",
                                             "15.01 - 20",
                                             "20.01 - 50",
                                             "50.01 - 100",
                                             "100.1 - 150",
                                             "150.01 - 200",
                                             "200.01 - 250",
                                             "250.01 - 300",
                                             "300.01 - 400"))

dfd$Distance <- dfd$Distance/1000

dist.dens.plot <- ggplot(dfd[!is.na(dfd$Density) & !is.na(dfd$Distance),], aes(x=Density, y=Distance)) + 
    geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, width=0.2) +
    xlab("Vessel density")
    

ggsave(
    paste(config$figs, "density_distance_boxplot.png", sep=""),
    dist.dens.plot,
    width=12, height=4,
    dpi="retina"
)