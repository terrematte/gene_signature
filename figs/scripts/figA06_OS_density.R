library(tidyverse)

mu <- ddply(data, "dataset", summarise, grp.mean=mean(obs.time))
mu$years <- c(7, 10) 
mu$text_years <- c("7-years", "10-years")
mu$y  <- c(0.15, 0.15)

p1 <- ggplot(data, aes(x = obs.time/365.2425, fill = dataset, colour = dataset)) + 
        geom_vline(data=mu, aes(xintercept=years, color=dataset),
             linetype="solid") +
        geom_text(mu, mapping=aes(x = years, y = y, label = text_years, group = dataset),  vjust = -0.5, size = 4) +
        geom_histogram(alpha = 0.5, position = "identity") +
        xlab("Overall Survival (years) after diagnose") +
        ylab("Density")

  
  p2 <- ggplot(data, aes(x= obs.time/365.2425, fill=dataset)) +
    geom_density(alpha=0.4) +
    geom_vline(data=mu, aes(xintercept=grp.mean/365.2425, color=dataset),
             linetype="dotted", show.legend = FALSE) +
    geom_vline(data=mu, aes(xintercept=years, color=dataset),
               linetype="solid", show.legend = FALSE) +
    geom_label(mu, mapping=aes(x = years, y = y, label = text_years, group = dataset),  vjust = -0.5, size = 4) +
    guides(
      fill = guide_legend(
        title = "Datasets",
        override.aes = aes(label = "")
      )
    ) +
    xlab("Overall Survival (years) after diagnose") +
    ylab("Density")
  
  library(cowplot)
  
  pdf(file="figs/figA06_OS_density.pdf", width = 8, height = 4)
  cowplot::plot_grid(
    p2
  )
  dev.off()

# pdf(file="figs/figA06_OS_distribuition.pdf", width = 12, height = 4)
# cowplot::plot_grid(
#   p1, p2,
#   nrow = 1,
#   labels = c("a","b")
#   #rel_widths = c(0.1, 1, 0.1)
# )
# dev.off()