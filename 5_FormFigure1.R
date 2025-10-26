library(ggplot2)
library(ggpubr)

fig1a=readRDS("PlotData/figure1a.RDS")
fig1b=readRDS("PlotData/figure1b.RDS")
fig1c=readRDS("PlotData/figure1c.RDS")
fig1d=readRDS("PlotData/figure1d.RDS")

fig1c=fig1c + theme(legend.position = "none") 
fig1d=fig1d + theme(legend.position = "none") 

# Arrange the first pair (top row)
bottom_row <- ggarrange(fig1c,  fig1d, ncol = 2, nrow = 1,
                     widths = c(2, 0.8), labels   = c("C", "D"),
                     common.legend = TRUE,  
                     legend = "right")

# Arrange the second pair (bottom row)
top_row <- ggarrange(fig1a,  fig1b, ncol = 2, nrow = 1,
            widths = c(2, 0.8), 
            labels   = c("A", "B"), common.legend = TRUE,
            legend = "right")

# Now, combine the two arranged plots
final_figure <- ggarrange(top_row, bottom_row, nrow = 2)


ggsave(final_figure,
       file="ManuscriptFigures/Figure1.tiff", width=12,
       height=10, bg="white", dpi=500)
