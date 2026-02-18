library(readxl)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
library(tidyr)


fig1a <- readRDS("PlotData/supp_figure1a.RDS") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 5, margin = margin(2,0,2,0))
  ) + labs(x = "Drug")

fig1b <- readRDS("PlotData/supp_figure1b.RDS") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 5, margin = margin(2,0,2,0))
  ) + labs(x = "Drug")

top_row <- ggarrange(fig1a, fig1b, 
                     ncol = 2, nrow = 1,
                     widths = c(2, 0.8), 
                     labels = c("A", "B"),
                     common.legend = TRUE,
                     legend = "bottom",
                     align = "h")

vene_comb <- read_excel("DataSets/bcd-23-0014_table_s9_suppst9.xlsx", col_types='text')
other_drugs <- sort(unique(vene_comb$drug[vene_comb$drug_type == "single agent" & vene_comb$drug != "VEN"]))
combos <- paste0("VEN + ", other_drugs)

interleaved_levels <- as.vector(rbind(other_drugs, combos))
drug_levels <- c("VEN", interleaved_levels)

plot_data <- vene_comb %>%
  filter(drug %in% drug_levels) %>%
  mutate(
    probit_auc = as.numeric(probit_auc),
    facet_group = case_when(
      drug == "VEN" ~ "Reference",
      grepl("VEN \\+", drug) ~ gsub("VEN \\+ ", "", drug),
      TRUE ~ drug
    ),
    group_color = case_when(
      drug == "VEN" ~ "Venetoclax",
      drug %in% other_drugs ~ "Single Agent",
      TRUE ~ "Combination"
    )
  )

plot_data$facet_group <- factor(plot_data$facet_group, levels = c("Reference", other_drugs))
plot_data$drug <- factor(plot_data$drug, levels = drug_levels)

safe_colors <- c(
  "Venetoclax" = "#56B4E9", 
  "Single Agent" = "#D55E00", 
  "Combination" = "#009E73"
)

plt_c <- ggplot(plot_data, aes(x = drug, y = probit_auc, fill = group_color)) +
  geom_boxplot(width = 0.75, outlier.size = 0.4) +
  facet_grid(. ~ facet_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = safe_colors) +
  theme_bw() +
  labs(y = "Area Under the Curve (AUC)", x = "Drug", fill = "Treatment Type") +
  theme(
    # Shrink header (strip) space
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    panel.spacing = unit(0.2, "lines"),
    # Shrink Axis Fonts
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7),
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    panel.grid.major.x = element_blank()
  )

final_fig <- ggarrange(top_row, plt_c, 
                       nrow = 2, 
                       labels = c("", "C"), 
                       heights = c(1, 1.2))
cairo_pdf("../Plots/Supp_figure1.pdf", width = 14, height = 10)
print(final_fig)
dev.off()
