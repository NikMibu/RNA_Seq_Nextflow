# ══════════════════════════════════════════════════════════
# Venn Diagram - Vergleich DEGs mit Paper
# ══════════════════════════════════════════════════════════

library(ggplot2)

# Deine gefundenen Gene (aus rMATS Analyse)
# Quelle: /mnt/d/RNA_seq/results/rmats/analysis/signifikante_events_table.txt
meine_gene <- c("ABCC5", "CRNDE", "ANKHD1")

# Paper Gene (Furney et al. - Table 2)
paper_gene <- c("ABCC5", "CRNDE", "UQCC", "GUSBP11", "ANKHD1", "ADAM12")

# Zählungen für 2-Set-Venn
nur_meine <- length(setdiff(meine_gene, paper_gene))   # 0
nur_paper <- length(setdiff(paper_gene, meine_gene))   # 3: UQCC, GUSBP11, ADAM12
beide     <- length(intersect(meine_gene, paper_gene)) # 3: ABCC5, CRNDE, ANKHD1

# Kreise als Polygon (Punkte auf Kreisumfang)
circle <- function(cx, cy, r, n = 100) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}

# Zwei sich überlappende Kreise (zentriert für 2-Set-Venn)
d1 <- circle(-0.6, 0, 1.1)
d1$set <- "Meine Analyse"
d2 <- circle( 0.6, 0, 1.1)
d2$set <- "Paper (Furney et al.)"

p <- ggplot() +
  geom_polygon(data = d1, aes(x, y, fill = "Meine Analyse"), alpha = 0.4) +
  geom_polygon(data = d2, aes(x, y, fill = "Paper (Furney et al.)"), alpha = 0.4) +
  annotate("text", x = -1.0, y = 0,   label = nur_meine, size = 6, fontface = "bold") +
  annotate("text", x =  1.0, y = 0,   label = nur_paper, size = 6, fontface = "bold") +
  annotate("text", x =  0,   y = 0,   label = beide,     size = 6, fontface = "bold") +
  scale_fill_manual(values = c("Meine Analyse" = "steelblue", "Paper (Furney et al.)" = "coral")) +
  labs(title = "Überlappung: Differenziell gesplicte Gene (FDR < 0.05)",
       subtitle = "Meine Analyse vs. Paper (Furney et al.)") +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  coord_fixed()

ggsave("/mnt/d/RNA_seq/results/Venn_Splicing_Genes.pdf", p, width = 6, height = 5)
