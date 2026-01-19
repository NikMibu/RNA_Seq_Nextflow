library(ggplot2)
library(dplyr)

setwd("/mnt/d/RNA_seq/results/rmats/rmats_output")

target_genes <- c('"ABCC5"', '"CRNDE"', '"UQCC"', '"UQCC1"', 
                  '"GUSBP11"', '"ANKHD1"', '"ADAM12"')

read_rmats <- function(file, event_type) {
  if (!file.exists(file)) return(NULL)
  data <- read.table(file, header = TRUE, sep = "\t", 
                     stringsAsFactors = FALSE, quote = "")
  data %>%
    filter(geneSymbol %in% target_genes, FDR < 0.05) %>%
    mutate(EventType = event_type, Gene = gsub('"', '', geneSymbol))
}

cat("Reading rMATS files...\n")
all_events <- bind_rows(
  read_rmats("SE.MATS.JC.txt", "SE"),
  read_rmats("A3SS.MATS.JC.txt", "A3SS"),
  read_rmats("A5SS.MATS.JC.txt", "A5SS"),
  read_rmats("MXE.MATS.JC.txt", "MXE"),
  read_rmats("RI.MATS.JC.txt", "RI")
) %>% filter(!is.na(Gene))

cat(paste0("Found ", nrow(all_events), " significant events\n\n"))

# Summary
summary_table <- all_events %>%
  group_by(Gene, EventType) %>%
  summarise(n_events = n(), mean_delta = mean(IncLevelDifference),
            min_FDR = min(FDR), .groups = "drop")
print(summary_table)

# Export table
write.table(all_events %>% select(Gene, EventType, ID, FDR, PValue, IncLevelDifference),
            "signifikante_events_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot 1: Overview
p1 <- ggplot(all_events %>% count(Gene, EventType), 
       aes(x = reorder(Gene, -n), y = n, fill = EventType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") +
  theme_minimal(base_size = 14) +
  labs(title = "Significant Splicing Events (FDR < 0.05)",
       subtitle = "SF3B1 mutant vs wildtype",
       x = "Gene", y = "Number of Events") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

ggsave("signifikante_events_count.pdf", p1, width = 8, height = 6)
cat("✅ signifikante_events_count.pdf\n")

# Plot 2: Gene-specific plots
for (gene in unique(all_events$Gene)) {
  gene_data <- all_events %>% 
    filter(Gene == gene) %>%
    mutate(event_id = paste(EventType, 1:n()))  # FIX: Create ID inside mutate
  
  p <- ggplot(gene_data, aes(x = event_id, y = IncLevelDifference, fill = EventType)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("FDR=%.2e", FDR)), 
              vjust = ifelse(gene_data$IncLevelDifference > 0, -0.5, 1.5), size = 3.5) +
    theme_minimal(base_size = 14) +
    labs(title = paste0(gene, " - Alternative Splicing"),
         subtitle = paste0(nrow(gene_data), " significant events"),
         x = "Event", y = "ΔPSI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  
  ggsave(paste0(gene, "_events.pdf"), p, width = 10, height = 6)
  cat(paste0("✅ ", gene, "_events.pdf\n"))
}

cat("\n🎉 Done!\n")
