library(Rtsne)
library(viridis)
library(RColorBrewer)
library(randomcoloR)
library(patchwork)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(readr)
library(data.table)
library(DESeq2)

setwd("/data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/")

all_metadata <- read_excel("/data2t_2/pathogen_TE_2025_New/01.new_raw_data/all_sample.xlsx")
Top_5_Up_TE_loci <- read.csv("/data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/Recurrent_TE_Loci_gt5_intergenic_n838_closest_genes.bed", 
                             sep = "\t", header = F)
colnames(Top_5_Up_TE_loci) <- c(
  "chr1", "start", "end", "TE_id", "UniqPathogens","TE_location" ,"class", "family","feature",
  "chr2", "gene_start", "gene_end", "gene_width", "strand",  "gene_info","distance"
)

datasets_to_analyze <- c("Virus_EBV_2", "Virus_EBV_5", "Fungi_Af_4", "Fungi_Af_6", "Virus_IAV_11", "Virus_RSV_1","Virus_HCV_1")

convert_time_to_hours <- function(time_point) {
  if (is.na(time_point) || time_point == "NA" || time_point == "0h" || time_point == "0min") {
    return("0h")
  }
  
  time_num <- as.numeric(gsub("min|h", "", time_point))
    if (grepl("h$", time_point)) {
    return(time_point)
  }
  
  if (time_num >= 60) {
    hours <- time_num / 60
    if (hours == floor(hours)) {
      return(paste0(floor(hours), "h"))
    } else {
      return(paste0(round(hours, 1), "h"))
    }
  } else {
    return(paste0(time_num, "min"))
  }
}


process_single_dataset <- function(dataset_name, all_metadata, Top_5_Up_TE_loci) {
  
  if ("Dataset" %in% colnames(all_metadata)) {
    metadata <- all_metadata[all_metadata$Dataset == dataset_name, ]
  } else {
    stop(paste("Cannot find metadata for dataset:", dataset_name))
  }
  
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata$Sample_Name
  if (dataset_name == "Virus_HCV_1") {
    if ("CellType_1" %in% colnames(metadata)) {
      metadata <- metadata[metadata$CellType_1 == "Huh-7", ]
      
      if (nrow(metadata) == 0) {
        warning(paste("No Huh-7 samples found in dataset:", dataset_name))
        return(NULL)
      }
      
      cat("Dataset", dataset_name, ": Selected", nrow(metadata), "Huh-7 samples\n")
    }
    

    if ("Condition" %in% colnames(metadata)) {
      condition_parts <- strsplit(metadata$Condition, "_")
      
      if (!"Infection_State" %in% colnames(metadata) || all(is.na(metadata$Infection_State))) {
        metadata$Infection_State <- sapply(condition_parts, function(x) x[1])
      }
      
      if (!"Time_Point" %in% colnames(metadata) || all(is.na(metadata$Time_Point))) {
        raw_time <- sapply(condition_parts, function(x) x[2])
        metadata$Time_Point <- ifelse(raw_time == "4320m", "72h", 
                                      ifelse(raw_time == "7200m", "120h", raw_time))
      }
    }
  }
  
  metadata$Time_Point[metadata$Time_Point == "NA" | is.na(metadata$Time_Point)] <- "0min"
  
  count_file <- paste0("/data2t_2/hushu/02.DESeq2_analysis_TE_loci/raw_data/", 
                       dataset_name, "_readscounts_matrix_TE_Loci.txt")
  
  if (!file.exists(count_file)) {
    warning(paste("Count file not found for dataset:", dataset_name))
    return(NULL)
  }
  
  count <- fread(count_file)
  count <- count[, c("Geneid", metadata$Sample_Name), with = FALSE]
  count <- as.data.frame(count)
  rownames(count) <- count$Geneid
  count_TE <- count[, -1, drop = FALSE] 
  count2_TE <- as.matrix(count_TE)
  count2_TE <- count2_TE[, rownames(metadata)]
  
  dds <- DESeqDataSetFromMatrix(countData = count2_TE, colData = metadata, design = ~ Condition)
  vsd <- vst(dds, blind = TRUE)
  normalized_counts_rmBE_TE <- as.data.frame(assay(vsd))
  exp_df_TE <- normalized_counts_rmBE_TE[rownames(normalized_counts_rmBE_TE) %in% Top_5_Up_TE_loci$TE_id, ]
  
  metadata <- metadata %>%
    mutate(order = match(rownames(metadata), colnames(exp_df_TE))) %>%
    arrange(order) %>%
    select(-order)
  
  TE_normalized_df <- as.data.frame(exp_df_TE) %>%
    tibble::rownames_to_column("TE_id") %>%
    pivot_longer(
      cols = -TE_id,
      names_to = "Sample_Name",
      values_to = "normalized_count"
    )
  
  TE_normalized_df <- TE_normalized_df %>%
    left_join(metadata %>% 
                tibble::rownames_to_column("Sample_ID") %>%
                dplyr::select(Sample_Name, Time_Point, Infection_State),
              by = "Sample_Name")
  
  TE_normalized_df$Time_Point_Hours <- sapply(TE_normalized_df$Time_Point, convert_time_to_hours)
  
  TE_normalized_df$Dataset <- dataset_name
  
  return(TE_normalized_df)
}

all_normalized_data <- list()

for (dataset in datasets_to_analyze) {
  cat("Processing dataset:", dataset, "\n")
  
  tryCatch({
    result <- process_single_dataset(dataset, all_metadata, Top_5_Up_TE_loci)
    if (!is.null(result)) {
      all_normalized_data[[dataset]] <- result
    }
  }, error = function(e) {
    cat("Error processing dataset", dataset, ":", conditionMessage(e), "\n")
  })
}

combined_data <- do.call(rbind, all_normalized_data)
combined_data$x_position <- NA
combined_data$Dataset_Time <- paste(combined_data$Dataset, combined_data$Time_Point_Hours, sep = "_")

dataset_positions <- list()
current_pos <- 1
x_positions <- c()
x_labels <- c()
dataset_boundaries <- c()

for (i in seq_along(datasets_to_analyze)) {
  dataset <- datasets_to_analyze[i]
  if (dataset %in% names(all_normalized_data)) {
    dataset_data <- all_normalized_data[[dataset]]
    if (dataset == "Virus_HCV_1") {
      dataset_data$Time_Point_Hours[dataset_data$Time_Point_Hours == "4320min"] <- "72h"
      dataset_data$Time_Point_Hours[dataset_data$Time_Point_Hours == "7200min"] <- "120h"
      
      dataset_data$time_infection_combo <- paste(dataset_data$Time_Point_Hours, 
                                                 dataset_data$Infection_State, sep = "_")
      
      unique_combos <- unique(dataset_data$time_infection_combo)
      combo_order <- c()
      if ("72h_Pre" %in% unique_combos) combo_order <- c(combo_order, "72h_Pre")
      if ("72h_Post" %in% unique_combos) combo_order <- c(combo_order, "72h_Post")
      if ("120h_Post" %in% unique_combos) combo_order <- c(combo_order, "120h_Post")
      
      dataset_positions[[dataset]] <- current_pos:(current_pos + length(combo_order) - 1)
      
      for (combo in combo_order) {
        x_positions <- c(x_positions, current_pos)
        time_label <- strsplit(combo, "_")[[1]][1]
        x_labels <- c(x_labels, time_label)
        
        time_point <- strsplit(combo, "_")[[1]][1]
        infection_state <- strsplit(combo, "_")[[1]][2]
        
        combined_data$x_position[combined_data$Dataset == dataset & 
                                   combined_data$Time_Point_Hours == time_point &
                                   combined_data$Infection_State == infection_state] <- current_pos
        
        current_pos <- current_pos + 1
      }
      
    } else {
      unique_times <- unique(dataset_data$Time_Point_Hours)
      
      non_zero_times <- unique_times[unique_times != "0h"]
      if (length(non_zero_times) > 0) {
        time_numbers <- as.numeric(gsub("[^0-9.]", "", non_zero_times))
        sorted_indices <- order(time_numbers)
        sorted_times <- non_zero_times[sorted_indices]
        time_order <- c("0h", sorted_times)
      } else {
        time_order <- "0h"
      }
      
      time_order <- time_order[time_order %in% unique_times]
      
      dataset_positions[[dataset]] <- current_pos:(current_pos + length(time_order) - 1)
      
      for (time_point in time_order) {
        x_positions <- c(x_positions, current_pos)
        x_labels <- c(x_labels, time_point)
        
        combined_data$x_position[combined_data$Dataset == dataset & 
                                   combined_data$Time_Point_Hours == time_point] <- current_pos
        
        current_pos <- current_pos + 1
      }
    }
    
    if (i < length(datasets_to_analyze)) {
      dataset_boundaries <- c(dataset_boundaries, current_pos - 0.5)
    }
  }
}

dataset_label_positions <- sapply(dataset_positions, function(pos) mean(pos))

for (dataset in names(all_normalized_data)) {
  dataset_data <- all_normalized_data[[dataset]]
  for (time_point in unique(dataset_data$Time_Point_Hours)) {
    for (infection_state in unique(dataset_data$Infection_State)) {
      subset_data <- dataset_data[dataset_data$Time_Point_Hours == time_point & 
                                    dataset_data$Infection_State == infection_state, ]
      if (nrow(subset_data) > 0) {
        values <- subset_data$normalized_count
        cat(sprintf("%s - %s - %s: n=%d, median=%.2f, range=[%.2f, %.2f], unique_values=%d\n", 
                    dataset, time_point, infection_state, 
                    length(values), median(values), min(values), max(values), 
                    length(unique(values))))
      }
    }
  }
}


for (dataset in names(all_normalized_data)) {
  dataset_data <- all_normalized_data[[dataset]]
  for (time_point in unique(dataset_data$Time_Point_Hours)) {
    for (infection_state in unique(dataset_data$Infection_State)) {
      subset_data <- dataset_data[dataset_data$Time_Point_Hours == time_point & 
                                    dataset_data$Infection_State == infection_state, ]
      if (nrow(subset_data) > 0) {
        values <- subset_data$normalized_count + 1
        q <- quantile(values, c(0, 0.25, 0.5, 0.75, 1))
        cat(sprintf("%s-%s-%s: Min=%.2f, Q1=%.2f, Median=%.2f, Q3=%.2f, Max=%.2f\n", 
                    dataset, time_point, infection_state, q[1], q[2], q[3], q[4], q[5]))
      }
    }
  }
}

p_linear_full <- ggplot(combined_data, aes(x = factor(x_position), y = normalized_count, fill = Infection_State, color = Infection_State)) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.05,
    outlier.alpha = 0.5,
    alpha = 0,  
    size = 0.3,   
    notch = FALSE,
    fatten = 1,  
    coef = 1.5   
  ) +
  geom_vline(xintercept = dataset_boundaries, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.5) +
  scale_y_continuous(
    trans = "log10",  
    breaks = c(1, 2, 5, 10, 20, 50, 100),
    labels = c("1", "2", "5", "10", "20", "50", "100"),  
    expand = expansion(mult = c(0.05, 0.05)) 
  ) +
  scale_x_discrete(breaks = as.character(x_positions), labels = x_labels) +
  labs(
    title = "TE Loci Expression", 
    x = "Time After Infection", 
    y = "log10 Normalized Count", 
    fill = "Condition",
    color = "Condition"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("Pre" = "#4c9bcf", "Post" = "#d00732")) +
  scale_color_manual(values = c("Pre" = "#4c9bcf", "Post" = "#d00732")) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 7, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    axis.line = element_line(color = "black", size = 0.3),
    legend.box = "horizontal"
  )

min_y <- min(combined_data$normalized_count, na.rm = TRUE)
for (dataset in names(dataset_label_positions)) {
  p_linear_full <- p_linear_full + 
    annotate("text", 
             x = dataset_label_positions[[dataset]], 
             y = min_y * 0.8+17,  
             label = dataset, 
             size = 1.5, 
             hjust = 0.5,
             color = "black",
             fontface = "bold")
}

print(p_linear_full)
ggsave("diff_Time_Point/Figure4_Recurrent_TE_Loci_Time_point_expression_plot.pdf", p_linear_full, width = 5, height = 2)

