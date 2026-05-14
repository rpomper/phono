

# load --------------------------------------------------------------------


setwd('/Users/pomp/Library/CloudStorage/OneDrive-MichiganStateUniversity/Code/phono web/prep code')
d <- read.csv('hml with homophones.csv')
a <- read.csv('AoA.csv')



# length, late8, aoa ------------------------------------------------------


library(dplyr)
library(stringr)

late8_phonemes <- c("S","Z","T","D","s","z","l","r")

# Build the complete raw_lexicon in one pipeline
raw_lexicon <- d %>%
  # 1. Merge the AoA/Frequency data onto the HML data by spelling
  left_join(a, by = "ortho") %>%
  
  # 1.5 FIX: Force frequency and AoA to be numeric
  mutate(
    Freq_pm = as.numeric(Freq_pm),
    Rating.Mean = as.numeric(Rating.Mean)
  ) %>%
  
  # 2. Group by the phonological form
  group_by(klatt) %>%
  
  # 3. Aggregate
  summarise(
    ortho = paste(unique(ortho), collapse = "; "),
    
    # If all values are NA, return NA. Otherwise, sum them.
    freq = if(all(is.na(Freq_pm))) NA_real_ else sum(Freq_pm, na.rm = TRUE),
    
    aoa = suppressWarnings(min(Rating.Mean, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  
  # 4. Add the static phonological variables & clean up Inf
  mutate(
    length = nchar(klatt),
    # FIX: Use str_count to count the actual number of late-8 phonemes per word
    late8 = as.integer(str_count(klatt, paste0("[", paste(late8_phonemes, collapse = ""), "]"))),
    aoa = na_if(aoa, Inf)
  )

# Check the results
head(raw_lexicon)

# full --------------------------------------------------------------------


library(igraph)
library(stringdist)

# Define the ages you want to loop through (e.g., 3 to 26)
ages <- 3:26
# ages <- 3:10

# Initialize two empty lists to hold just our lightweight summaries
master_norms_list <- list()
master_bins_list <- list()

for (current_age in ages) {
  
  # 1. Filter the master lexicon for this specific age
  # (Keeping words acquired at or before current_age)
  age_lexicon <- raw_lexicon %>%
    filter(aoa <= current_age)
  
  # Extract the unique phonological forms
  klatt_words <- age_lexicon$klatt
  
  # this is the original process, which breaks around age 19
  # because the lexicon becomes too large
  # and so the levenshtein distance comparisons exceeds the vector memory limit (16.0 gb)
  
  # # 2. Calculate the string distance matrix (Levenshtein distance)
  # # This finds the number of additions, deletions, or substitutions 
  # # needed to turn one word into another.
  # dist_matrix <- stringdistmatrix(klatt_words, klatt_words, method = "lv")
  # 
  # # 3. Create an adjacency matrix (1 if distance is exactly 1, else 0)
  # adj_matrix <- (dist_matrix == 1) * 1
  # rownames(adj_matrix) <- klatt_words
  # colnames(adj_matrix) <- klatt_words
  
  # Create an empty list to store our edges
  edge_list_chunks <- list()
  
  # Find all the unique word lengths in our lexicon (e.g., 1 to 15 phonemes)
  lengths <- sort(unique(nchar(klatt_words)))
  
  # 2. Loop through each word length
  for (L in lengths) {
    
    words_L <- klatt_words[nchar(klatt_words) == L]
    
    # A: Compare words of length L against OTHER words of length L (Substitutions)
    if (length(words_L) > 1) {
      dist_mat <- stringdistmatrix(words_L, words_L, method = "lv")
      edges <- which(dist_mat == 1, arr.ind = TRUE)
      edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE] # Remove duplicates
      
      if (nrow(edges) > 0) {
        edge_list_chunks[[paste0("sub_", L)]] <- data.frame(
          from = words_L[edges[, 1]],
          to = words_L[edges[, 2]],
          stringsAsFactors = FALSE
        )
      }
    }
    
    # B: Compare words of length L against words of length L + 1 (Insertions/Deletions)
    # (We don't need L-1 because the previous loop iteration already caught it!)
    words_L_plus_1 <- klatt_words[nchar(klatt_words) == (L + 1)]
    
    if (length(words_L) > 0 && length(words_L_plus_1) > 0) {
      dist_mat <- stringdistmatrix(words_L, words_L_plus_1, method = "lv")
      edges <- which(dist_mat == 1, arr.ind = TRUE)
      
      if (nrow(edges) > 0) {
        edge_list_chunks[[paste0("ins_", L)]] <- data.frame(
          from = words_L[edges[, 1]],
          to = words_L_plus_1[edges[, 2]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # 3. Combine all chunks into one master edge list
  edge_list <- bind_rows(edge_list_chunks)
  
  # 4. Build the undirected network graph
  # g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = data.frame(name = klatt_words))
  
  # 5. Calculate the network statistics
  # Degree (number of edges)
  node_degrees <- degree(g)
  
  # Clustering Coefficient (local transitivity)
  # igraph natively returns NaN for nodes with degree < 2. 
  # 'isolates = "zero"' forces it to return 0 instead of NA/NaN.
  # Let igraph return NaN for nodes with degree < 2, then convert to NA
  node_c <- transitivity(g, type = "local")
  node_c[is.nan(node_c)] <- NA_real_
  
  # Neighbors (a semicolon-separated string of the neighboring klatt strings)
  # adj_list <- as_adj_list(g)
  # node_neighbors <- sapply(adj_list, function(x) paste(names(x), collapse = ";"))
  
  # 6. Store in a dataframe for this specific age
  age_stats <- data.frame(
    ortho = age_lexicon$ortho,
    klatt = klatt_words,
    length = age_lexicon$length,
    late8 = age_lexicon$late8,
    freq = as.numeric(age_lexicon$freq), # Ensuring numeric for math
    aoa = age_lexicon$aoa, 
    degree = as.integer(node_degrees),
    c = round(node_c, 3),             
    stringsAsFactors = FALSE
  )
  
  # 7. Calculate aggregate statistics for current age vocabulary
  
  # 7a. Calculate Norms (Averages & SDs) for this age
  master_norms_list[[as.character(current_age)]] <- age_stats %>%
    summarise(
      age = current_age,
      phonemes_avg = round(mean(length, na.rm = TRUE), 2),
      phonemes_sd  = round(sd(length, na.rm = TRUE), 2),
      late8_avg    = round(mean(late8, na.rm = TRUE), 2),
      late8_sd     = round(sd(late8, na.rm = TRUE), 2),
      freq_avg     = round(mean(freq, na.rm = TRUE), 2),
      freq_sd      = round(sd(freq, na.rm = TRUE), 2),
      freq_median  = round(median(freq, na.rm = TRUE), 2),
      aoa_avg      = round(mean(aoa, na.rm = TRUE), 2),
      aoa_sd       = round(sd(aoa, na.rm = TRUE), 2),
      degree_avg   = round(mean(degree, na.rm = TRUE), 2),
      degree_sd    = round(sd(degree, na.rm = TRUE), 2),
      c_avg        = round(mean(c, na.rm = TRUE), 2),
      c_sd         = round(sd(c, na.rm = TRUE), 2)
    )
  
  # 7b. Calculate Histogram Bin Counts directly matching your app's logic
  aoa_bins <- age_stats %>%
    mutate(val = round(aoa, 0)) %>%
    count(val, name = "n") %>%
    mutate(metric = "aoa", age = current_age)
  
  phoneme_bins <- age_stats %>%
    count(length, name = "n") %>%
    rename(val = length) %>%
    mutate(metric = "phonemes", age = current_age)
  
  late8_bins <- age_stats %>%
    count(late8, name = "n") %>%
    rename(val = late8) %>%
    mutate(metric = "late8", age = current_age)
  
  degree_bins <- age_stats %>%
    count(degree, name = "n") %>%
    rename(val = degree) %>%
    mutate(metric = "degree", age = current_age)
  
  c_bins <- age_stats %>%
    filter(!is.na(c)) %>%
    mutate(val = as.numeric(as.character(cut(c, seq(-.05, 1.0, .05), seq(0, 1, .05))))) %>%
    count(val, name = "n") %>%
    mutate(metric = "c", age = current_age)
  
  # Stack all the bins for this age together
  master_bins_list[[as.character(current_age)]] <- bind_rows(
    aoa_bins, phoneme_bins, late8_bins, degree_bins, c_bins
  )
  
  cat("Finished aggregating age:", current_age, "| Words processed:", length(klatt_words), "\n")
  
  }

# # 8. Combine into our two final, lightweight tables at the very end
# final_norms <- bind_rows(master_norms_list)
# final_bins  <- bind_rows(master_bins_list)


# # 9. EXPORT FILES (CSV for humans, RDS for the Shiny app)
# 
# setwd('/Users/rpomper/Library/CloudStorage/OneDrive-MichiganStateUniversity/Code/phono web/data')
# 
# # Save CSVs for independent investigation and sharing
# write.csv(final_norms, "child_norms_summary.csv", row.names = FALSE)
# write.csv(final_bins,  "child_hist_bins.csv",     row.names = FALSE)
# 
# # Save RDS files for lightning-fast loading inside the Shiny app
# saveRDS(final_norms, "child_norms_summary.rds")
# saveRDS(final_bins,  "child_hist_bins.rds")
# 
# cat("Success! Exported both CSV and RDS formats.\n")



# add in norms for the full lexicon ---------------------------------------

# i.e., not filtering to only include words with AoA data


current_age <- "adult"
# 1. Filter the master lexicon for this specific age
# (Keeping words acquired at or before current_age)
age_lexicon <- raw_lexicon

# Extract the unique phonological forms
klatt_words <- age_lexicon$klatt

# this is the original process, which breaks around age 19
# because the lexicon becomes too large
# and so the levenshtein distance comparisons exceeds the vector memory limit (16.0 gb)

# # 2. Calculate the string distance matrix (Levenshtein distance)
# # This finds the number of additions, deletions, or substitutions 
# # needed to turn one word into another.
# dist_matrix <- stringdistmatrix(klatt_words, klatt_words, method = "lv")
# 
# # 3. Create an adjacency matrix (1 if distance is exactly 1, else 0)
# adj_matrix <- (dist_matrix == 1) * 1
# rownames(adj_matrix) <- klatt_words
# colnames(adj_matrix) <- klatt_words

# Create an empty list to store our edges
edge_list_chunks <- list()

# Find all the unique word lengths in our lexicon (e.g., 1 to 15 phonemes)
lengths <- sort(unique(nchar(klatt_words)))

# 2. Loop through each word length
for (L in lengths) {
  
  words_L <- klatt_words[nchar(klatt_words) == L]
  
  # A: Compare words of length L against OTHER words of length L (Substitutions)
  if (length(words_L) > 1) {
    dist_mat <- stringdistmatrix(words_L, words_L, method = "lv")
    edges <- which(dist_mat == 1, arr.ind = TRUE)
    edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE] # Remove duplicates
    
    if (nrow(edges) > 0) {
      edge_list_chunks[[paste0("sub_", L)]] <- data.frame(
        from = words_L[edges[, 1]],
        to = words_L[edges[, 2]],
        stringsAsFactors = FALSE
      )
    }
  }
  
  # B: Compare words of length L against words of length L + 1 (Insertions/Deletions)
  # (We don't need L-1 because the previous loop iteration already caught it!)
  words_L_plus_1 <- klatt_words[nchar(klatt_words) == (L + 1)]
  
  if (length(words_L) > 0 && length(words_L_plus_1) > 0) {
    dist_mat <- stringdistmatrix(words_L, words_L_plus_1, method = "lv")
    edges <- which(dist_mat == 1, arr.ind = TRUE)
    
    if (nrow(edges) > 0) {
      edge_list_chunks[[paste0("ins_", L)]] <- data.frame(
        from = words_L[edges[, 1]],
        to = words_L_plus_1[edges[, 2]],
        stringsAsFactors = FALSE
      )
    }
  }
}

# 3. Combine all chunks into one master edge list
edge_list <- bind_rows(edge_list_chunks)

# 4. Build the undirected network graph
# g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = data.frame(name = klatt_words))

# 5. Calculate the network statistics
# Degree (number of edges)
node_degrees <- degree(g)

# Clustering Coefficient (local transitivity)
# igraph natively returns NaN for nodes with degree < 2. 
# 'isolates = "zero"' forces it to return 0 instead of NA/NaN.
# Let igraph return NaN for nodes with degree < 2, then convert to NA
node_c <- transitivity(g, type = "local")
node_c[is.nan(node_c)] <- NA_real_

# Neighbors (a semicolon-separated string of the neighboring klatt strings)
# adj_list <- as_adj_list(g)
# node_neighbors <- sapply(adj_list, function(x) paste(names(x), collapse = ";"))

# 6. Store in a dataframe for this specific age
age_stats <- data.frame(
  ortho = age_lexicon$ortho,
  klatt = klatt_words,
  length = age_lexicon$length,
  late8 = age_lexicon$late8,
  freq = as.numeric(age_lexicon$freq), # Ensuring numeric for math
  aoa = age_lexicon$aoa, 
  degree = as.integer(node_degrees),
  c = round(node_c, 3),             
  stringsAsFactors = FALSE
)

# 7. Calculate aggregate statistics for current age vocabulary

# 7a. Calculate Norms (Averages & SDs) for this age
# Adult Norms
master_norms_list[[as.character(current_age)]] <- age_stats %>%
  summarise(
    age = current_age,
    phonemes_avg = round(mean(length, na.rm=TRUE), 2), phonemes_sd  = round(sd(length, na.rm=TRUE), 2),
    late8_avg    = round(mean(late8, na.rm=TRUE), 2),  late8_sd     = round(sd(late8, na.rm=TRUE), 2),
    freq_avg     = round(mean(freq, na.rm=TRUE), 2),   freq_sd      = round(sd(freq, na.rm=TRUE), 2), freq_median  = round(median(freq, na.rm=TRUE), 2),
    aoa_avg      = round(mean(aoa, na.rm=TRUE), 2),    aoa_sd       = round(sd(aoa, na.rm=TRUE), 2),
    degree_avg   = round(mean(degree, na.rm=TRUE), 2), degree_sd    = round(sd(degree, na.rm=TRUE), 2),
    c_avg        = round(mean(c, na.rm=TRUE), 2),      c_sd         = round(sd(c, na.rm=TRUE), 2)
  )

# Adult Bins
aoa_bins     <- age_stats %>% filter(!is.na(aoa)) %>% mutate(val = round(aoa, 0)) %>% count(val, name = "n") %>% mutate(metric = "aoa", age = current_age)
phoneme_bins <- age_stats %>% count(length, name = "n") %>% rename(val = length) %>% mutate(metric = "phonemes", age = current_age)
late8_bins   <- age_stats %>% count(late8, name = "n") %>% rename(val = late8) %>% mutate(metric = "late8", age = current_age)
degree_bins  <- age_stats %>% count(degree, name = "n") %>% rename(val = degree) %>% mutate(metric = "degree", age = current_age)
c_bins       <- age_stats %>% filter(!is.na(c)) %>% mutate(val = as.numeric(as.character(cut(c, seq(-.05, 1.0, .05), seq(0, 1, .05))))) %>% count(val, name = "n") %>% mutate(metric = "c", age = current_age)

master_bins_list[[as.character(current_age)]] <- bind_rows(aoa_bins, phoneme_bins, late8_bins, degree_bins, c_bins)
cat("Finished adult lexicon! Combining and exporting tables...\n")

# 8. Combine into our two final, lightweight tables at the very end
# final_norms <- bind_rows(master_norms_list)
# final_bins  <- bind_rows(master_bins_list)



# export ------------------------------------------------------------------

# EXPORT FILES (CSV for humans, RDS for the Shiny app)

# 1. Force the 'age' column in every dataframe inside the lists to be a character string
master_norms_list <- lapply(master_norms_list, function(df) {
  df$age <- as.character(df$age)
  return(df)
})

master_bins_list <- lapply(master_bins_list, function(df) {
  df$age <- as.character(df$age)
  return(df)
})

# 2. Now combine and export safely!
final_norms <- bind_rows(master_norms_list)
final_bins  <- bind_rows(master_bins_list)

setwd('/Users/pomp/Library/CloudStorage/OneDrive-MichiganStateUniversity/Code/phono web/data')

# Save CSVs for independent investigation and sharing
write.csv(final_norms, "norms_summary.csv", row.names = FALSE)
write.csv(final_bins,  "hist_bins.csv",     row.names = FALSE)

# Save RDS files for lightning-fast loading inside the Shiny app
saveRDS(final_norms, "norms_summary.rds")
saveRDS(final_bins,  "hist_bins.rds")

cat("Success! Exported both CSV and RDS formats.\n")

write.csv(raw_lexicon,"full_lexicon.csv")
saveRDS(raw_lexicon,"full_lexicon.rds")

ortho_lookup <- d %>% select(ortho, klatt)

# Save it to your app's data folder
saveRDS(ortho_lookup, "ortho_lookup.rds")

# degree, c ---------------------------------------------------------------

library(igraph)
library(stringdist)
library(dplyr)

# 1. Filter for Age 26
age_lexicon <- raw_lexicon %>% filter(aoa <= 4)
klatt_words <- age_lexicon$klatt

# Create an empty list to store our edges
edge_list_chunks <- list()

# Find all the unique word lengths in our lexicon (e.g., 1 to 15 phonemes)
lengths <- sort(unique(nchar(klatt_words)))

# 2. Loop through each word length
for (L in lengths) {
  
  words_L <- klatt_words[nchar(klatt_words) == L]
  
  # A: Compare words of length L against OTHER words of length L (Substitutions)
  if (length(words_L) > 1) {
    dist_mat <- stringdistmatrix(words_L, words_L, method = "lv")
    edges <- which(dist_mat == 1, arr.ind = TRUE)
    edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE] # Remove duplicates
    
    if (nrow(edges) > 0) {
      edge_list_chunks[[paste0("sub_", L)]] <- data.frame(
        from = words_L[edges[, 1]],
        to = words_L[edges[, 2]],
        stringsAsFactors = FALSE
      )
    }
  }
  
  # B: Compare words of length L against words of length L + 1 (Insertions/Deletions)
  # (We don't need L-1 because the previous loop iteration already caught it!)
  words_L_plus_1 <- klatt_words[nchar(klatt_words) == (L + 1)]
  
  if (length(words_L) > 0 && length(words_L_plus_1) > 0) {
    dist_mat <- stringdistmatrix(words_L, words_L_plus_1, method = "lv")
    edges <- which(dist_mat == 1, arr.ind = TRUE)
    
    if (nrow(edges) > 0) {
      edge_list_chunks[[paste0("ins_", L)]] <- data.frame(
        from = words_L[edges[, 1]],
        to = words_L_plus_1[edges[, 2]],
        stringsAsFactors = FALSE
      )
    }
  }
}

# 3. Combine all chunks into one master edge list
edge_list <- bind_rows(edge_list_chunks)

# 4. Build the undirected network graph directly from the edges
# We include all 'klatt_words' as vertices so words with 0 neighbors are still included
g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = data.frame(name = klatt_words))

# 5. Calculate the network statistics
node_degrees <- degree(g)
node_c <- transitivity(g, type = "local", isolates = "zero")
adj_list <- as_adj_list(g)
node_neighbors <- sapply(adj_list, function(x) paste(names(x), collapse = ";"))

# 6. Store in a dataframe
age_stats <- data.frame(
  klatt = klatt_words,
  age = 26,
  degree = as.integer(node_degrees),
  c = round(node_c, 3),
  neighbors = node_neighbors,
  stringsAsFactors = FALSE
)

# Check the results!
head(age_stats)

