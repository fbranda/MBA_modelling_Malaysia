### Load library ----

#library everywhere needed
library(tidyverse) #manipulate datasets
library(readr) #load the data

#library factor analysis
library(corrplot) #corrplot
library(psych) #parallel

#library DPMM model section
library(mvtnorm)
library(clusterSim) #cluster models validation
library(fpc) #cluster.stat
library(viridis)
library(gridExtra)

#library for maps and spatial distribution sections
library(sf)
library(ggspatial)  
library(INLA) #spatial
library(sp)
library(glue) #manipulate characters
library(inlabru) #spatial

setwd("")

### Load the data ----
data_MBA_sim <- read_csv("MBA modelling/Input Data/data_MBA_modelling.csv")

### Factor analysis ----

#select the mfi 
MFI<- data_MBA_sim %>% dplyr::select(2:13)

#summaries
#summaries<-malaysiaMFI|> map_df(summary)
#rownames(summaries) <- colnames(malaysiaMFI)

#create the log
log_mfi <- log(MFI+1)

#summaries log
#summaries_log<-log_mfi|> map_df(summary)
#rownames(summaries_log) <- colnames(log_mfi)

#create the scale log
slog_mfi<-scale(log_mfi)

#summaries slog
summaries_slog<-summary(slog_mfi)

##correlation matrix
corrplot(cor(slog_mfi),method = "color",addCoef.col = "black", title="Correlation matrix scaled log-MFI values",mar=c(0,0,1,0), type = "upper")

###### giardia ----
KMO(slog_mfi[,2:3])
cortest.bartlett(slog_mfi[,2:3]) #barlett test

#decide the number of factors:
fa.parallel(slog_mfi[,2:3], fa = "fa", n.iter = 100) #1
eigen(cor(slog_mfi[,2:3]))$values #1

fa_giardia <- fa(slog_mfi[,2:3], nfactors = 1, rotate = "none")

giardia <- fa_giardia$scores

###### lf ----
KMO(slog_mfi[,4:7])
cortest.bartlett(slog_mfi[,4:7]) #barlett test

#decide the number of factors:
fa.parallel(slog_mfi[,4:7], fa = "fa", n.iter = 100) #1
eigen(cor(slog_mfi[,4:7]))$values #1

fa_lf <- fa(slog_mfi[,4:7], nfactors = 1, rotate = "none")

lf <- fa_lf$scores

###### trachoma ----
KMO(slog_mfi[,8:9])
cortest.bartlett(slog_mfi[,8:9]) #barlett test

#decide the number of factors:
fa.parallel(slog_mfi[,8:9], fa = "fa", n.iter = 100) #1
eigen(cor(slog_mfi[,8:9]))$values #1

fa_trachoma <- fa(slog_mfi[,8:9], nfactors = 1, rotate = "none")

trachoma <- fa_trachoma$scores

###### yaws ----
KMO(slog_mfi[,11:12])
cortest.bartlett(slog_mfi[,11:12]) #barlett test

#decide the number of factors:
fa.parallel(slog_mfi[,11:12], fa = "fa", n.iter = 100) #1
eigen(cor(slog_mfi[,11:12]))$values #1

fa_yaws <- fa(slog_mfi[,11:12], nfactors = 1, rotate = "none")

yaws<- fa_yaws$scores

##### DB with diseases ----

toxo<- slog_mfi[,1]
strongy<- slog_mfi[,10]

dismfi_slog<- cbind(toxo, giardia, lf, trachoma, strongy, yaws)

dismfi_slog<- as.data.frame(dismfi_slog)
colnames(dismfi_slog)<-c("Toxoplasmosis", "Giardiasis", "Lymphatic Filariasis", "Trachoma", "Strongyloidiasis", "Yaws")

dismfi_slog_c<- dismfi_slog
dismfi_slog_c$subjects<-seq(1,dim(dismfi_slog)[1],1)

dismfi_slog_c<-as.data.frame(dismfi_slog_c) %>%
  pivot_longer("Toxoplasmosis":"Yaws", names_to="antigen", values_to="slogMFI") %>% dplyr::select(antigen,slogMFI,subjects)

ggplot(dismfi_slog_c, aes(slogMFI)) + 
  geom_density() +
  facet_wrap(~ antigen, ncol=2)+
  theme_minimal()

## merge the databases
data_MBA_modelling<-cbind(data_MBA_sim, dismfi_slog)

### DPMM models with different combinations of alpha and sigma_y  ----
set.seed(42)

alpha_vals <- c(0.01, 1.00, 2.00, 3.00)
sigma_y_vals <- c(0.50, 1.00, 1.5, 3.00)

c_init <- rep(1, nrow(dismfi_slog))
mu0 <- matrix(rep(0, 6), ncol = 6, byrow = TRUE)
sigma0 <- diag(6) * 1

# List to save the results locally
all_results <- list()

### PLEASE RUN THE "crp_gibbs_fun.R" FILE!!!!!! ###

for (i in seq_along(alpha_vals)) {
  for (j in seq_along(sigma_y_vals)) {
    
    alpha_val <- alpha_vals[i]
    sigma_val <- sigma_y_vals[j]
    
    # Esegui il modello
    results <- crp_gibbs(
      data = dismfi_slog,
      alpha = alpha_val,
      mu0 = mu0,
      sigma0 = sigma0,
      sigma_y = diag(6) * sigma_val,
      c_init = c_init,
      maxIters = 1000
    )
    
    # Salva i risultati in lista (opzionale)
    name <- paste0("alpha_", alpha_val, "_sigma_", sigma_val)
    all_results[[name]] <- results
    
    # Salva il file sul disco
    saveRDS(results, file = paste0("results_", name, ".rds"))
    
    # Calcola il cluster più frequente per ogni osservazione
    tab <- apply(results, 1, function(x) {
      cluster_counts <- table(x)
      most_frequent <- names(cluster_counts[which.max(cluster_counts)])
      return(most_frequent)
    })
    
    # Output su console
    cat("alpha = ", alpha_val, "sigma_y = ", sigma_val, "\n")
    print(table(tab))
  }
}

#### sensitivity analysis ----
var_names <- ls(pattern = "^results_alpha_[0-9.]+_sigma_[0-9.]+$")

# Create a named list from the variables
results_list <- setNames(lapply(var_names, get), var_names)

save(results_list,file = "alphasigma_listM.RData")

###### load the saved sensitivity analysis results ----
load("alphasigma_listM.RData")

# Apply your function to each
mode_function <- function(mat) {
  apply(mat, 1, function(x) {
    freq <- table(x)
    as.numeric(names(freq)[which.max(freq)])
  })
}

processed_results <- lapply(results_list, mode_function)

#number of clusters and units by clusters for each combination
purrr::map(processed_results,table)

## compute metrics
# Compute distance matrix once
d <- dist(dismfi_slog)

# Function to compute metrics
compute_four_indices  <- function(d, clustering, data) {
  if (length(unique(clustering)) < 2) {
    return(c(Silhouette = NA, Dunn = NA, DB = NA, CH = NA))
  }
  sil <- index.S(d, clustering, singleObject = 0)
  dunn <- cluster.stats(d, clustering)$dunn
  db <- index.DB(data, clustering)$DB
  ch <- cluster.stats(d, clustering)$ch
  return(c(Silhouette = sil, Dunn = dunn, DB = db, CH = ch))
}

# list to save the results in
summary_list <- vector("list", length(results_list))
names(summary_list) <- names(results_list)

for (nm in names(results_list)) {
  mat_res <- results_list[[nm]] 
  K_per_iter <- apply(mat_res, 2, function(col) length(unique(col)))
  K_mean     <- mean(K_per_iter)
  K_sd       <- sd(K_per_iter)
  
  part_map <- mode_function(mat_res)   
  
  idx4 <- compute_four_indices(d, part_map, as.matrix(dismfi_slog))

  summary_list[[nm]] <- c(
    alpha       = as.numeric(sub(".*alpha_([0-9.]+)_sigma_.*", "\\1", nm)),
    sigma       = as.numeric(sub(".*sigma_([0-9.]+)$",        "\\1", nm)),
    K_mean      = K_mean,
    K_sd        = K_sd,
    Silhouette  = idx4["Silhouette"],
    Dunn        = idx4["Dunn"],
    DB          = idx4["DB"],
    CH          = idx4["CH"]
  )
}

df_summary <- bind_rows(
  lapply(names(summary_list), function(nm) {
    as.data.frame(as.list(summary_list[[nm]]))
  }),
  .id = "comb"
)

df_summary <- df_summary %>%
  rename(
    Silhouette = `Silhouette.Silhouette`,
    Dunn       = `Dunn.Dunn`,
    DB         = `DB.DB`,
    CH         = `CH.CH`
  ) %>%
  mutate(
    alpha       = as.numeric(alpha),
    sigma       = as.numeric(sigma),
    K_mean      = as.numeric(K_mean),
    K_sd        = as.numeric(K_sd),
    Silhouette  = as.numeric(Silhouette),
    Dunn        = as.numeric(Dunn),
    DB          = as.numeric(DB),
    CH          = as.numeric(CH)
  )

df_summary_ord <- df_summary %>%
  arrange(desc(Silhouette))

print(df_summary_ord)

###### clusters plots ----
malaysia_dpmm<-data_MBA_modelling

#the following line depend on the best combination given by the sensitivity analysis.
malaysia_dpmm$cluster<-processed_results[["results_alpha_3_sigma_1.5"]]

malaysia_dpmm<-malaysia_dpmm|> mutate(
  occupation_grouped=as.factor(occupation),
  cluster=as.factor(cluster),
  stratum=as.factor(veg.stratum),
  wealth=as.factor(wealth),
  gender=as.factor(gender)
)

# Convert to long format
malaysia_long <- malaysia_dpmm %>%
  pivot_longer(
    cols = c(Toxoplasmosis, Giardiasis, `Lymphatic Filariasis`,
             Trachoma, Strongyloidiasis, Yaws),
    names_to = "disease",
    values_to = "value"
  )

# Summarise means and 95% CI
summary_df <- malaysia_long %>%
  group_by(cluster, disease) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = if_else(n > 0, sd / sqrt(n), NA_real_),
    ci = 1.96 * se,
    lower = mean - ci,
    upper = mean + ci,
    .groups = "drop"
  )

# Position dodge for spacing between clusters
pd <- position_dodge(width = 0.6)

# Plot with shape 21 and fill color by cluster
p0<-ggplot(summary_df, aes(x = disease, y = mean, fill = cluster)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, color = cluster),
    width = 0.4,
    position = pd,
    size = 0.8
  ) +
  geom_point(
    shape = 21, size = 4, stroke = 1,
    position = pd, color = "black"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "",
    x = "Disease",
    y = "Mean scaled log-MFI",
    fill = "Cluster",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    axis.text = element_text(size = 12)
  )
p0

#####continuous
# Position dodge so points / CIs don't overlap
pd <- position_dodge(width = 0.6)

# Summary with mean and 95% CI for age and dem
summary_df <- malaysia_dpmm %>%
  pivot_longer(cols = c(age, pop.density), names_to = "variable", values_to = "value") %>%
  group_by(cluster, variable) %>%
  summarise(
    mean  = mean(value, na.rm = TRUE),
    sd    = sd(value, na.rm = TRUE),
    n     = sum(!is.na(value)),
    se    = if_else(n > 0, sd / sqrt(n), NA_real_),
    ci    = 1.96 * se,
    lower = mean - ci,
    upper = mean + ci,
    .groups = "drop"
  )

# Friendly labels for facets
var_labels <- c(age = "Age", dem = "Pop. density")

# Plot
p5_ci <- ggplot(summary_df, aes(x = factor(cluster), y = mean, fill = factor(cluster))) +
  # error bars colored by cluster (no separate color legend)
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(cluster)),
                width = 0.4, position = pd, size = 0.9) +
  # points shape 21: fill mapped to cluster, black outline
  geom_point(shape = 21, size = 4, stroke = 1, position = pd, color = "black") +
  facet_wrap(~ variable, scales = "free_y", ncol = 2,
             labeller = as_labeller(var_labels)) +
  scale_fill_brewer(palette = "Dark2", name = "Cluster") +
  scale_color_brewer(palette = "Dark2", guide = "none") + # hide duplicate legend for color
  labs(
    title = "",
    x = "Cluster",
    y = "Mean Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print
p5_ci

## categorical variables
#library(colorspace)  # for lighten()
#library(binom)

plot_dodged_with_ci <- function(varname, title,
                                font_size = 14,
                                lighten_amount = 0.15) {
  # prepare data (keeps non-syntactic names like "wealth class")
  df <- malaysia_dpmm %>%
    dplyr::select(gender, stratum, wealth, occupation_grouped, cluster) %>%
    dplyr::select(
      gender = 1,
      `vegetation stratum` = 2,
      `wealth class` = 3,
      occupation = 4,
      cluster = 5
    ) %>%
    count(cluster, !!sym(varname)) %>%
    group_by(cluster) %>%
    mutate(total = sum(n),
           prop = n / total) %>%
    ungroup()
  
  # CIs (Wilson)
  cis <- binom::binom.confint(x = df$n, n = df$total, methods = "wilson")
  df <- df %>% mutate(lower = cis$lower * 100, upper = cis$upper * 100, prop100 = prop * 100)
  
  # Construct viridis palette and lighten the first color
  n_levels <- length(unique(df[[varname]]))
  pal <- viridis::viridis(n_levels, option = "D")    # you can change option
  pal[1] <- colorspace::lighten(pal[1], amount = lighten_amount)
  
  # Ensure mapping order matches the discrete values found in df
  # (this helps when varname contains spaces/non-syntactic names)
  values_order <- unique(df[[varname]])
  names(pal) <- values_order
  
  # Plot
  ggplot(df, aes(x = cluster, y = prop100, fill = !!sym(varname))) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.9),
                  width = 0.25, color = "black") +
    scale_fill_manual(values = pal) +
    labs(title = title,
         y = "Percentage", x = "Cluster", fill = varname) +
    theme_minimal(base_size = font_size) +
    theme(
      plot.title = element_text(size = font_size + 2, face = "bold"),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size - 2),
      legend.title = element_text(size = font_size),
      legend.text = element_text(size = font_size - 2)
    )
}

p1 <- plot_dodged_with_ci("gender", "a)", font_size = 15, lighten_amount = 0.2)
p2 <- plot_dodged_with_ci("wealth class", "b)", font_size = 15, lighten_amount = 0.2)
p3 <- plot_dodged_with_ci("vegetation stratum", "c)", font_size = 15, lighten_amount = 0.2)
p4 <- plot_dodged_with_ci("occupation", "d)", font_size = 15, lighten_amount = 0.2)

grid.arrange(p1, p2, p3, p4, nrow = 2)

###### statistical tests----
set.seed(42)

results_pvalues <- data.frame(
  Cluster = integer(), Variable = character(),
  P_Value = numeric(), Test_Used = character(),
  stringsAsFactors = FALSE
)

cat_vars <- c("gender", "veg.stratum", "wealth", "occupation")
cont_vars <- c("age", "pop.density")

for (k in 2:4) {
  temp_data <- malaysia_dpmm %>% mutate(is_target = ifelse(cluster == k, 1, 0))
  
  for (var in cat_vars) {
    tbl <- table(temp_data[[var]], temp_data$is_target)
    test_res <- fisher.test(tbl, simulate.p.value = TRUE, B = 10000)  # B aumentato per stabilità
    results_pvalues <- rbind(results_pvalues, data.frame(
      Cluster = k, Variable = var, P_Value = test_res$p.value,
      Test_Used = "Fisher Exact Test (simulated)"
    ))
  }
  
  for (var in cont_vars) {
    test_res <- wilcox.test(temp_data[[var]] ~ temp_data$is_target, exact = FALSE)
    results_pvalues <- rbind(results_pvalues, data.frame(
      Cluster = k, Variable = var, P_Value = test_res$p.value,
      Test_Used = "Wilcoxon Rank Sum Test"
    ))
  }
}

# Bonferroni per variabile, sui 3 confronti effettivamente fatti (cluster 2,3,4 - non 4)
results_pvalues <- results_pvalues %>%
  group_by(Variable) %>%
  mutate(P_Bonf = p.adjust(P_Value, method = "bonferroni")) %>%
  ungroup()

print(results_pvalues)

### maps----
# Read Malaysia shapefile
malaysia_map <- st_read("gadm41_MYS_shp/gadm41_MYS_2.shp")

# Filter to Sabah regions
sabah <- malaysia_map %>%
  filter(NAME_2 %in% c("Kudat", "Ranau", "Pitas", "Kota Marudu"))  %>%
  st_transform(crs = 32650)

# Convert your NTD dataframe into sf object (points in EPSG:32650)
NTD_sf <- st_as_sf(malaysia_dpmm, coords = c("x", "y"), crs = 32650)

##study area
sabah_base <- malaysia_map %>%
  filter(NAME_1 == "Sabah") %>%
  st_transform(crs = 32650)

ggplot() +
  geom_sf(data = sabah_base,
          fill  = "grey90",
          color = "black",
          size  = 0.2) +
  geom_sf(data = sabah,
          aes(fill = NAME_2),
          color = "black",
          size  = 0.2) +
  scale_fill_brewer(palette = "Set3", name = "District") +
  
  coord_sf(crs   = st_crs(32650),
           datum = NA,
           expand = FALSE) +
  
  # --- scale bar (metrics because CRS is in metres) ---
  annotation_scale(location = "br",    # "bl" = bottom-left; try "br","tl","tr"
                   width_hint = 0.25, # fraction of plot width the bar should take
                   bar_cols = c("grey60","white"),
                   text_cex = 0.8) +
  
  # --- north arrow ---
  annotation_north_arrow(location = "tr",    # "tr" = top-right
                         which_north = "true", # "true" or "grid"
                         pad_x = unit(0.5, "cm"),
                         pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering()) +
  
  theme_minimal() +
  labs(
    title = "",
    x     = "",
    y     = ""
  )

# Filter NTD_sf to only include points from cluster 2
NTD_cluster2 <- NTD_sf %>%
  filter(cluster == 2)

# Plot only cluster 2 points over the Sabah map
ggplot() +
  geom_sf(data = sabah, fill = "white", color = "black") + # Plot background map
  geom_sf(data = NTD_cluster2, color = "#D95F02", size = 2) + # Points for cluster 3
  theme_minimal() +
  labs(title = "Cluster 2 Locations over Sabah Regions", x = "Longitude", y = "Latitude")


#### Spatial distribution ----

###### df con variabili dummy  
dummy_matrix <- model.matrix(~ cluster-1, data = malaysia_dpmm)
malaysiadpmm_dummies <- cbind(malaysia_dpmm, dummy_matrix)

###### load a mesh 
mesh_sabah<-readRDS("mesh_sabah.rds")

###### define spde
spde <- inla.spde2.pcmatern(
  mesh = mesh_sabah,
  prior.range = c(20, 0.01),   # P(range < 20) = 0.01 
  prior.sigma = c(1, 0.01)    # P(sigma > 1) = 0.01
)

###### create a prediction grid 
# Create a grid
# Define the resolution of your grid (in meters, since CRS is EPSG:32650)
grid_spacing <- 1000  # e.g., 1 km between points

sabah_union <- st_union(sabah)

# Get the bounding box of Sabah
sabah_bbox <- st_bbox(sabah_union)

# Create a regular grid of points over the bounding box
x_coords <- seq(sabah_bbox$xmin, sabah_bbox$xmax, by = grid_spacing)
y_coords <- seq(sabah_bbox$ymin, sabah_bbox$ymax, by = grid_spacing)
grid <- expand.grid(x = x_coords, y = y_coords)

# Convert the grid to an sf object
grid_points <- st_as_sf(grid, coords = c("x", "y"), crs = 32650)

# Keep only points within the Sabah region
grid_within_sabah <- grid_points[st_within(grid_points, sabah_union, sparse = FALSE), ]

# Plot to check
ggplot() +
  geom_sf(data = sabah, fill = "white", color = "black") +
  geom_sf(data = grid_within_sabah, color = "blue", size = 0.5, alpha = 0.6) +
  labs(title = "1 km × 1 km prediction grid within Sabah")+
  theme_minimal()
  

# Convert sf grid to SpatialPointsDataFrame
grid_sp <- as(grid_within_sabah, "Spatial")
# Create matrix of coordinates
grid_sp@data$coordinates <- coordinates(grid_sp)

#### plot the predicted probabilities of belonging to a cluster

# Lists to store results
models_list <- vector("list", 4)
predictions_list <- vector("list", 4)

# Loop over 1 to 4
models_and_preds <- map(1:4, function(i) {
  # Dynamic variable names
  cluster_col <- glue("cluster{i}")
  cluster_var <- malaysiadpmm_dummies[, c("x", "y", cluster_col)]
  
  # Convert to sf and reproject
  cluster_sf <- st_as_sf(cluster_var, coords = c("x", "y"), crs = 32650)
  cluster_sp <- as(cluster_sf, "Spatial")
  
  # Define formula dynamically
  formula <- as.formula(paste0(cluster_col, " ~ Intercept(1) + spatial_field(coordinates, model = spde)"))
  
  # Fit the model
  model_fit <- bru(
    components = formula,
    family = "binomial",
    data = cluster_sp,
    options = list(control.compute = list(dic = TRUE, waic = TRUE))
  )
  
  # Predict
  prediction <- predict(
    model_fit,
    newdata = grid_sp,
    formula = ~ 1 / (1 + exp(-(Intercept + spatial_field)))
  )
  
  # Return model and prediction
  list(model = model_fit, prediction = prediction)
})

# Separate into two lists
models_list <- map(models_and_preds, "model")
save(models_list, file = "m0_allcluster_dpmm1.RData")

predictions_list <- map(models_and_preds, "prediction")
save(predictions_list, file = "m0pred_allcluster_dpmm.RData")

load("m0pred_allcluster_dpmm.RData")

#plot
grid_coords <- coordinates(grid_sp) %>%
  as.data.frame()
colnames(grid_coords) <- c("X", "Y")

# Combine predictions into one dataframe
all_preds_df <- map2_dfr(
  .x = predictions_list,
  .y = 1:4,
  .f = function(pred, i) {
    tibble(
      X = grid_coords$X,
      Y = grid_coords$Y,
      pred = pred$mean,
      cluster = glue("Cluster {i}")
    )
  }
)

all_preds_df$cluster <- factor(all_preds_df$cluster, levels = paste0("Cluster ", 1:4))

# Plot 
ggplot(all_preds_df, aes(X, Y, fill = pred)) +
  geom_tile() +
  geom_sf(data = sabah, inherit.aes = FALSE, fill = NA, color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", name = "Predicted\nProbability") +
  coord_sf() +
  theme_minimal(base_size = 10) +
  labs(title = "Estimated Probability of Belonging to Clusters") +
  facet_wrap(~ cluster, ncol = 3, nrow = 2)

# Generate individual plots per cluster
plots <- map(1:4, function(i) {
  cluster_df <- all_preds_df %>% filter(cluster == paste0("Cluster ", i))
  
  ggplot(cluster_df, aes(X, Y, fill = pred)) +
    geom_tile() +
    geom_sf(data = sabah, inherit.aes = FALSE, fill = NA, color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "magma", name = "Predicted\nProbability") +
    coord_sf() +
    theme_minimal(base_size = 10) +
    labs(title = paste("Estimated Probability of Belonging to Clusters", i))
})

p1<- plots[[1]]
p2<- plots[[2]]
p3<- plots[[3]]
p4<- plots[[4]]
gridExtra::grid.arrange(p1, p3, p4, nrow=1)
p2

