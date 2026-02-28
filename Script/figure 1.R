#####################################
## @Description: Generate Figure 1 (NiV outbreaks, reservoir distribution,
##               phylogeny) and derive reservoir risk indices and
##               seasonal susceptibility model.
## @version: 1.0
## @Author: Li Kangguo
## @Date: 2026-02-05 22:24:56
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-02-28 11:03:43
#####################################

library(sf)
library(patchwork)
library(ggtree)
library(treeio)
library(tidyverse)
library(ggrepel)
library(ggtreeExtra) 
library(phytools)
library(ggthemes)
library(ggspatial)
library(maptiles)
p_deaths <- ggplot(data_panel_A) +
     geom_col(aes(x = Year, y = deaths, fill = country),
              show.legend = T) +
     scale_x_continuous(limits = c(1997, 2027),
                        expand = expansion(add = c(0.5, 0.5)),
                        breaks = seq(1998, 2026, by = 5)) +
     scale_y_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 150),
                        breaks = scales::pretty_breaks(n = 5)) +
     scale_fill_manual(values = fill_color,
                       breaks = names(fill_color),
                       drop = FALSE,
                       name = 'Country') +
     labs(title = NULL,
          x = 'Year',
          fill = 'Country',
          y = 'Number of deaths') +
     theme_bw() +
     theme(legend.position = 'bottom',
           panel.grid = element_blank())+
     guides(fill = guide_legend(nrow = 1))

# panel B -----------------------------------------------------------------

basemap <- st_read(globalmap_shp, quiet = TRUE)
china_border <- st_read(chinamap_shp, quiet = TRUE)

data_panel_B <- read.delim(gbif_dir, sep = "\t", header = TRUE, quote = "")

species_mapping <- data_panel_B |> 
     group_by(species) |>
     summarise(count = n(),
               .groups = 'drop') |>
     arrange(desc(count)) |> 
     mutate(species_clean = if_else(species == '', 'Unknown', species),
            species_clean = if_else(count > 1000, species_clean, 'Other species'))

species_levels <- species_mapping |>
     group_by(species_clean) |>
     summarise(total = sum(count)) |>
     arrange(desc(total)) |>
     pull(species_clean)

species_levels <- c(species_levels[!species_levels %in% c('Other species', 'Unknown')], 'Other species', 'Unknown')

data_panel_B <- data_panel_B |>
     left_join(species_mapping, by = 'species') |> 
     mutate(species_clean = factor(species_clean, levels = species_levels))

# manual fallbacks for countries not found in basemap
fallbacks <- tibble(
     country = c('Bangladesh','India','Malaysia','Singapore','Philippines','Thailand','Madagascar','China','Indonesia','Australia'),
     lon = c(90.4125,77.2090,101.6869,103.8198,120.9842,100.5018,47.5079,104.1954,106.8456,130),
     lat = c(23.8103,28.6139,5,1,14.5995,13.7563,-18.8792,35.8617,-6.2088,-35.2809)
)

# fetch Satellite tiles for background (Esri.WorldImagery)
osm_bbox <- st_bbox(c(xmin = 30, xmax = 160, ymin = -50, ymax = 50), crs = st_crs(4326))
osm_area <- st_as_sfc(osm_bbox)

osm_tiles <- get_tiles(
     x = osm_area,
     provider = "Esri.WorldImagery",
     zoom = 4,
     crop = TRUE,
     project = FALSE
)

pB <- ggplot() +
     layer_spatial(osm_tiles) +
     # geom_sf(data = basemap, fill = 'grey95', color = 'grey60') +
     # geom_sf(data = china_border, fill = 'white', color = 'black', size = 0.1) +
     # add country labels (no basemap shown)
     # geom_point(data = fallbacks, aes(x = lon, y = lat), inherit.aes = FALSE, colour = 'black', size = 1) +
     geom_text(data = fallbacks, aes(x = lon, y = lat, label = country), inherit.aes = FALSE,
               size = 4, color = "white", fontface = "bold") +
     geom_point(data = data_panel_B,
               aes(x = decimalLongitude, y = decimalLatitude, color = species_clean),
               size = 0.5, alpha = 0.6) +
     coord_sf(xlim = c(30, 160), ylim = c(-50, 50), crs = 4326, expand = FALSE) +
     scale_x_continuous(breaks = seq(30, 160, by = 20),
                        expand = c(0, 0),
                        labels = function(x) paste0(x, "°E")) +
     scale_y_continuous(breaks = seq(-50, 50, by = 20),
                        expand = c(0, 0),
                        labels = function(y) ifelse(y > 0, paste0(y, "°N"), ifelse(y < 0, paste0(abs(y), "°S"), "0°"))) +
     scale_color_viridis_d()+
     labs(title = 'B', x = "Longitude", y = "Latitude") +
     theme_bw()+
     theme(legend.text = element_text(face = 'italic'),
           legend.position = 'right',
           plot.title.position = 'plot',
           panel.grid.major = element_line(color = "white", linewidth = 0.1),
           panel.grid.minor = element_blank())+
     guides(color = guide_legend(title = 'Species', override.aes = list(size = 4)))

# panel C -----------------------------------------------------------------

tr <- read.tree(ncbi_tree)
meta <- read_tsv(ncbi_meta, show_col_types = FALSE)

meta2 <- meta |>
     transmute(label = as.character(Accession),
               Country = as.character(Country)) |>
     distinct(label, .keep_all = TRUE) |>
     mutate(Country = if_else(is.na(Country) | Country == "", "Unknown", Country))

tip_df <- tibble(label = tr$tip.label,
                 acc = str_replace(label, "\\.[0-9]+$", "")) |> 
     left_join(meta2, by = c("acc" = "label")) |>
     mutate(Country = if_else(is.na(Country) | Country == "", "Unknown", Country))

tip_df2 <- tip_df |>
     left_join(meta |> transmute(acc = as.character(Accession),
                                  Collection_Date = as.character(Collection_Date)),
               by = "acc") |>
     mutate(
          year = str_extract(Collection_Date, "(19|20)[0-9]{2}"),
          year = as.integer(year),
          acc = if_else(!Country %in% c("Bangladesh", "India", "Malaysia"), 
                        paste0(acc, " (", Country, ")"),
                        acc)
     )

tr2 <- tr
tr2 <- unroot(tr2)
tr2 <- midpoint.root(tr2) 

p <- ggtree(tr2) %<+% tip_df2 +
     geom_tiplab(aes(color = Country, label = acc), size = 2, show.legend = F) +
     scale_color_manual(values = fill_color, na.value = "grey60") +
     theme_tree2() +
     labs(color = "Country") +
     theme(legend.position = "none") +
     scale_x_continuous(trans = "sqrt", expand = expansion(mult = c(0.01, 0.02)))

heat_df <- tip_df2 |>
     select(label, year) |>
     distinct(label, .keep_all = TRUE) |>
     tibble::column_to_rownames("label")

pC <- gheatmap(p,
               heat_df,
               offset = 0.03,
               width  = 0.1,
               colnames = F,
               colnames_angle = 90,
               font.size = 3) +
     scale_fill_viridis_c(na.value = "grey85",
                          breaks = seq(1998, 2026, by = 3),
                          guide = guide_colorbar(title = "Collection Year", barwidth = 1, barheight = 12),
                          name = "Year") +
     theme(legend.position = "inside",
           legend.position.inside = c(0.3, 1),
           legend.justification = c(0, 1),
           plot.margin = margin(5.5, 30, 5.5, 5.5))+
     labs(title = 'C')

# Combine panels ----------------------------------------------------------

final_plot <- cowplot::plot_grid(
     free(p_cases / p_deaths) + pB + plot_layout(ncol = 1, heights = c(1, 1)),
     pC,
     ncol = 2,
     rel_widths = c(1.5, 1)
)

ggsave(filename = './Outcome/fig1.png',
       plot = final_plot,
       width = 12,
       height = 8,
       dpi = 300)

# Risk Calculation (Eco-Viral Weighted Index) -----------------------------

# Logic: Combine "Potential Host Presence" (GBIF) vs "Confirmed Virus Activity" (NCBI)
# - GBIF Data: Provides ecological potential (Bat abundance/diversity).
# - NCBI Data: Provides confirmed viral spillover evidence (Sequences from outbreaks).

# Define Country Map
country_iso_map <- c(
     "BD" = "Bangladesh", "IN" = "India", "MY" = "Malaysia", 
     "SG" = "Singapore", "PH" = "Philippines", "TH" = "Thailand"
)

# 1. Host Score (Ecological Potential from GBIF)
gbif_risk <- data_panel_B |>
     filter(countryCode %in% names(country_iso_map)) |>
     mutate(Country = country_iso_map[as.character(countryCode)]) |>
     group_by(Country) |>
     summarise(
          GBIF_Recs = n(),
          GBIF_Spp = n_distinct(species),
          .groups = "drop"
     )

# 2. Virus Score (Confirmed Spillover Evidence from NCBI)
# Using 'tip_df2' generated for Panel C
ncbi_risk <- tip_df2 |>
     filter(Country %in% unname(country_iso_map)) |> # Ensure we match project countries
     group_by(Country) |>
     summarise(
          NCBI_Seqs = n(),
          Last_Seq_Year = max(year, na.rm = TRUE),
          .groups = "drop"
     )

# 3. Composite Weighted Risk Calculation (Robust Z-Score Method)
combined_risk <- tibble(Country = unname(country_iso_map)) |>
     left_join(gbif_risk, by = "Country") |>
     left_join(ncbi_risk, by = "Country") |>
     replace_na(list(
          GBIF_Recs = 0, 
          GBIF_Spp = 0, 
          NCBI_Seqs = 0, 
          Last_Seq_Year = 2000
    )) |>
     mutate(
          # [Cor.1] Log-transform using consistent base (natural log)
          log_Host  = log1p(GBIF_Recs),
          log_Virus = log1p(NCBI_Seqs),
          
          # [Cor.2] Robust Z-Score Calculation using scale()
          z_Host_val  = as.numeric(scale(log_Host)),
          z_Virus_val = as.numeric(scale(log_Virus)),
          
          # If all values are identical (SD=0), scale() produces NaN. Replace with 0.
          z_Host = if_else(is.nan(z_Host_val), 0, z_Host_val),
          z_Virus = if_else(is.nan(z_Virus_val), 0, z_Virus_val),

          # [Cor.3] Weighted Combination (NCBI=2x confirmed evidence)
          Raw_Z_Score = (1 * z_Host + 2 * z_Virus) / 3,

          # [Cor.4] Normalize to [0.01, 1]
          Reservoir_Risk_Calculated = scales::rescale(Raw_Z_Score, to = c(0.01, 1))
     ) |>
     arrange(desc(Reservoir_Risk_Calculated))

write.csv(combined_risk, file = "./Outcome/risk_assessment.csv", row.names = FALSE)

## Fit Biological Seasonality Function: Bootstrap resampling + GAM with cyclic splines
library(mgcv)
library(purrr)

# Prepare monthly-aggregated cases
cases_monthly <- data_panel_A |>
     filter(!is.na(Date)) |>
     mutate(Year = year(Date), Month = month(Date)) |>
     group_by(Year, Month) |>
     summarise(Cases = sum(cases, na.rm = TRUE), .groups = 'drop')

# Ensure months 1..12 appear at least in summary table
monthly_avg <- cases_monthly |>
     group_by(Month) |>
     summarise(AvgCases = mean(Cases, na.rm = TRUE),
               SDcases = sd(Cases, na.rm = TRUE), .groups = 'drop') |>
     complete(Month = 1:12, fill = list(AvgCases = 0, SDcases = 0)) |>
     arrange(Month)

# Bootstrap parameters
B <- 1000
set.seed(20260228)
months <- 1:12

# Container for bootstrap predictions
boot_preds <- matrix(NA, nrow = B, ncol = length(months))

for (b in seq_len(B)) {
     # resample months (rows) with replacement
     idx <- sample(nrow(cases_monthly), replace = TRUE)
     boot_data <- cases_monthly[idx, , drop = FALSE]

     # ensure there are enough unique month values to fit a cyclic spline
     unique_months <- length(unique(boot_data$Month))
     if (unique_months < 3) {
          # not enough unique months to fit a cyclic spline; record NA and continue
          boot_preds[b, ] <- NA_real_
          next
     }

     # choose k adaptively: cannot exceed number of unique months
     k_b <- min(6, unique_months)

     # Fit GAM with cyclic spline; use Negative Binomial, fallback to Poisson
     fit_b <- tryCatch({
          gam(Cases ~ s(Month, bs = "cc", k = k_b), data = boot_data, family = nb(), method = "REML")
     }, error = function(e) {
          tryCatch({
               gam(Cases ~ s(Month, bs = "cc", k = k_b), data = boot_data, family = poisson(), method = "REML")
          }, error = function(e2) {
               return(NULL)
          })
     })

     if (is.null(fit_b)) {
          boot_preds[b, ] <- NA_real_
          next
     }

     # Predict expected counts for months 1..12 (if model fit succeeded)
     pred_b <- predict(fit_b, newdata = data.frame(Month = months), type = "response")
     boot_preds[b, ] <- as.numeric(pred_b)
}

# Summarize bootstrap predictions: mean and 95% CI
pred_mean  <- apply(boot_preds, 2, mean, na.rm = TRUE)
pred_lower <- apply(boot_preds, 2, quantile, probs = 0.025, na.rm = TRUE)
pred_upper <- apply(boot_preds, 2, quantile, probs = 0.975, na.rm = TRUE)

# Prediction wrapper returns bootstrap mean by month
seasonal_predict <- function(month_vec) {
     pred_mean[((month_vec - 1) %% 12) + 1]
}

seasonal_model <- list(
     boot_preds = boot_preds,
     pred_mean = pred_mean,
     pred_lower = pred_lower,
     pred_upper = pred_upper,
     monthly_avg = monthly_avg,
     predict = seasonal_predict
)

save(seasonal_model, seasonal_predict, file = "./Outcome/seasonal_model.RData")
