#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-02-05 22:24:56
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-02-06 11:34:25
#####################################

library(sf)
library(patchwork)
library(ggtree)
library(treeio)
library(tidyverse)
library(ggrepel)
library(ggtreeExtra) 
library(phytools)

# Output directory for figures
out_dir <- "./Outcome/"
outbreaks_path <- "./Data/Outbreak/outbreaks_nipah.csv"
gbif_dir <- "./Data/GBIF/0013514-260129131611470.csv"
globalmap_shp <- "./Data/globalmap/globalmap.shp"
chinamap_shp <- "./Data/globalmap/china_border.shp"
ncbi_tree <- "./Data/NCBI/nipah_80.aln.trim.fasta.treefile"
ncbi_meta <- "./Data/NCBI/sequences.tsv"

fill_color <- c("#462255FF", "#FF8811FF", "#9DD9D2FF", "#046E8FFF", "#D44D5CFF", "grey50")
names(fill_color) <- c("Bangladesh", "India", "Malaysia", "Singapore", "Philippines", "Thailand")

# panel A -----------------------------------------------------------------

data_panel_A <- read.csv(outbreaks_path) |> 
     mutate(cases = as.numeric(gsub('[^0-9]', '', cases)),
            deaths = as.numeric(gsub('[^0-9]', '', deaths)),
            # extract content before -
            date = gsub('-.*$', '', started),
            Date = mdy(date),
            Year = year(Date),
            id = ifelse(id == '', NA_character_, id),
            country = factor(country, levels = names(fill_color)))

p_cases <- ggplot(data_panel_A) +
     geom_col(aes(x = Year, y = cases, fill = country)) +
     scale_x_continuous(limits = c(1997, 2027),
                        expand = expansion(add = c(0.5, 0.5)),
                        breaks = seq(1998, 2026, by = 5)) +
     scale_y_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 300),
                        breaks = scales::pretty_breaks(n = 5)) +
     scale_fill_manual(values = fill_color,
                       name = 'Country') +
     labs(title = 'A',
          x = NULL,
          y = 'Number of Cases') +
     theme_bw() +
     theme(legend.position = 'none',
           plot.title.position = 'plot',
           panel.grid = element_blank(),
           axis.text.x = element_blank())

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
          y = 'Number of Deaths') +
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
     lon = c(90.4125,77.2090,101.6869,103.8198,120.9842,100.5018,47.5079,104.1954,106.8456,149.1300),
     lat = c(23.8103,28.6139,5,1,14.5995,13.7563,-18.8792,35.8617,-6.2088,-35.2809)
)

pB <- ggplot(data = data_panel_B) +
     # geom_sf(data = basemap, fill = 'grey95', color = 'grey60') +
     # geom_sf(data = china_border, fill = 'white', color = 'black', size = 0.1) +
     # add country labels (no basemap shown)
     # geom_point(data = fallbacks, aes(x = lon, y = lat), inherit.aes = FALSE, colour = 'black', size = 1) +
     geom_text(data = fallbacks, aes(x = lon, y = lat, label = country), inherit.aes = FALSE,
               size = 4) +
     geom_point(aes(x = decimalLongitude, y = decimalLatitude, color = species_clean), size = 0.05, alpha = 0.3)+
     coord_sf(xlim = c(35, 150), ylim = c(-50, 50)) +
     scale_fill_viridis_d()+
     labs(title = 'B') +
     theme_void()+
     theme(legend.text = element_text(face = 'italic'))+
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
                          name = "Year") +
     theme(legend.position = "inside",
           legend.position.inside = c(0, 1),
           legend.justification = c(0, 1),
           plot.margin = margin(5.5, 30, 5.5, 5.5))+
     labs(title = 'C')

# Combine panels ----------------------------------------------------------

final_plot <- cowplot::plot_grid(
     free(p_cases / p_deaths) + pB + plot_layout(ncol = 1, heights = c(1, 1)),
     pC,
     ncol = 2,
     rel_widths = c(1.2, 1)
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
     mutate(Country = country_iso_map[countryCode]) |>
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

# 3. Composite Weighted Risk Calculation
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
          # A. Base Score: Ecological Presence (Log-dampened to avoid sampling bias)
          # Weight: 1.0
          Score_Host = log10(GBIF_Recs + 1),
          
          # B. Active Score: Viral Confirmation (Root-dampened but higher impact)
          # Weight: 2.0 (Viral evidence > Host evidence)
          # Explanation: 1 sequence proves presence more definitively than 100 bat sightings.
          Score_Virus = sqrt(NCBI_Seqs) * 2, 
          
          # C. Final Synthesis
          Raw_Score = Score_Host + Score_Virus,
          
          # D. Scaling to 1.0 - 5.0 (Model Input Range)
          Reservoir_Risk_Calculated = scales::rescale(Raw_Score, to = c(1, 5))
     ) |>
     arrange(desc(Reservoir_Risk_Calculated))

write.csv(combined_risk, file = "./Outcome/risk_assessment.csv", row.names = FALSE)

# Fit Biological Seasonality Function -------------------------------------

# We fit a simple harmonic regression on monthly-aggregated outbreak cases
# using the historical outbreak table `data_panel_A`. The model is fit on
# the log(1 + mean monthly cases) to stabilize variance. We save the fitted
# model (and a small wrapper predict function) to `./Outcome/seasonal_model.rds`.

# Prepare monthly-aggregated cases (Year-Month grid)
cases_monthly <- data_panel_A |>
     filter(!is.na(Date)) |>
     mutate(Year = year(Date), Month = month(Date)) |>
     group_by(Year, Month) |>
     summarise(Cases = sum(cases, na.rm = TRUE), .groups = 'drop')

# Compute average cases per calendar month across years
monthly_avg <- cases_monthly |>
     group_by(Month) |>
     summarise(AvgCases = mean(Cases, na.rm = TRUE),
               SDcases = sd(Cases, na.rm = TRUE), .groups = 'drop') |>
     arrange(Month)

# Fit harmonic regression on log1p(AvgCases)
# Model: log1p(AvgCases) ~ sin(2*pi*Month/12) + cos(2*pi*Month/12)
harmonic_fit <- tryCatch({
     lm(log1p(AvgCases) ~ sin(2 * pi * Month / 12) + cos(2 * pi * Month / 12), data = monthly_avg)
}, error = function(e) {
     # Fallback to simple mean model if fit fails
     lm(log1p(AvgCases) ~ 1, data = monthly_avg)
})

# Create a small prediction wrapper that returns expected monthly multiplicative factor
seasonal_predict <- function(month_vec) {
     # month_vec can be numeric vector (1..12)
     newdf <- data.frame(Month = as.integer((month_vec - 1) %% 12 + 1))
     pred_log <- predict(harmonic_fit, newdata = newdf)
     # convert back from log1p
     pmax(0, exp(pred_log) - 1)
}

# Save model and supporting data
seasonal_model <- list(
     fit = harmonic_fit,
     monthly_avg = monthly_avg,
     predict = seasonal_predict
)

save(seasonal_model, harmonic_fit, seasonal_predict, monthly_avg, file = "./Outcome/seasonal_model.RData")
