#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-02-05 22:24:56
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-02-05 22:25:01
#####################################

library(sf)
library(patchwork)
library(ggtree)
library(treeio)
library(tidyverse)
library(ggtreeExtra) 

# Output directory for figures
out_dir <- "./Outcome/"
outbreaks_path <- "./Data/Outbreak/outbreaks_nipah.csv"
gbif_dir <- "./Data/GBIF/0013514-260129131611470.csv"
globalmap_shp <- "./Data/globalmap/globalmap.shp"
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

pB <- ggplot(data = data_panel_B) +
     geom_sf(data = basemap, fill = 'grey95', color = 'grey60') +
     geom_point(aes(x = decimalLongitude, y = decimalLatitude, color = species_clean), size = 0.05, alpha = 0.3)+
     coord_sf(xlim = c(35, 150), ylim = c(-50, 50)) +
     scale_fill_viridis_d()+
     labs(title = 'B') +
     theme_void()+
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
          year = as.integer(year)
     )

p <- ggtree(tr) %<+% tip_df2 +
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

pC <- gheatmap(p, heat_df,
               offset = 0.02,
               width  = 0.1,
               colnames = F,
               colnames_angle = 90,
               font.size = 3) +
     scale_fill_viridis_c(na.value = "grey85",
                          name = "Year") +
     theme(legend.position = "inside",
           legend.position.inside = c(0.9, 1),
           legend.justification = c(1, 1),
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



