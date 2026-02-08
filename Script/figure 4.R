#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-02-04 15:04:57
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-02-04 18:01:07
#####################################

library(openxlsx)
library(tidyverse)
library(forecast)
library(patchwork)
library(scales)

# Merge Global context into Countries
load("./Outcome/figure2_workspace.RData")

Sys.setlocale("LC_TIME", "English")

scientific_10 <- function(x) {
     ifelse(x == 0, 0, parse(text = gsub("[+]", "", gsub("e", "%*%10^", scales::scientific_format()(x)))))
}

# 2. Risk Analysis (Probabilistic Exposure Model)
# -----------------------------------------------------------------------------
# Combine three components to compute an exposure score for each origin-month:
#  - Inbound_Volume: passenger flows from origin.
#  - Epidemiological context: per-country `Reservoir_Risk`, `Spillover_Activity`,
#    and `Detection_Capacity` (from ./Outcome/fig3_risk_metadata.csv).
#  - Seasonal_Index: estimated monthly susceptibility derived from the
#    previously fitted seasonal model (./Outcome/seasonal_model.RData).
#
# Exposure_Score = Inbound_Volume * Reservoir_Risk * Spillover_Activity * Detection_Capacity * Seasonal_Index

# Load previously computed epidemiological parameters
# We use the per-country values calculated earlier and exported to
# `./Outcome/fig3_risk_metadata.csv`. That file contains:
#   Country, Reservoir_Risk, Spillover_Activity, Detection_Capacity
# We map `Country` -> `Origin_Country` for joining with `df_final`.
risk_metadata <- read.csv("./Outcome/fig3_risk_metadata.csv", stringsAsFactors = FALSE) |>
     rename(Origin_Country = Country)

# Load estimated seasonal model and create a 1..12 seasonal index.
# `seasonal_model$predict(1:12)` returns expected monthly averages;
# we rescale these values to the index range used in the model (1..3).
if (file.exists("./Outcome/seasonal_model.RData")) {
     load("./Outcome/seasonal_model.RData")
     if (exists("seasonal_model") && is.list(seasonal_model) && is.function(seasonal_model$predict)) {
          seasonal_pred_vals <- seasonal_model$predict(1:12)
          seasonal_index_vec <- scales::rescale(seasonal_pred_vals, to = c(1, 3))
     } else {
          seasonal_index_vec <- rep(1, 12)
     }
} else {
     seasonal_index_vec <- rep(1, 12)
}

df_analyzed <- df_final |>
     left_join(risk_metadata, by = "Origin_Country") |>
     mutate(
          # Fall back to conservative defaults when metadata is missing
          Reservoir_Risk = replace_na(Reservoir_Risk, 1.0),
          Spillover_Activity = replace_na(Spillover_Activity, 0.01),
          Detection_Capacity = replace_na(Detection_Capacity, 0.5),

          # Seasonal factor: use the estimated monthly index
          Seasonal_Index = seasonal_index_vec[((Month - 1) %% 12) + 1],

          # Composite exposure score using Detection_Capacity (replaces prior Endemic_Connectivity)
          Exposure_Score = Inbound_Volume * Reservoir_Risk * Spillover_Activity * Detection_Capacity * Seasonal_Index
     )

# 3. Risk Grading System (Statistical Anomaly Standard)
# ==============================================================================
# Objective Grading Methodology:
# Instead of forcing a "bell curve" where top countries are always "High/Critical"
# (Relative Grading), we use Statistical Process Control (SPC) principles.
#
# We establish a "Baseline" using the distribution of scores.
# - Levels are defined by Standard Deviations (SD) from the Log-Mean.
# - This ensures "High" & "Critical" labels exist theoretically, but are only
#   triggered by significant deviations (Outliers), satisfying the observation
#   that routine inputs should be "Low/Moderate".
# ==============================================================================

# Calculate distribution parameters (Using all history + forecast context)
# We use Log-transformation to normalize the skewed volume data
risk_stats <- df_analyzed |>
     filter(Exposure_Score > 0) |>
     summarise(
          Log_Mean = mean(log10(Exposure_Score + 1)),
          Log_SD = sd(log10(Exposure_Score + 1))
     )

thresh_base <- risk_stats$Log_Mean
thresh_step <- risk_stats$Log_SD

df_analyzed <- df_analyzed |>
     mutate(
          # Z-Score helps quantify "How extreme is this month?"
          Log_Score = log10(Exposure_Score + 1),
          Risk_Level = case_when(
               # > 3 SD from mean (Extreme Anomaly) -> Critical
               Log_Score > (thresh_base + 3 * thresh_step) ~ "Level 5: Critical",
               
               # > 2 SD from mean (Rare High) -> High
               Log_Score > (thresh_base + 2 * thresh_step) ~ "Level 4: High",
               
               # > 1 SD from mean (Seasonal Peak/High Vol) -> Moderate
               Log_Score > (thresh_base + 1 * thresh_step) ~ "Level 3: Moderate",
               
               # > Mean -> Low
               Log_Score > thresh_base ~ "Level 2: Low",
               
               # Following Mean -> Very Low
               TRUE ~ "Level 1: Very Low"
          )
     )

# Order factors explicitly (Restoring full risk spectrum)
df_analyzed$Risk_Level <- factor(df_analyzed$Risk_Level, levels = c(
     "Level 1: Very Low", "Level 2: Low", "Level 3: Moderate", "Level 4: High", "Level 5: Critical"
))

# Calculate Share of Global
df_analyzed <- df_analyzed |>
     mutate(Share_of_Global = Inbound_Volume / Global_Inbound_Volume)

# 5. Visualization 2: Risk Assessment
#####################################

plot_start_year <- 2024
df_plot_risk <- df_analyzed |> filter(Year >= 2022)

fill_colors <- RColorBrewer::brewer.pal(n = 5, name = "RdYlGn")[5:1] # Reverse for proper mapping
names(fill_colors) <- levels(df_plot_risk$Risk_Level)

# Plot Risk Trend (Line + Points determined by Color)
# Y-axis = Inbound_Volume (Flow), Color = Risk_Level
p_risk_trend <- ggplot(df_plot_risk, aes(x = Date_Obj, y = Inbound_Volume, group = Origin_Country)) +
     annotate("rect", xmin = as.Date(last_date), xmax = max(df_plot_A$Date_Obj),
              ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
     geom_line(color = "grey70", linewidth = 0.5) +
     geom_point(aes(color = Risk_Level), size = 2.5, show.legend = T) +
     facet_wrap(~Origin_Country, scales = "free_y", ncol = 1) +
     # Ensure all risk levels appear in the legend even if unused
     scale_color_manual(values = fill_colors, drop = FALSE) +
     scale_y_continuous(labels = scientific_10,
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.15))) +
     scale_x_date(date_labels = "%b %Y",
                  expand = c(0, 0),
                  date_breaks = "6 months") +
     geom_vline(xintercept = as.numeric(last_date), linetype = "dashed", color = "red") +
     labs(
          title = 'A',
          x = NULL, 
          y = "Monthly volume",
          color = "Risk Level"
     ) +
     theme_bw() +
     theme(legend.position = "bottom",
           plot.title.position = "plot")

# Plot Risk Forecast Bar (Details)
df_forecast_only <- df_plot_risk |> filter(Type == "Forecast")

p_risk_bar <- ggplot(df_forecast_only, aes(x = as.factor(format(Date_Obj, "%b %Y")), y = Exposure_Score, fill = Risk_Level)) +
     geom_col(show.legend = F) +
     facet_wrap(~Origin_Country, scales = "free_y", ncol = 1) +
     # Keep full color scale available for all risk levels
     scale_fill_brewer(palette = "RdYlGn", direction = -1, drop = FALSE) +
     scale_y_continuous(labels = scientific_10,
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.15))) +
     labs(
          title = "B",
          x = "Date", y = "Composite risk score"
     ) +
     theme_bw() +
     theme(
          legend.position = "none",
          plot.title.position = "plot"
     )

# Compose Figure 2
fig2 <- p_risk_trend + p_risk_bar + plot_layout(widths = c(2, 1), guides = 'collect') &
     theme(legend.position = 'bottom')

ggsave("./Outcome/fig4.png",
       fig2,
       width = 16, 
       height = 10, 
       bg="white")

