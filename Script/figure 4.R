#####################################
## @Description: Generate Figure 4 (probabilistic importation risk
##               analysis and SPC-based grading) using outputs from
##               Figures 2 and 3.
## @version: 1.0
## @Author: Li Kangguo
## @Date: 2026-02-04 15:04:57
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-02-28 11:45:12
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
# Exposure_Score = Volume_Score * Reservoir_Risk * Spillover_Activity * Detection_Capacity * Seasonal_Index
#   (all five factors are normalized to [0.01, 1]; Volume_Score = rescale(log1p(Inbound_Volume)))

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
load("./Outcome/seasonal_model.RData")
# Use bootstrap mean predictions if present
seasonal_pred_vals <- if (!is.null(seasonal_model$pred_mean)) seasonal_model$pred_mean else seasonal_model$predict(1:12)
# With bootstrap GAM output (raw expected counts), normalize to [0.1, 1] relative amplitude
seasonal_index_vec <- scales::rescale(seasonal_pred_vals, to = c(0.1, 1))

df_analyzed <- df_final |>
     left_join(risk_metadata, by = "Origin_Country") |>
     mutate(
          # Fall back to neutral mid-range defaults when risk metadata is missing.
          # 0.5 represents unknown-but-non-negligible risk (not zero, not maximum).
          Reservoir_Risk      = replace_na(Reservoir_Risk, 0.5),
          Spillover_Activity  = replace_na(Spillover_Activity, 0.1),
          Detection_Capacity  = replace_na(Detection_Capacity, 0.5),

          # Seasonal factor: monthly index normalized to [0.1, 1]
          Seasonal_Index = seasonal_index_vec[((Month - 1) %% 12) + 1],

          # Normalize inbound volume to [0.01, 1]:
          #   Step 1 — Ensure Inbound_Volume is non-negative (Holt-Winters artifact safety)
          Inbound_Volume = pmax(Inbound_Volume, 0),
          #   Step 2 — log1p transformation compresses the large variance across countries.
          Log_Volume   = log1p(Inbound_Volume),
          #   Step 3 — min-max rescaling to [0.01, 1].
          #   Note: We ungroup first or use global max/min to ensure scale is global.
          Volume_Score = scales::rescale(Log_Volume, to = c(0.01, 1)),

          # Composite Importation Risk Score — all five components are [0,1].
          # The multiplicative structure reflects serial epidemiological dependency:
          # each factor conditions the next (travel → reservoir → spillover →
          # undetected → importation). If any link is near zero, overall risk is near zero.
          Exposure_Score = Volume_Score * Reservoir_Risk * Spillover_Activity * Detection_Capacity * Seasonal_Index
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

eps <- 1e-6

risk_stats <- df_analyzed |>
     filter(Exposure_Score > 0) |>
     mutate(
          Safe_Score = pmin(pmax(Exposure_Score, eps), 1 - eps)
     ) |>
     summarise(
          Logit_Mean = mean(log(Safe_Score / (1 - Safe_Score))),
          Logit_SD = sd(log(Safe_Score / (1 - Safe_Score)))
     )

spc_base <- risk_stats$Logit_Mean
spc_step <- risk_stats$Logit_SD

df_analyzed <- df_analyzed |>
     mutate(Safe_Score = pmin(pmax(Exposure_Score, eps), 1 - eps),
            Logit_Score = log(Safe_Score / (1 - Safe_Score)),
            
            Risk_Level_SPC = case_when(Logit_Score > (spc_base + 3 * spc_step) ~ "Level 5: Critical",   # > 3SD (极端罕见)
                                       Logit_Score > (spc_base + 2 * spc_step) ~ "Level 4: High",       # > 2SD (显著异常)
                                       Logit_Score > (spc_base + 1 * spc_step) ~ "Level 3: Moderate",   # > 1SD (需关注)
                                       Logit_Score > (spc_base + 0.5 * spc_step) ~ "Level 2: Low",      # > 0.5SD (轻微偏高)
                                       TRUE ~ "Level 1: Very Low")) |>
     select(-Safe_Score) |> 
     mutate(Risk_Level_Q = case_when(Exposure_Score > 0.95 ~ "Level 5: Critical",
                                     Exposure_Score > 0.90 ~ "Level 4: High",
                                     Exposure_Score > 0.75 ~ "Level 3: Moderate",
                                     Exposure_Score > 0.50 ~ "Level 2: Low",
                                     TRUE ~ "Level 1: Very Low"))

df_analyzed$Risk_Level <- factor(df_analyzed$Risk_Level_SPC, levels = c(
     "Level 1: Very Low", "Level 2: Low", "Level 3: Moderate", "Level 4: High", "Level 5: Critical"
))

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
     # Let ggplot compute per-facet breaks (scales = "free_y") but format labels in millions
     scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 0.01),
                            limits = c(0, NA),
                            expand = expansion(mult = c(0, 0.15))) +
     scale_x_date(date_labels = "%b %Y",
                  expand = c(0, 0),
                  date_breaks = "6 months") +
     geom_vline(xintercept = as.numeric(last_date), linetype = "dashed", color = "red") +
     labs(
          title = 'A',
          x = NULL, 
              y = "Monthly volume (millions)",
          color = "Risk Level"
     ) +
     theme_bw() +
     theme(legend.position = "bottom",
           plot.title.position = "plot")

# Plot Risk Forecast Bar (Details)
df_forecast_only <- df_plot_risk |> filter(Type == "Forecast")

# Use actual Date_Obj on x and format as dates so ordering is chronological;
# rotate labels to avoid overlap and show monthly ticks for the forecast horizon.
p_risk_bar <- ggplot(df_forecast_only, aes(x = Date_Obj, y = Exposure_Score, fill = Risk_Level)) +
     geom_col(show.legend = F) +
     facet_wrap(~Origin_Country, scales = "free_y", ncol = 1) +
     # Keep full color scale available for all risk levels
     scale_fill_brewer(palette = "RdYlGn", direction = -1, drop = FALSE) +
     # Remove scientific notation; use simple decimal for the [0,1] score
     scale_y_continuous(labels = scales::number_format(accuracy = 0.0001),
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.2))) +
     scale_x_date(date_labels = "%b %Y",
                  date_breaks = "1 month",
                  expand = c(0, 0)) +
     labs(
          title = "B",
          x = "Date", y = "Importation risk score (SPC)"
     ) +
     theme_bw() +
     theme(
          legend.position = "none",
          plot.title.position = "plot"
     )

df_forecast_only$Risk_Level_Q <- factor(df_forecast_only$Risk_Level_Q, levels = c(
     "Level 1: Very Low", "Level 2: Low", "Level 3: Moderate", "Level 4: High", "Level 5: Critical"
))
p_risk_bar_q <- ggplot(df_forecast_only, aes(x = Date_Obj, y = Exposure_Score, fill = Risk_Level_Q)) +
     geom_col(show.legend = F) +
     facet_wrap(~Origin_Country, scales = "free_y", ncol = 1) +
     scale_fill_brewer(palette = "RdYlGn", direction = -1, drop = FALSE) +
     scale_y_continuous(labels = scales::number_format(accuracy = 0.0001),
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.2))) +
     scale_x_date(date_labels = "%b %Y",
                  date_breaks = "1 month",
                  expand = c(0, 0)) +
     labs(
          title = "C",
          x = "Date", y = "Importation risk score (Quantile)"
     ) +
     theme_bw() +
     theme(
          legend.position = "none",
          plot.title.position = "plot"
     )

# Compose Figure 2
fig2 <- p_risk_trend + p_risk_bar + p_risk_bar_q + plot_layout(widths = c(2, 1, 1), guides = 'collect') &
     theme(legend.position = 'bottom')

ggsave("./Outcome/fig4.png",
       fig2,
       width = 20, 
       height = 10,
       bg="white")
