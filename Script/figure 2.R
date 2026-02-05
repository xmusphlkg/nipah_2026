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
library(sf) # Use sf for shapefile handling
sf_use_s2(FALSE) # Disable S2 to avoid topology errors on centroids

# 1. Load Data (5 Countries + Global Context)
data_folder <- "./Data/OAG/"
files <- list.files(data_folder, pattern = "\\.xlsx$", full.names = TRUE)
file_global_hist <- "./Data/全球入境流量统计.xlsx"

# A. Read 5 Countries Data
read_country_data <- function(filepath) {
     df <- read.xlsx(filepath)
     # Columns: 出发国家 (Origin), 到达国家 (Dest), 旅客数 (Volume), 时间 (Time: YYYYMM)
     names(df) <- c("Origin_Country", "Dest_Country", "Inbound_Volume", "Time")
     
     df <- df |>
          mutate(
               Time_Str = as.character(Time),
               Year = as.numeric(substr(Time_Str, 1, 4)),
               Month = as.numeric(substr(Time_Str, 5, 6)),
               Date_Obj = as.Date(paste(Year, Month, "01", sep = "-"))
          ) |>
          select(Origin_Country, Dest_Country, Inbound_Volume, Date_Obj, Year, Month)
     
     return(df)
}

df_countries <- map_dfr(files, read_country_data)

# B. Read Global Inbound Data (For context/comparison)
df_global_raw <- read.xlsx(file_global_hist)
# Inspect structure: usually Time, Global
names(df_global_raw) <- c("Time", "Global_Inbound_Volume")

df_global <- df_global_raw |>
     select(Time, Global_Inbound_Volume) |>
     mutate(
          Time_Str = as.character(Time),
          Year = as.numeric(substr(Time_Str, 1, 4)),
          Month = as.numeric(substr(Time_Str, 5, 6)),
          Date_Obj = as.Date(paste(Year, Month, "01", sep = "-"))
     ) |>
     select(Date_Obj, Global_Inbound_Volume)


# 1.1 Forecast (Holt-Winters)

# Determine global last date to set forecast start
last_date <- max(df_countries$Date_Obj)
forecast_horizon <- 6 

# --- Function to Forecast Single Series ---
forecast_series <- function(ts_data, name_label) {
     model_hw <- tryCatch({
          hw(ts_data, h = forecast_horizon, seasonal = "multiplicative")
     }, error = function(e) {
          hw(ts_data, h = forecast_horizon, seasonal = "additive")
     })
     
     future_dates <- seq(last_date, by = "month", length.out = forecast_horizon + 1)[-1]
     
     data.frame(
          Date_Obj = future_dates,
          Forecast_Value = as.numeric(model_hw$mean),
          Type = "Forecast"
     )
}

# Forecast 5 Countries
get_country_forecast <- function(df_subset) {
     cntry <- df_subset$Origin_Country[1]
     df_subset <- df_subset |> arrange(Date_Obj)
     
     start_y <- df_subset$Year[1]
     start_m <- df_subset$Month[1]
     ts_vol <- ts(df_subset$Inbound_Volume, start = c(start_y, start_m), frequency = 12)
     
     df_fc <- forecast_series(ts_vol, cntry) |>
          mutate(
               Origin_Country = cntry, 
               Dest_Country = "China",
               Inbound_Volume = Forecast_Value
          ) |> select(-Forecast_Value)
     
     # Add Year/Month
     df_fc <- df_fc |>
          mutate(
               Year = as.numeric(format(Date_Obj, "%Y")),
               Month = as.numeric(format(Date_Obj, "%m"))
          )
     
     bind_rows(df_subset |> mutate(Type="History"), df_fc)
}

df_countries_all <- df_countries |>
     group_split(Origin_Country) |>
     map_dfr(get_country_forecast)

# Forecast Global Volume
ts_global <- ts(df_global$Global_Inbound_Volume, start = c(2016, 1), frequency = 12) # Approximation of start, adjusting later if needed
# Better to extract start from data
ts_global <- ts(df_global$Global_Inbound_Volume, 
                start = c(as.numeric(format(min(df_global$Date_Obj), "%Y")), 
                          as.numeric(format(min(df_global$Date_Obj), "%m"))), 
                frequency = 12)

df_global_fc <- forecast_series(ts_global, "Global") |>
     rename(Global_Inbound_Volume_Fc = Forecast_Value)

df_global_combined <- bind_rows(
     df_global |> mutate(Type = "History"),
     df_global_fc |> rename(Global_Inbound_Volume = Global_Inbound_Volume_Fc)
) |> select(Date_Obj, Global_Inbound_Volume, Type)

# Merge Global context into Countries
df_final <- df_countries_all |>
     left_join(df_global_combined, by = c("Date_Obj", "Type"))

# 2. Risk Analysis (Improved Probabilistic Exposure Model)
# ==============================================================================
# Methodology: Multi-Factorial Importation Risk Model
# ------------------------------------------------------------------------------
# We calculate biological risk based on three core components:
#
# 1. Volume (Inbound_Volume):
#    - The sheer number of travelers from the origin country.
#
# 2. Epidemiological Context (Risk_Metadata):
#    - Reservoir_Risk: Evidence of virus in natural hosts (Pteropus bats).
#    - Spillover_Activity: Rate of active transmission events.
#    - Endemic_Connectivity: The probability that a traveler comes from the specific 
#      *endemic region* within the country. 
#      (Crucial for India: Outbreaks are in Kerala (remote), while flights are from Delhi/Mumbai.
#       This factor dilutes the risk for large countries with localized outbreaks.)
#
# 3. Seasonal Susceptibility (Seasonal_Index):
#    - Biological seasonality peaking in Winter/Spring (Feb).
#
# Formula: 
# Exposure_Score = Inbound_Volume * Reservoir_Risk * Spillover_Activity * Endemic_Connectivity * Seasonal_Index
# ==============================================================================

# Define Epidemiological Parameters
risk_metadata <- tibble::tribble(
  ~Origin_Country, ~Reservoir_Risk, ~Spillover_Activity, ~Endemic_Connectivity, ~Notes,
  "India",         5.0,             1.0,                 0.1, "Localized outbreaks (Kerala) far from China-bound hubs (High dilution)",
  "Bangladesh",    5.0,             1.0,                 0.8, "Small geography, outbreaks near transport hubs (Low dilution)",
  "Malaysia",      2.0,             0.05,                0.5, "Good connectivity, but reservoir spillover path blocked",
  "Philippines",   2.0,             0.1,                 0.2, "Archipelago limit: outbreaks often on isolated islands",
  "Singapore",     0.5,             0.01,                1.0, "City-state, perfect connectivity but no reservoir risk"
)

df_analyzed <- df_final |>
  left_join(risk_metadata, by = "Origin_Country") |>
  mutate(
    # Handle missing countries (Default to low risk if unknown)
    Reservoir_Risk = replace_na(Reservoir_Risk, 1.0),
    Spillover_Activity = replace_na(Spillover_Activity, 0.01),
    Endemic_Connectivity = replace_na(Endemic_Connectivity, 0.1),
    
    # A. Seasonal Factor: Smooth curve peaking in Feb (Winter/Spring) for Nipah
    # Using Cosine function to simulate biological seasonality (Range: 1.0 to 3.0)
    # Peak at Month 2 (Feb), Trough at Month 8 (Aug)
    Seasonal_Index = 1 + 2 * (0.5 * (cos((Month - 2) / 6 * pi) + 1)),
    
    # B. Calculated Composite Risk Score
    # Added Endemic_Connectivity to account for geographic dilution
    Exposure_Score = Inbound_Volume * Reservoir_Risk * Spillover_Activity * Endemic_Connectivity * Seasonal_Index
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

# 4. Visualization 1: Volume, Share and Map
################################

# A. Global Trend vs 5-Countries Trend
# Prepare Global + 5-Country Sum
df_5_sum <- df_analyzed |>
     group_by(Date_Obj, Type) |>
     summarise(Volume_5_Countries = sum(Inbound_Volume), .groups = 'drop')

df_plot_A <- df_global_combined |>
     left_join(df_5_sum, by = c("Date_Obj", "Type")) |>
     pivot_longer(cols = c(Global_Inbound_Volume, Volume_5_Countries),
                  names_to = "Category", values_to = "Volume") |>
     mutate(Category = factor(Category, levels = c("Global_Inbound_Volume", "Volume_5_Countries"),
                              labels = c("Global", "5 Countries")))

df_plot_B <- df_plot_A |> 
     pivot_wider(names_from = Category, values_from = Volume) |> 
     mutate(Share_5_Countries = `5 Countries` / Global)

p_global <- ggplot(df_plot_A, aes(x = Date_Obj, y = Volume, color = Category)) +
     # Forecast Background
     annotate("rect", xmin = as.Date(last_date), xmax = max(df_plot_A$Date_Obj),
              ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
     geom_vline(xintercept = as.numeric(last_date), linetype = "dashed", col = "red") +
     geom_line() +
     scale_color_manual(values = c("Global" = "black", "5 Countries" = "#119DA4FF")) +
     scale_y_continuous(labels = scales::scientific,
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.15))) +
     scale_x_date(date_labels = "%Y",
                  expand = c(0, 0),
                  date_breaks = "1 year") +
     labs(title = "A. Inbound Volume", x = NULL, y = "Monthly Volume")+
     theme_bw() +
     theme(legend.position = c(0.4, 0.99),
           legend.justification = c(0, 1),
           legend.title = element_blank(), 
           legend.direction = "horizontal",
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           legend.background = element_rect(fill="white", color="gray90"))

p_global_bar <- ggplot(df_plot_B, aes(x = Date_Obj, y = Share_5_Countries)) +
     annotate("rect", xmin = as.Date(last_date), xmax = max(df_plot_B$Date_Obj),
              ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
     geom_vline(xintercept = as.numeric(last_date), linetype = "dashed", col = "red") +
     geom_area(fill = "#119DA4FF") +
     scale_y_continuous(labels = scales::percent,
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.15))) +
     scale_x_date(date_labels = "%Y",
                  expand = c(0, 0),
                  date_breaks = "1 year") +
     labs(x = 'Date', y = "Proportion")+
     theme_bw()

p_global <- p_global + p_global_bar + 
     plot_layout(ncol = 1, heights = c(5, 1))

# B. 5 Countries Trend
p_countries <- ggplot(df_analyzed, aes(x = Date_Obj, y = Inbound_Volume, color = Origin_Country)) +
     annotate("rect", xmin = as.Date(last_date), xmax = max(df_plot_A$Date_Obj),
              ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
     geom_line() +
     geom_vline(xintercept = as.numeric(last_date), linetype = "dashed", col = "red") +
     scale_y_continuous(labels = scales::scientific,
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.15))) +
     scale_x_date(date_labels = "%Y",
                  expand = c(0, 0),
                  date_breaks = "1 year") +
     scale_color_manual(name = "Country",
                        values = c("#462255FF", "#FF8811FF", "#9DD9D2FF", "#046E8FFF", "#D44D5CFF")) +
     labs(title = "C. Individual trend", x = 'Date', y = "Volume", color = "Country") +
     theme_bw() +
     theme(legend.position = "none")

# C. Yearly Share Stacked Bar
df_yearly_share <- df_analyzed |>
     group_by(Year, Origin_Country) |>
     summarise(Avg_Share = mean(Share_of_Global, na.rm = TRUE), .groups = "drop")

p_share <- ggplot(df_yearly_share, aes(x = as.character(Year), y = Avg_Share, fill = Origin_Country)) +
     geom_col(position = "fill") +
     scale_fill_manual(name = "Country",
                       values = c("#462255FF", "#FF8811FF", "#9DD9D2FF", "#046E8FFF", "#D44D5CFF")) +
     scale_y_continuous(labels = percent,
                        limits = c(0, 1),
                        expand = expansion(mult = c(0, 0))) +
     labs(x = "Year", y = "Share of global volume", title = "D. Share of global volume") +
     theme_bw() +
     theme(legend.position = "bottom")

# D. Map with Bubbles (Using globalmap.shp)
# Load Shapefile
shp_path <- "./Data/globalmap/globalmap.shp"
world_shp <- st_read(shp_path, quiet = TRUE)

# Normalize names for matching
world_shp <- world_shp |>
     mutate(NAME_UPPER = toupper(NAME_ENG))

# Prepare Centroid Data
# We match our 5 countries + China
# Data names: India, Bangladesh, Singapore, Philippines, Malaysia
# Shapefile names: Need to ensure they match. 
# "PHILIPPINES", "MALAYSIA", "SINGAPORE", "BANGLADESH", "INDIA", "CHINA" usually exist.

# Extract China Geometry for highlighting
china_sf <- world_shp |> filter(NAME_UPPER == "CHINA" | NAME_CHN == "中国")

# Extract Source Countries Geometries & Calculate Centroids
source_names <- toupper(unique(df_analyzed$Origin_Country))
source_countries_sf <- world_shp |> 
     filter(NAME_UPPER %in% source_names)

# Calculate centroids
# Use st_centroid on spherical geometry might warn, but fine for mapping
source_centroids <- st_centroid(source_countries_sf)
source_coords <- st_coordinates(source_centroids)

# Combine with name for joining
source_centroids <- source_centroids |>
     bind_cols(as.data.frame(source_coords))

# Prepare Volume Data for bubbles
recent_vol_agg <- df_analyzed |>
     filter(Type == "History", Year >= (max(Year)-1)) |>
     group_by(Origin_Country) |>
     summarise(Volume = mean(Inbound_Volume)) |>
     mutate(NAME_UPPER = toupper(Origin_Country))

# Join
map_points <- source_centroids |>
     left_join(recent_vol_agg, by = "NAME_UPPER")

p_map <- ggplot() +
     # Base Map (All countries)
     geom_sf(data = world_shp, fill = "cornsilk", color = "white", linewidth = 0.1) +
     # Highlight China
     geom_sf(data = china_sf, fill = "gray80", color = "white", linewidth = 0.2) +
     # Source Country Bubbles
     geom_point(data = map_points, aes(x = X, y = Y, size = Volume),
                color = "red", alpha = 0.6) +
     scale_size_continuous(range = c(3, 12), name = "Volume", labels = comma) +
     # Labels
     geom_text(data = map_points, aes(x = X, y = Y, label = Origin_Country),
               size = 3, vjust = -1.5, fontface = "bold", color = "black") +
     # Zoom to relevant area (Asia-Pacific)
     coord_sf(xlim = c(60, 150), ylim = c(-10, 50)) +
     labs(title = "B. Source Volume Distribution") +
     theme_bw() +
     theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
     )

# Compose Figure 1
fig1 <- free(p_global) + p_countries + free(p_map, side = 'l') + free(p_share) +
     plot_layout(ncol = 2, widths = c(1, 1)) 

ggsave("./Outcome/fig2.png",
       fig1,
       width = 14, height = 10, bg="white")


# 5. Visualization 2: Risk Assessment
#####################################

plot_start_year <- 2024
df_plot_risk <- df_analyzed |> filter(Year >= plot_start_year)

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
     # Force display of unused expected levels (Level 4/5) in logic
     scale_color_manual(values = fill_colors, drop = FALSE) +
     scale_y_continuous(labels = scales::comma) +
     geom_vline(xintercept = as.numeric(last_date), linetype = "dashed", color = "red") +
     labs(
          title = 'A',
          x = NULL, 
          y = "Inbound Passenger Volume",
          color = "Risk Level"
     ) +
     theme_bw() +
     theme(legend.position = "bottom")

# Plot Risk Forecast Bar (Details)
df_forecast_only <- df_plot_risk |> filter(Type == "Forecast")

p_risk_bar <- ggplot(df_forecast_only, aes(x = as.factor(format(Date_Obj, "%Y-%m")), y = Exposure_Score, fill = Risk_Level)) +
     geom_col(show.legend = F) +
     facet_wrap(~Origin_Country, scales = "free_y", ncol = 1) +
     # Force display of unused expected levels (Level 4/5) in logic
     scale_fill_brewer(palette = "RdYlGn", direction = -1, drop = FALSE) +
     labs(
          title = "B",
          x = "Forecast Month", y = "Composite Risk Score"
     ) +
     theme_bw() +
     theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none" # Shared legend implied or redundant
     )

# Compose Figure 2
fig2 <- p_risk_trend + p_risk_bar + plot_layout(widths = c(2, 1), guides = 'collect') &
     theme(legend.position = 'bottom')

ggsave("./Outcome/fig3.png",
       fig2,
       width = 16, 
       height = 10, 
       bg="white")

