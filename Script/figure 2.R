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

# load data ---------------------------------------------------------------

data_folder <- "./Data/OAG/"
files <- list.files(data_folder, pattern = "\\.xlsx$", full.names = TRUE)
file_global_hist <- "./Data/全球入境流量统计.xlsx"

# A. Read 6 Countries Data
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

# Forecast (Holt-Winters) -------------------------------------------------

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

rm(df_countries, df_global, df_global_fc, df_countries_all)

# visual ------------------------------------------------------------------

fill_color <- c("#462255FF", "#FF8811FF", "#9DD9D2FF", "#046E8FFF", "#D44D5CFF", "grey50")
names(fill_color) <- c("Bangladesh", "India", "Malaysia", "Singapore", "Philippines", "Thailand")

# A. Global Trend vs 5-Countries Trend
# Prepare Global + 5-Country Sum
df_5_sum <- df_final |>
     group_by(Date_Obj, Type) |>
     summarise(Volume_5_Countries = sum(Inbound_Volume), .groups = 'drop')

df_plot_A <- df_global_combined |>
     left_join(df_5_sum, by = c("Date_Obj", "Type")) |>
     pivot_longer(cols = c(Global_Inbound_Volume, Volume_5_Countries),
                  names_to = "Category", values_to = "Volume") |>
     mutate(Category = factor(Category, levels = c("Global_Inbound_Volume", "Volume_5_Countries"),
                              labels = c("Global", "6 Countries")))

df_plot_B <- df_plot_A |> 
     pivot_wider(names_from = Category, values_from = Volume) |> 
     mutate(Share_5_Countries = `6 Countries` / Global)

p_global <- ggplot(df_plot_A, aes(x = Date_Obj, y = Volume, color = Category)) +
     # Forecast Background
     annotate("rect", xmin = as.Date(last_date), xmax = max(df_plot_A$Date_Obj),
              ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
     geom_vline(xintercept = as.numeric(last_date), linetype = "dashed", col = "red") +
     geom_line() +
     scale_color_manual(values = c("Global" = "black", "6 Countries" = "#119DA4FF")) +
     scale_y_continuous(labels = scales::scientific,
                        limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.15))) +
     scale_x_date(date_labels = "%Y",
                  expand = c(0, 0),
                  date_breaks = "1 year") +
     labs(title = "A", x = NULL, y = "Monthly Volume")+
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
p_countries <- ggplot(df_final, aes(x = Date_Obj, y = Inbound_Volume, color = Origin_Country)) +
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
                        values = fill_color)+
     labs(title = "B", x = 'Date', y = "Volume", color = "Country") +
     theme_bw() +
     theme(legend.position = "bottom") +
     guides(color = guide_legend(nrow = 1, byrow = TRUE))

# Compose Figure 1
fig1 <- free(p_global) | p_countries

ggsave("./Outcome/fig2.png",
       fig1,
       width = 14, height = 6, bg="white")

save.image("./Outcome/figure2_workspace.RData")
