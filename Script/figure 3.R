
#####################################
## @Description: Generate Figure 3 (component risks: reservoir,
##               spillover, detection, and seasonality) and export
##               fig3_risk_metadata.csv for the importation model.
## @version: 1.0
## @Author: Li Kangguo
## @Date: 2026-02-04 15:04:57
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-02-28 12:00:00
#####################################

library(tidyverse)
library(openxlsx)
library(patchwork)

# panel A -----------------------------------------------------------------

pA <- ggplot() +
     xlim(0,1) + ylim(0,1) +
     theme_void() +
     annotate("text",
              x = 0.5, y = 0.6,
              label = "Risk[i,t] == Volume[i,t] %*% (R[reservoir] %*% A[spillover] %*% D[detact]) %*% S[season,t]",
              parse = TRUE,
              size = 8) +
     labs(title = "A")

# panel B -----------------------------------------------------------------

data_panel_B <- read.csv("./Outcome/risk_assessment.csv")

pB_1 <- ggplot(data_panel_B,
               aes(y = Country, x = GBIF_Recs)) +
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
     labs(x = "Host occurrences (GBIF)", y = NULL, title = 'A') +
     theme_bw()

pB_2 <- ggplot(data_panel_B,
               aes(y = Country, x = NCBI_Seqs)) +
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
     labs(x = "Viral sequences (NCBI)", y = NULL) +
     theme_bw()+
     theme(axis.text.y = element_blank())

pB_3 <- ggplot(data_panel_B,
               aes(y = Country, x = Reservoir_Risk_Calculated)) + 
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        breaks = seq(0, 1, by = 0.2),
                        limits = c(0, 1.05)) +
     labs(x = "Composite reservoir risk (0-1)", y = NULL) +
     theme_bw()+
     theme(axis.text.y = element_blank())

pB <- pB_1 + pB_2 + pB_3 + plot_layout(nrow = 1, widths = c(1,1,1)) &
     theme(plot.title.position = "plot")
     
# panel C -----------------------------------------------------------------

data_panel_C <- read.csv("./Data/Risk/data.csv")

data_panel_C <- data_panel_C |>
     mutate(
          # Normalize special characters: 
          # Remove punctuation spaces (U+2008) from Area and fix middle-dot decimals in Pop
          Area_raw = str_replace_all(Area.at.risk..km2., "[ \\s]", ""),
          Pop_raw  = str_replace_all(Population.at.risk..millions., "·", "."),
          
          # Identify Region headers (rows where Area is empty) to fill downwards
          Region_Label = ifelse(Area.at.risk..km2. == "" | is.na(Area.at.risk..km2.), Empty.Cell, NA)
     ) |>
     fill(Region_Label, .direction = "down") |>
     # Remove the header/summary rows
     filter(Area.at.risk..km2. != "" & !is.na(Area.at.risk..km2.)) |>
     
     # 2. Extract numeric values for Estimates and CIs
     mutate(
          Country = str_remove_all(Country, "[*†]") |> trimws(),
          
          # --- Area Statistics ---
          Area_km2     = str_extract(Area_raw, "^[0-9.]+") |> as.numeric(),
          Area_CI_Low  = str_extract(Area_raw, "(?<=\\()[0-9.]+") |> as.numeric(),
          Area_CI_High = str_extract(Area_raw, "(?<=–)[0-9.]+(?=\\))") |> as.numeric(),
          
          # --- Population Statistics ---
          Pop_millions = str_extract(Pop_raw, "^[0-9.]+") |> as.numeric(),
          Pop_CI_Low   = str_extract(Pop_raw, "(?<=\\()[0-9.]+") |> as.numeric(),
          Pop_CI_High  = str_extract(Pop_raw, "(?<=–)[0-9.]+(?=\\))") |> as.numeric()
     ) |>
     # 3. Final selection for a tidy format
     select(
          Region = Region_Label, 
          Country, 
          Area_km2, Area_CI_Low, Area_CI_High,
          Pop_millions, Pop_CI_Low, Pop_CI_High
     ) |> 
     filter(Country %in% data_panel_B$Country) |> 
     mutate(
          spillover = scales::rescale(Pop_millions, to = c(0.01, 1))
     )

pC_1 <- ggplot(data_panel_C,
                  aes(y = Country))+
     geom_point(aes(x = Pop_millions), color = "#119DA4FF", size = 3) +
     geom_errorbarh(aes(xmin = Pop_CI_Low, xmax = Pop_CI_High), height = 0.2, color = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 120)) +
     labs(x = "Population at risk (millions)", y = NULL, title = 'B') +
     theme_bw()+
     theme(plot.title.position = "plot")

pC_2 <- ggplot(data_panel_C,
                  aes(y = Country))+
     geom_col(aes(x = spillover), fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        breaks = seq(0, 1, by = 0.2),
                        limits = c(0, 1)) +
     labs(x = "Normalized spillover risk", y = NULL) +
     theme_bw()+
     theme(axis.text.y = element_blank(),
           plot.margin = margin(l = 20))

pC <- pC_1 + pC_2 + plot_layout(nrow = 1, widths = c(1,1))

# panel D -----------------------------------------------------------------

data_panel_D <- read.xlsx("./Data/GHS/2021.xlsx")

data_panel_D <- data_panel_D |>
     select(Country, Detect) |>
     mutate(
          # Normalize GHS score (0-100) to 0-1, then invert so that high capacity = low risk.
          Detect_Normalized = scales::rescale(Detect, to = c(0, 0.99)),
          Detect_score = 1 - Detect_Normalized,
          # Ensure strict lower bound to prevent pure zero in multiplicative model
          Detect_score = pmax(Detect_score, 0.01)
     )

pD_1 <- ggplot(data_panel_D,
                  aes(y = Country, x = Detect))+
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 100)) +
     labs(x = "Estimated detection capacity", y = NULL, title = 'C') +
     theme_bw()+
     theme(plot.title.position = "plot")

pD_2 <- ggplot(data_panel_D,
                  aes(y = Country, x = Detect_score))+
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        breaks = seq(0, 1, by = 0.2),
                        limits = c(0, 1)) +
     labs(x = "Normalized detection risk", y = NULL) +
     theme_bw()+
     theme(axis.text.y = element_blank(),
           plot.margin = margin(l = 20))

pD <- pD_1 + pD_2 + plot_layout(nrow = 1, widths = c(1,1))

# panel E -----------------------------------------------------------------

load("./Outcome/seasonal_model.RData")

# Read raw outbreak records and aggregate to Year-Month observations
library(lubridate)
outbreaks_raw <- read.csv("./Data/Outbreak/outbreaks_nipah.csv", stringsAsFactors = FALSE)

obs_monthly <- outbreaks_raw |>
     mutate(start_parsed = str_extract(started, "^[^-]+"),
            Date = mdy(str_trim(start_parsed)),
            cases = as.numeric(cases)) |>
     filter(!is.na(Date)) |>
     mutate(Year = year(Date), Month = month(Date)) |>
     group_by(Year, Month) |>
     summarise(Cases = sum(cases, na.rm = TRUE), .groups = 'drop')

# Extract bootstrap summaries from seasonal_model
season_pred <- tibble(
     Month = 1:12,
     Predicted = seasonal_model$pred_mean,
     Lower = seasonal_model$pred_lower,
     Upper = seasonal_model$pred_upper
)

# Plot: show raw Year-Month observations as jittered points + fitted line and CI
pE <- ggplot() +
     geom_jitter(data = obs_monthly, aes(x = Month, y = Cases, color = 'Observed'), width = 0.15, height = 0, size = 2, alpha = 0.8) +
     geom_ribbon(data = season_pred, aes(x = Month, ymin = Lower, ymax = Upper), fill = "#D44D5CFF", alpha = 0.2) +
     geom_line(data = season_pred, aes(x = Month, y = Predicted, color = 'Fitted'), size = 1) +
     scale_color_manual(values = c('Observed' = "#046E8FFF", 'Fitted' = "#D44D5CFF")) +
     scale_x_continuous(breaks = 1:12, labels = month.abb) +
     labs(x = "Month", y = "Monthly cases", title = 'D') +
     theme_bw() +
     theme(legend.position = 'inside',
           legend.direction = "horizontal",
           legend.title = element_blank(),
           legend.position.inside = c(0.5, 0.99),
           legend.justification = c(0.5, 1))
     
# Combine all panels ------------------------------------------------------

fig <- cowplot::plot_grid(pB, pC,
                          pD, pE,
                          rel_widths = c(1, 0.6),
                          ncol = 2, align = "v")

ggsave("./Outcome/fig3.png",
       fig,
       width = 16, 
       height = 7, 
       bg="white")

risk_metadata <- data_panel_B |> 
     select(Country, Reservoir_Risk = Reservoir_Risk_Calculated) |> 
     left_join(data_panel_C |> select(Country, Spillover_Activity = spillover), by = "Country") |> 
     left_join(data_panel_D |> select(Country, Detection_Capacity = Detect_score), by = "Country")

write.csv(risk_metadata,
          "./Outcome/fig3_risk_metadata.csv",
          row.names = F)
