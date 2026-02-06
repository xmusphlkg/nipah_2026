
library(tidyverse)
library(openxlsx)

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
               aes(y = Country, x = Score_Host)) +
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 4)) +
     labs(x = "Host risk score", y = NULL, title = 'B') +
     theme_bw()

pB_2 <- ggplot(data_panel_B,
               aes(y = Country, x = Score_Virus)) +
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 16)) +
     labs(x = "Virus risk score", y = NULL) +
     theme_bw()+
     theme(axis.text.y = element_blank())

pB_3 <- ggplot(data_panel_B,
               aes(y = Country, x = Reservoir_Risk_Calculated)) + 
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 5)) +
     labs(x = "Reservoir risk score", y = NULL) +
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
     labs(x = "Population at risk (millions)", y = NULL, title = 'C') +
     theme_bw()+
     theme(plot.title.position = "plot")

pC_2 <- ggplot(data_panel_C,
                  aes(y = Country))+
     geom_col(aes(x = spillover), fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 1)) +
     labs(x = "Normalized spillover risk", y = NULL) +
     theme_bw()+
     theme(axis.text.y = element_blank())

pC <- pC_1 + pC_2 + plot_layout(nrow = 1, widths = c(1,1))

# panel D -----------------------------------------------------------------

data_panel_D <- read.xlsx("./Data/GHS/2021.xlsx")

data_panel_D <- data_panel_D |>
     select(Country, Detect) |> 
     mutate(Detect_score = scales::rescale(-Detect, to = c(0.01, 1)))

pD_1 <- ggplot(data_panel_D,
                  aes(y = Country, x = Detect))+
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 100)) +
     labs(x = "Estimated detection capacity", y = NULL, title = 'D') +
     theme_bw()+
     theme(plot.title.position = "plot")

pD_2 <- ggplot(data_panel_D,
                  aes(y = Country, x = Detect_score))+
     geom_col(fill = "#119DA4FF") +
     scale_x_continuous(expand = expansion(mult = c(0, 0)),
                        limits = c(0, 1)) +
     labs(x = "Normalized detection risk", y = NULL) +
     theme_bw()+
     theme(axis.text.y = element_blank())

pD <- pD_1 + pD_2 + plot_layout(nrow = 1, widths = c(1,1))

# panel E -----------------------------------------------------------------

load("./Outcome/seasonal_model.RData")

print(summary(harmonic_fit))

# Plot fitted seasonal curve vs observed monthly averages
season_pred <- tibble(Month = 1:12,
                      Predicted = seasonal_predict(1:12))

# Combine with monthly statistics for plotting
plot_season_df <- monthly_avg |>
     left_join(season_pred, by = "Month") |>
     mutate(MonthLabel = factor(month.abb[Month], levels = month.abb))

pE <- ggplot(plot_season_df, aes(x = Month)) +
     geom_point(aes(y = AvgCases, color = 'Average'), size = 3) +
     geom_line(aes(y = Predicted, color = 'Fitted'), size = 1) +
     scale_color_manual(values = c('Average' = "#046E8FFF", 'Fitted' = "#D44D5CFF")) +
     scale_x_continuous(breaks = 1:12, labels = month.abb) +
     labs(x = "Month", y = "Average monthly cases", title = 'E') +
     theme_bw() +
     theme(legend.position = 'inside',
           legend.direction = "horizontal",
           legend.title = element_blank(),
           legend.position.inside = c(0.5, 0.99),
           legend.justification = c(0.5, 1))
     
# Combine all panels ------------------------------------------------------

design <- "
AF
BD
CE
"

fig <- pA + pB + pC + pD + free(pE) + patchwork::plot_spacer() +
     plot_layout(design = design, widths = c(1, 0.6), heights = c(0.15, 1, 1)) &
     theme(plot.title.position = "plot")

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
