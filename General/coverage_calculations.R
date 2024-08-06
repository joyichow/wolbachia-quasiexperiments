# Calculating average coverage across time
library(readxl)
library(dplyr)
library(lubridate)

# Malaysia 
release_num <- "/your/local/file/path/wolb_interventions/Release Numbers - Frequency of Wolbachia - Ovitrap Index.xlsx" # change to your local file path

release_sheets <- excel_sheets(release_num)

release_data <- lapply(release_sheets, function(sheet_name) {
  site_data <- read_excel(release_num, sheet = sheet_name)
  site_data <- mutate(site_data,
                      Date = as.Date(paste0("20", substr(`Year/Week`, 1,2), "W", substr(`Year/Week`, 4, 5)), format = "%Y W%U"),
                      Site = sheet_name)
  site_data <- mutate(site_data, Months_since_start = interval(min(Date), Date) / months(1))
  
  site_data <- site_data %>%
    select(YearWeek = 1, Continuous_FreqWolb_for_graph = 14, Date, Months_since_start, Site)
  
  return(site_data)
})

combined_release <- bind_rows(release_data)

release_aggregated <- combined_release  %>%
  arrange(YearWeek) %>%
  group_by(YearWeek, Site) %>%
  summarize(FreqWolb_Avg = mean(Continuous_FreqWolb_for_graph, na.rm = TRUE)) 

release_wide <- release_aggregated %>% 
  pivot_wider(
    names_from = Site, 
    values_from = FreqWolb_Avg,
    values_fill = list(FreqWolb_Avg = NA) # Fill missing values with NA
  )

colnames(release_wide) <- c("YearWeek", "AL", "SF", "SL", "AF", "MH", "SC")

intervention_start <- data.frame(
  Site = c("MH", "SF", "SL", "SC", "AL", "AF"),
  StartYearWeek = c("17/31", "17/18", "17/18", "17/46", "17/12", "17/36")
)

combined_release <- combined_release %>%
  left_join(intervention_start, by = "Site")

intervention_start <- intervention_start %>%
  mutate(
    Year = as.numeric(substr(StartYearWeek, 1, 2)) + 2000, 
    Week = as.numeric(substr(StartYearWeek, 4, 5)),        
    StartDate = as.Date(paste(Year, Week, 1, sep = "-"), format = "%Y-%U-%u") 
  ) %>%
  select(-Year, -Week)

release_long <- release_wide %>%
  gather(Site, Coverage, -YearWeek) %>%
  mutate(YearWeek = as.Date(paste0("20", substr(YearWeek, 1, 2), "-W", substr(YearWeek, 4, 5), "-1"), format = "%Y-W%U-%u"))

release_merged <- release_long %>%
  left_join(intervention_start, by = "Site")

release_merged <- release_merged %>%
  mutate(
    Months_since_start = as.numeric(difftime(YearWeek, StartDate, units = "days")) / 30.44, # Average days in a month
    Time_interval = case_when(
      Months_since_start >= 0.97 & Months_since_start < 6  ~ "1-6 months",
      Months_since_start >= 6 & Months_since_start < 12 ~ "7-12 months",
      Months_since_start >= 12 & Months_since_start < 18 ~ "13-18 months",
      Months_since_start >= 18 & Months_since_start < 24 ~ "19-24 months",
      TRUE ~ NA_character_
    )
  )

average_coverage <- release_merged %>%
  filter(!is.na(Time_interval)) %>%
  group_by(Time_interval) %>%
  summarize(Average_Coverage = sprintf("%.4f", mean(Coverage, na.rm = TRUE)))
            

# Brazil
brazil_coverage <- read_xlsx("/Users/joyichow/Desktop/wolb_interventions/Niteroi_ZoneData_May2021_Shared.xlsx")
brazil_coverage <- brazil_coverage[c(1,2,8)]

brazil_coverage <- brazil_coverage %>%
  mutate(
    YearMonthWithDay = paste(yearmonth, "01", sep = "_"),
    Date = as.Date(YearMonthWithDay, format = "%Y_%m_%d"),
    Zone_Start_Date = case_when(
      zone == "Zone 1" ~ as.Date("2017-04-01"),
      zone == "Zone 2" ~ as.Date("2017-08-01"),
      zone == "Zone 3" ~ as.Date("2017-11-01"),
      zone == "Zone 4" ~ as.Date("2019-09-01"),
      TRUE ~ NA_Date_
    )
  ) %>%
  select(-YearMonthWithDay)


brazil_coverage <- brazil_coverage %>%
  mutate(
    Days_since_start = as.numeric(difftime(Date, Zone_Start_Date, units = "days")),
    Months_since_start = Days_since_start / 30.44, # Average days in a month
    Time_interval = case_when(
      Months_since_start >= 0.97 & Months_since_start < 6  ~ "1-6 months",
      Months_since_start >= 6 & Months_since_start < 12 ~ "7-12 months",
      Months_since_start >= 12 & Months_since_start < 18 ~ "13-18 months",
      Months_since_start >= 18 & Months_since_start < 24 ~ "19-24 months",
      TRUE ~ NA_character_
    )
  )

average_coverage <- brazil_coverage %>%
  filter(!is.na(Time_interval)) %>%
  group_by(Time_interval) %>%
  summarize(Average_Coverage = sprintf("%.2f", mean(wmel_pc, na.rm = TRUE))) # Format to two decimal places

print(average_coverage)
