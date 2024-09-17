library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)
library(exactextractr)
library(broom)
library(stringr)
library(stringdist)

standardize_name <- function(names) {
  sapply(names, function(name) {
    if (is.na(name)) return(NA)
    name %>%
      str_to_lower() %>%
      str_replace_all("[ёЁ]", "е") %>%
      str_replace_all("[^[:alnum:]]", " ") %>%
      str_replace_all("республика|область|край|автономный округ|ао", "") %>%
      str_trim() %>%
      str_replace_all("\\s+", " ")
  })
}

find_best_match <- function(names, name_vector) {
  sapply(names, function(name) {
    if (is.na(name)) return(NA)
    distances <- stringdist(name, name_vector, method = "lv")
    name_vector[which.min(distances)]
  })
}

# 1. Load and process night lights data
tiff_files <- list.files(pattern = "composite_\\d{4}\\.tif$", full.names = TRUE, recursive=TRUE)

provinces <- gadm("RUS", level=1, path="tmp")
provinces <- st_transform(st_as_sf(provinces), crs = "EPSG:4326")

results_list <- list()

for (tiff_file in tiff_files) {
  year <- gsub(".*composite_(\\d{4})\\.tif$", "\\1", tiff_file)
  cat("Processing year:", year, "\n")
  
  raster_data <- try(rast(tiff_file), silent = TRUE)
  if (inherits(raster_data, "try-error")) {
    warning("Failed to load raster file:", tiff_file)
    next
  }
  
  provinces_proj <- if (!compareCRS(raster_data, provinces)) {
    st_transform(provinces, crs = crs(raster_data))
  } else {
    provinces
  }
  
  stats <- try(exact_extract(raster_data, st_as_sf(provinces_proj), c("mean", "median", "sum", "stdev")), silent = TRUE)
  if (inherits(stats, "try-error") || nrow(stats) == 0) {
    warning("Failed to extract statistics for year:", year)
    next
  }
  
  results_list[[year]] <- cbind(provinces_proj, stats) %>%
    as.data.frame() %>%
    mutate(year = as.integer(year))
}

all_years_stats <- do.call(rbind, results_list)
print("Data combined successfully.")

# 2. Process population data
pop_raster <- rast("rus_ppp_2020_constrained.tif")
provinces_proj <- st_transform(provinces, crs = crs(pop_raster))
pop_stats <- exact_extract(pop_raster, st_as_sf(provinces_proj), 'sum')
provinces_proj$population <- pop_stats

# 3. Calculate night lights change
nl_change <- all_years_stats %>%
  filter(year %in% c(2021, 2023)) %>%
  dplyr::select(NL_NAME_1, year, sum) %>%
  pivot_wider(names_from = year, values_from = sum, names_prefix = "sum_") %>%
  mutate(percent_change = (sum_2023 - sum_2021) / sum_2021 * 100)

# 4. Process war losses data
nl_change$std_name <- standardize_name(nl_change$NL_NAME_1)
Russian_losses_in_Ukraine$std_name <- standardize_name(Russian_losses_in_Ukraine$region)

manual_corrections <- c(
  "санкт петербург горсовет" = "санкт петербург",
  "чечено ингушская" = "чеченская",
  "карачаево черкессия" = "карачаево черкесия",
  "северная осетия алани" = "северная осетия алания",
  "магадан магаданская" = "магаданская",
  "пермская" = "пермский",
  "еврейская" = "еврейская автономная",
  "камчатская" = "камчатский",
  "чукотский" = "чукотский автономный округ"
)

nl_change$std_name <- ifelse(nl_change$std_name %in% names(manual_corrections),
                             manual_corrections[nl_change$std_name],
                             nl_change$std_name)

# 5. Match and merge datasets
matching_df <- data.frame(
  nl_name = nl_change$std_name,
  losses_name = unname(find_best_match(nl_change$std_name, Russian_losses_in_Ukraine$std_name)),
  stringsAsFactors = FALSE
) %>% filter(!is.na(nl_name))

merged_data <- nl_change %>%
  left_join(matching_df, by = c("std_name" = "nl_name")) %>%
  left_join(Russian_losses_in_Ukraine, by = c("losses_name" = "std_name")) %>%
  left_join(st_drop_geometry(provinces_proj), by = "NL_NAME_1") %>%
  mutate(deaths_per_100k = (value / population))

# 6. Analysis and visualization
# Scatterplot
ggplot(merged_data, aes(x = percent_change, y = deaths_per_100k)) +
  geom_point() +
  labs(x = "Percentage Change in Night Time Lights (2021-2023)",
       y = "Deaths per 100,000 inhabitants",
       title = "Relationship between Change in Night Time Lights and War Losses (per capita)") +
  theme_minimal()

# Linear regression
model <- lm(percent_change ~ deaths_per_100k, data = merged_data)

# Get model summary
model_summary <- summary(model)

# Create a data frame with model statistics
model_stats <- data.frame(
  R_squared = round(model_summary$r.squared, 3),
  Adj_R_squared = round(model_summary$adj.r.squared, 3),
  p_value = format.pval(model_summary$coefficients[2, 4], digits = 3)
)

# Create the scatter plot with regression line
ggplot(merged_data, aes(x = deaths_per_100k, y = percent_change)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = F) +
  labs(
    x = "Deaths per 100,000 inhabitants",
    y = "Percentage Change in Night Time Lights (2021-2023)",
    title = "Relationship between War Losses and Change in Night Time Lights",
    subtitle = paste("R² =", model_stats$R_squared, 
                     "| Adj R² =", model_stats$Adj_R_squared,
                     "| p-value =", model_stats$p_value)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_text(
    data = merged_data,
    aes(label = NL_NAME_1),
    check_overlap = TRUE,
    vjust = -0.5,
    hjust = 0.5,
    size = 3
  )

# Print model summary
print(model_summary)

# Optional: Create a residual plot
ggplot(augment(model), aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Fitted values",
    y = "Residuals",
    title = "Residual Plot"
  ) +
  theme_minimal()
# Maps
russia_data <- provinces_proj %>%
  left_join(merged_data, by = "NL_NAME_1")

ggplot() +
  geom_sf(data = russia_data, aes(fill = percent_change), color = "white", size = 0.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "% Change in Night Lights") +
  theme_minimal() +
  labs(title = "Percentage Change in Night Time Lights (2021-2023) by Region")

ggplot() +
  geom_sf(data = russia_data, aes(fill = deaths_per_100k), color = "white", size = 0.2) +
  scale_fill_viridis_c(option = "magma", name = "Deaths per 100,000") +
  theme_minimal() +
  labs(title = "War Losses (Deaths per 100,000 inhabitants) by Region")

# 7. Summary statistics
summary_stats <- merged_data %>%
  summarise(
    total_population = sum(population, na.rm = TRUE),
    total_deaths = sum(value, na.rm = TRUE),
    avg_deaths_per_100k = mean(deaths_per_100k, na.rm = TRUE),
    median_deaths_per_100k = median(deaths_per_100k, na.rm = TRUE)
  )

print(summary_stats)

# Check for unmatched regions
unmatched_regions <- merged_data %>%
  filter(is.na(value) | is.na(population)) %>%
  dplyr::select(NL_NAME_1, std_name)

print("Unmatched regions:")
print(unmatched_regions)
