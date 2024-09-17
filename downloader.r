#Check packages, probably imported unneccessary too.
library(httr)
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(geodata)
library(future)
library(future.apply)
library(progressr)
library(readr)
library(hdf5r)

# Configuration
base_url <- "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/VNP46A4/"
download_dir <- "black_marble_data"
years <- 2019:2023
countries <- c("RUS", "ITA", "DEU", "FRA")
Sys.setenv(EDL_TOKEN = "")
# Insert the token..
edl_token <- Sys.getenv("EDL_TOKEN")
if (edl_token == "") {
  stop("EDL token not found. Please set it as an environment variable 'EDL_TOKEN'.")
}

# Create directory structure
dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)
sapply(years, function(year) {
  dir.create(file.path(download_dir, as.character(year)), showWarnings = FALSE)
})

# Load country boundaries
load_country_boundaries <- function(countries, download_dir) {
  gadm_shapes <- lapply(countries, function(country) {
    gadm(country, level = 0, path = download_dir) %>% st_as_sf()
  })
  do.call(rbind, gadm_shapes) %>% st_transform(crs = 4326)
}

# Load tile shapefile
load_tiles <- function(filepath) {
  tiles <- st_read(filepath)
  st_crs(tiles) <- 4326
  return(tiles)
}

# Get intersecting tiles
get_intersecting_tiles <- function(country_shape, tiles) {
  sf::sf_use_s2(FALSE)
  intersecting <- st_intersects(tiles, country_shape, sparse = FALSE)
  tiles$TileID[apply(intersecting, 1, any)]
}

# Improved download function with retry if error
download_file <- function(url, destfile, token, max_retries = 3) {
  for (i in 1:max_retries) {
    tryCatch({
      response <- GET(url, add_headers(Authorization = paste("Bearer", token)))
      if (status_code(response) == 200) {
        writeBin(content(response, "raw"), destfile)
        return(TRUE)
      }
    }, error = function(e) {
      warning("Error downloading file: ", e$message)
    })
    Sys.sleep(2^i)  # Exponential backoff
  }
  FALSE
}

# Download data for a specific year and tile
download_data_for_tile <- function(year, tile, edl_token, base_url, download_dir) {
  year_dir <- file.path(download_dir, as.character(year))
  csv_url <- paste0(base_url, year, "/001.csv")
  csv_destfile <- file.path(year_dir, "001.csv")
  
  if (!file.exists(csv_destfile)) {
    if (!download_file(csv_url, csv_destfile, edl_token)) {
      return(NULL)
    }
  }
  
  csv_data <- read_csv(csv_destfile, col_names = FALSE, show_col_types = FALSE)
  colnames(csv_data) <- c("name")
  
  tile_pattern <- paste0("VNP46A4.A", year, "001\\.", tile, "\\..*\\.h5")
  file_names <- csv_data$name[grepl(tile_pattern, csv_data$name)]
  
  if (length(file_names) == 0) {
    message("No files found for tile ", tile, " in year ", year)
    return(NULL)
  }
  
  lapply(file_names, function(file_name) {
    destfile <- file.path(year_dir, file_name)
    if (!file.exists(destfile)) {
      url <- paste0(base_url, year, "/001/", file_name)
      if (download_file(url, destfile, edl_token)) {
        message("Downloaded: ", file_name)
      } else {
        message("Failed to download: ", file_name)
      }
    } else {
      message("File already exists: ", file_name)
    }
  })
}

# Improved file_to_raster function
file_to_raster <- function(h5_file, variable, quality_flag_rm) {
  h5_data <- h5file(h5_file, "r+")
  on.exit(h5_data$close_all())
  
  if (h5_file %>% str_detect("VNP46A4")) {
    lat <- h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/lat"]][]
    lon <- h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/lon"]][]
    out <- h5_data[[paste0("HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/", variable)]][,]
    
    if (length(quality_flag_rm) > 0) {
      qf_name <- paste0(str_replace_all(variable, c("_Num" = "", "_Std" = "")), "_Quality")
      if (qf_name %in% h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields"]]$names) {
        qf <- h5_data[[paste0("HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/", qf_name)]][,]
        out[qf %in% quality_flag_rm] <- NA
      }
    }
    
    out <- matrix(as.numeric(out), ncol = ncol(out))
    out <- t(out)
    
    r <- terra::rast(out, crs = "EPSG:4326", extent = c(min(lon), max(lon), min(lat), max(lat)))
    r <- remove_fill_value(r, variable)
    r <- apply_scaling_factor(r, variable)
    
    return(r)
  } else {
    stop("Unsupported file type")
  }
}

# Process data for a specific year
process_data_for_year <- function(year, download_dir, variable = "NearNadir_Composite_Snow_Free", quality_flag_rm = NULL, all_countries) {
  year_dir <- file.path(download_dir, as.character(year))
  pattern <- paste0("VNP46A4.A", year, "001\\..*\\.h5$")
  h5_files <- list.files(year_dir, pattern = pattern, full.names = TRUE)
  
  if (length(h5_files) == 0) {
    message("No HDF5 files found for year ", year)
    return(NULL)
  }
  
  tiff_file <- file.path(year_dir, paste0("composite_", year, ".tif"))
  if (file.exists(tiff_file)) {
    message("Composite TIFF file already exists: ", tiff_file)
    return(NULL)
  }
  
  r_list <- lapply(h5_files, function(h5_file) {
    file_to_raster(h5_file, variable, quality_flag_rm)
  })
  
  composite_raster <- if (length(r_list) > 1) do.call(terra::mosaic, c(r_list, fun = "mean")) else r_list[[1]]
  
  countries_vect <- vect(all_countries)
  if (!terra::same.crs(composite_raster, countries_vect)) {
    countries_vect <- terra::project(countries_vect, crs(composite_raster))
  }
  composite_raster_cropped <- terra::crop(composite_raster, countries_vect)
  
  terra::writeRaster(composite_raster_cropped, tiff_file, overwrite = TRUE)
  message("Created composite TIFF file: ", tiff_file)
}

# Main execution
main <- function() {
  all_countries <- load_country_boundaries(countries, download_dir)
  tiles_shp <- load_tiles("~/Russia_ntl/BlackMarbleTiles/BlackMarbleTiles.shp")
  intersecting_tiles <- get_intersecting_tiles(all_countries, tiles_shp)
  
  plan(multisession, workers = parallel::detectCores() - 1)
  options(future.rng.onMisuse = "ignore")
  
  total_tasks <- length(years) * length(intersecting_tiles)
  
  with_progress({
    p <- progressor(steps = total_tasks)
    
    future_mapply(
      function(year, tile) {
        download_data_for_tile(year, tile, edl_token, base_url, download_dir)
        p(sprintf("Completed download for year %d, tile %s", year, tile))
      },
      year = rep(years, each = length(intersecting_tiles)),
      tile = rep(intersecting_tiles, times = length(years)),
      future.seed = TRUE
    )
  })
  
  with_progress({
    p <- progressor(steps = length(years))
    
    future_lapply(years, function(year) {
      process_data_for_year(year, download_dir, all_countries = all_countries)
      p(sprintf("Processed data for year %d", year))
    }, future.seed = TRUE)
  })
}

# Run the main function
main()

