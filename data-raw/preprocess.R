# preprocess data -- read from data-raw, preprocess, output as .RData files to data

library(lubridate)

package_loc <- "C:/Reich/dengue-ssr-prediction"

for(loc in c("San_Juan", "Iquitos")) {
	# read in
	dat <- read.csv(file.path(package_loc, "data-raw", paste0(loc, "_Training_data.csv")))
	
	# convert dates
	dat$week_start_date <- ymd(dat$week_start_date)
	
	# save (with name same as file name)
	ds_name <- paste0(loc, "_train")
	assign(ds_name, dat)
	save(list = ds_name, file = file.path(package_loc, "data", paste0(ds_name, ".RData")))
}
