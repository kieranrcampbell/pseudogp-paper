
## Analysis of p-values


## get the paths right
base_dir <- ""
system <- devtools::session_info()$platform$system
if(length(grep("darwin", system)) > 0) {
  # we're on the mac
  base_dir <- "~/mount/"
} else {
  # we're on linux
  base_dir <- "/net/isi-scratch/kieran/"
}

print(paste("Using base_dir", base_dir))

devtools::load_all(paste0(base_dir, "switch/sctools"))

source(paste0(base_dir, "GP/pseudogp2/stan/gbio/multi/diffexpr/common.R"))

outputfile <- paste0(base_dir, "GP/pseudogp2/stan/gbio/varygamma/diffexpr/g2.pdf")
h5file <- paste0(base_dir, "GP/pseudogp2/data/varygamma/g2.hdf5")
pstfile <- paste0(base_dir, "GP/pseudogp2/data/varygamma_traces.h5")

fdrfile <- paste0(base_dir, "GP/pseudogp2/stan/gbio/varygamma/diffexpr/g2.txt")

##------- monocle ONLY
source(paste0(base_dir, "GP/pseudogp2/stan/diffexpr/monocle/prep_data.R"))
sce <- load_data()
## end


sigList <- pvalsFromHDF5(h5file)

pst <- h5read(pstfile, "g2/pst")

## need sce & pst
generatePlots(sigList, sce, pst, outputfile, fdrfile)
