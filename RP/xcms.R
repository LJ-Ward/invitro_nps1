library(stringr)
library(tibble)
library(dplyr)
library(xcms)
library(Spectra)
library(MsExperiment)
library(CAMERA)
library(openxlsx)
library(SummarizedExperiment)

##############################################################################
#
# Don't touch anything below. 
#
# All settings are defined in settings.yaml
#
# This script will output the same data but in two different formats.
# 1. Rds file with XcmsExperiment object (also exported as excel file)
# 2. Rds file with SummarisedExperiment object (also exported as excel file)
# Each excel files will have 3 sheets. 
#  - One sheet with the expression matrix (feature intensities) and 
#  - One sheet with data about the rows (fatures)
#  - One sheet with data about the columns (samples)
#
##############################################################################

x <- commandArgs(trailingOnly = TRUE)
if (length(x) > 0) {
  settings_file <- x[1]
}  else {
  settings_file <- "settings.yaml" 
}

logfile <- paste("LOG_", settings_file, ".log", sep="")

log <- function(...) {
  cat(as.character(round(Sys.time())), ..., "\n", file=logfile, append=T)
  cat(as.character(round(Sys.time())), ..., "\n")
}

log("#")
log("#")
log("Welcome! Starting new run!")
#
# Read settings
#
if (!file.exists(settings_file)) { 
  log("Settings file does not exist (", settings_file, ")")
  stop("Settings file does not exist (", settings_file, ")")
}
settings <- yaml::read_yaml(file = settings_file)
log("Successfully read settings from ", settings_file)

print(settings)

#
# Read metadata
#
file_metadata <- settings$general$sample_metadata 
if (!file.exists(file_metadata)) { 
  stop("Metadata file does not exist") 
}

# Check various formats
if(str_detect(file_metadata, pattern = "\\.xlsx$")) {
  pd <- openxlsx::read.xlsx(file_metadata) %>% as_tibble()
} else {
  pd <- read_delim(file_metadata, progress=F, show_col_types=F)
}
log("Successfully read metadata from", file_metadata)

#
# Test for filename of mzml
if (!"filename_mzml" %in% names(pd)) {
  log("Didn't find column 'filename_mzlm' in metadata, so I will try and find the files...")
  pd <- pd %>% mutate(filename_mzml = NA)
  if(is.null(settings$general$input_path)) {
    log("Error: missing input_path - if you mean current directory, please put '.' as input path")
    stop("Oups! Missing input path in settings!!") 
  }
  for (i in 1:nrow(pd)) {
    pattern <- pd$mzml_id[i]
    x <- list.files(pattern = paste(pattern, "*.mzml",sep=""), 
               path = settings$general$input_path, 
               full.names = T, recursive = T, ignore.case = T)
    if(length(x)==0) {
      log("Could not find input file for: ", i,pattern)
      stop("Could not find input file for: ", i,pattern) 
    }
    if(length(x)>1) {
      message <- paste("Problem with too many files found for:", i,pattern, paste(x))
      log(message)
      stop(message) 
    }
    # No problems
    log(paste(i, pattern, "-->", x))
    pd$filename_mzml[i] <- x
    
  }  
}

#
# Test for sample_group - if not present we make it a single groupChromPeaks
if (!"sample_group" %in% names(pd)) {
  log("Didn't find any sample groups defined, setting all to the same")
  pd <- pd %>% mutate(sample_group = "sample_group_1")
}

#
# Test input files present
for (infile in pd$filename_mzml) {
  if (!file.exists(infile)) { 
  }
}
log("All good, all", nrow(pd), "infiles are present.")

#
# Creating output directories
dir.create(settings$general$output_path, showWarnings = F, recursive = T)
log("Output folder", settings$general$output_path, "created")

#
# Defining output filenames
pieces <- str_split_1(basename(settings$general$sample_metadata),fixed("."))
output_file_base <- paste(settings$general$output_path, "/xcms_output_",pieces[1:length(pieces)-1], "_",  settings$general$run_id,sep="")

#
# Parallelization
ncores     <- settings$general$ncores     # Number of CPUs you want to use in parallel
chunk_size <- settings$general$chunk_size # Number of files to peakfind simultaneously

if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(ncores)))
} else {
  register(bpstart(SnowParam(ncores)))
}

#
# Read all mzml files etc.
#
log("Now reading", nrow(pd), "mzml files for xcms processing")
rawData <- readMsExperiment(spectraFiles = pd$filename_mzml, sampleData = pd)
log("Hooray!! I just read", length(rawData),"files!")
print(rawData)

#rawData <- readMsExperiment(spectraFiles = pd$filename_mzml, sampleData = pd, backend=MsBackendMzR())
# Development
#x <- settings$xcms_parameters$prefiltering
#f1 <- Spectra::Spectra(pd$filename_mzml[1],backend= MsBackendMzR())
#f2 <- Spectra::Spectra(pd$filename_mzml[2],backend= MsBackendMzR())
# rawData <- MsExperiment(experimentFiles = MsExperimentFiles(mzmls =pd$filename_mzml), sampleData = pd, spectra = c(f1,f2))

#
# Filter for specified rt and mz ranges
x <- settings$xcms_parameters$prefiltering
log("Filtering on retention time to only keep rt values between", x$rtmin, "and", x$rtmax)
rawData <- filterRt(rawData, rt = c(x$rtmin,x$rtmax))
log("Done")
log("Filtering on mz range to only keep mz values between", x$mzmin,"and", x$mzmax)
rawData <- filterMzRange(rawData, mz = c(x$mzmin,x$mzmax))
log("Done")

###############################
#
# XCMS workflow
#
###############################

#
# Finding peaks
log("Running findChromPeaks")
x <- settings$xcms_parameters$centwave
cwp <- CentWaveParam(peakwidth = eval(parse(text = x$peakwidth)),
                     snthresh  = x$snthresh,
                     ppm       = x$ppm,
                     prefilter = eval(parse(text = x$prefilter)),
                     mzdiff    = x$mzdiff)
#
# Chunksize set to -1 for parallel reading of mzML 
# https://github.com/sneumann/xcms/issues/781
processedData <- findChromPeaks(rawData, param = cwp, chunkSize=-1) 
log("Done")
#
# Refining
log("Running refineChromPeaks")
mpp <- MergeNeighboringPeaksParam(expandRt = 4)
processedData <- refineChromPeaks(processedData, mpp, chunkSize=chunk_size)
log("Done,found", nrow(processedData@chromPeaks), "peaks in MS1")
#
# Grouping peaks
log("Running groupChromPeaks")
x <- settings$xcms_parameters$peak_grouping1
pdp <- xcms::PeakDensityParam(sampleGroups = sampleData(processedData)$sample_group,
                              maxFeatures  = x$maxFeatures,
                              bw           = x$bw,
                              minFraction  = x$minFraction,
                              binSize      = x$binSize)
processedData <- groupChromPeaks(processedData, param = pdp) 
log("Done. After refinement we have", nrow(processedData@chromPeaks), "peaks.")
#
# Rt alignment
log("Running adjustRtime")
x <- settings$xcms_parameters$alignment
pdp <- xcms::PeakGroupsParam(smooth      = x$smooth,
                             minFraction = x$minFraction,
                             span        = x$span,
                             extraPeaks  = x$extraPeaks,
                             family      = x$family)
processedData <- adjustRtime(processedData, param= pdp, chunkSize=chunk_size)
log("Done adjusting", nrow(processedData@chromPeaks), "peaks.")
#
# Group peaks
log("Running groupChromPeaks")
x <- settings$xcms_parameters$peak_grouping2
pdp <- xcms::PeakDensityParam(sampleGroups = sampleData(processedData)$sample_group,
                              maxFeatures  = x$maxFeatures,
                              bw           = x$bw,
                              minFraction  = x$minFraction,
                              binSize      = x$binSize)
processedData <- groupChromPeaks(processedData, param = pdp) 
log("Done. After grouping we have ", nrow(processedData@featureDefinitions), "features.")
#
# Fill peaks
log("Running fillChromPeaks")
processedData <- fillChromPeaks(processedData, param = ChromPeakAreaParam(), chunkSize=chunk_size)
log("Done. After filling we have", 
    nrow(processedData@chromPeaks),
    "peaks and",
    nrow(processedData@featureDefinitions), "features.")

#
# CAMERA adduct detection
log("Running CAMERA annotation")
x <- settings$CAMERA_parameters
rules     <- read.csv(x$camera_rules)
xset_integrated <- as(processedData, "xcmsSet")
sampclass(xset_integrated) <- sampleData(processedData)$sample_group
sampnames(xset_integrated) <- sampleData(processedData)$mzml_id
xsa      <- xsAnnotate(xset_integrated)
xsaF     <- groupFWHM(xsa, perfwhm=0.3,intval="into")
xsaI     <- findIsotopes(xsaF,mzabs=0.001,ppm=3,intval="into",maxcharge=8,maxiso=5,minfrac=0.2)
xsaIC    <- groupCorr(xsaI,cor_eic=0.75)
xsaFA    <- findAdducts(xsaIC,polarity=x$polarity,rules=rules, multiplier=2,ppm=3,mzabs=0.001)
camera_results <- getPeaklist(xsaFA)[,c("isotopes","adduct","pcgroup")]
log("Done")

#
# Convert to cool object (smaller size)
# https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html#27_Final_result
log("Creating summarized experiment object")
res <- quantify(processedData, value = "into", method = "sum")
log("Done")

#
# Add camera output to res object
log("Adding CAMERA annotation to results")
x <- cbind(rowData(res), camera_results)
rowData(res) <- x
x <- cbind(featureDefinitions(processedData), camera_results)
featureDefinitions(processedData) <- x
log("Done")

#
# Export to excel
log("Exporting xcms experiment Excel file")

outfile <- paste0(output_file_base, "_xcmsexp.xlsx")

openxlsx::write.xlsx(
  list(
    sampleData = sampleData(processedData),
    featureDefinitions = featureDefinitions(processedData),
    featureValues = featureValues(processedData)
  ),
  file = outfile
)

log("Done")
log("All done. Thank you for your time.")

# Copy log to result dir to keep it with output
file.copy(from = logfile, to = settings$general$output_path)
