args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = "") {
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}

has_flag <- function(name) any(args == name)

course_root <- arg_value("--course-root", "")
skip_cran <- has_flag("--skip-cran")
skip_bioc <- has_flag("--skip-bioc")
skip_github <- has_flag("--skip-github")
skip_local <- has_flag("--skip-local")
force <- has_flag("--force")

options(repos = c(CRAN = "https://cran.r-project.org"))
options(timeout = max(3600, getOption("timeout")))

message("R version: ", paste(R.version$major, R.version$minor, sep = "."))
if (nzchar(course_root)) message("Course root: ", course_root)

install_if_missing <- function(pkgs, source = "CRAN") {
  pkgs <- unique(pkgs[nzchar(pkgs)])
  missing <- pkgs[force | !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) {
    message("All ", source, " packages are already installed.")
    return(invisible(TRUE))
  }
  message("Installing ", source, " packages: ", paste(missing, collapse = ", "))
  if (source == "CRAN") {
    install.packages(missing, dependencies = TRUE)
  } else if (source == "Bioconductor") {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
  invisible(TRUE)
}

cran_packages <- c(
  "remotes","devtools","RANN","RColorBrewer","cowplot","Seurat","SeuratObject",
  "fitdistrplus","future","future.apply","ggplot2","ggrepel","ggridges","httr",
  "ica","irlba","jsonlite","leiden","lmtest","matrixStats","miniUI","patchwork",
  "pbapply","plotly","png","progressr","purrr","Rcpp","RcppAnnoy","reticulate",
  "rlang","ROCR","Rtsne","scales","scattermore","sctransform","shiny",
  "spatstat.explore","spatstat.geom","tibble","uwot","RcppEigen","RcppProgress",
  "CCA","clustree","aplot","circlize","data.table","doParallel","doRNG","dplyr",
  "ggfun","gghalves","ggplotify","ggsci","magrittr","Matrix","msigdbr","pagoda2",
  "plyr","pointr","RcppML","readr","reshape2","RMTstat","R.utils",
  "RobustRankAggreg","roxygen2","stringr","tidyr","tidyselect","tidytree","VAM",
  "harmony","hdf5r","assertthat","randomcoloR","NMF","AnnoProbe","clustermole",
  "fastDummies","qs","igraph","scCustomize","paletteer","miscTools","SCINA",
  "phylogram","ape","car","dendextend","e1071","fastSave","foreach","ggalluvial",
  "ggExtra","ggforce","ggpubr","ggthemes","gridExtra","ridge","vioplot","WGCNA",
  "tidyverse"
)

bioc_packages <- c(
  "DOSE","clusterProfiler","enrichplot","org.Hs.eg.db","limma","edgeR",
  "maftools","ComplexHeatmap","preprocessCore","DESeq2","TCGAbiolinks","sva",
  "GenomicFeatures","GEOquery","affy","impute","affyPLM","annotate",
  "affycoretools","graph","gcrma","CLL","AnnotationDbi","oligo","singscore",
  "quantiseqr","GSVA","monocle","SCpubr","UCell","AUCell","BiocParallel",
  "decoupleR","fgsea","ggtree","GSEABase","Nebulosa","scde",
  "SummarizedExperiment","viper","sparseMatrixStats","celldex","tricycle",
  "GeneOverlap","TOAST","infercnv","metapod","bluster","scran","celda",
  "glmGamPoi","SingleR","wk","units","s2","classInt","warp","biglm","sf",
  "spData","slider","ResidualMatrix","speedglm","spdep","rsample","pscl","grr",
  "batchelor","scDblFinder","SingleCellExperiment","scater"
)

github_specs <- list(
  list(repo = "GfellerLab/EPIC"),
  list(repo = "Moonerss/CIBERSORT"),
  list(repo = "dviraran/xCell"),
  list(repo = "ebecht/MCPcounter", subdir = "Source"),
  list(repo = "cansysbio/ConsensusTME"),
  list(repo = "cit-bioinfo/mMCP-counter"),
  list(repo = "omnideconv/immunedeconv"),
  list(repo = "Hy4m/linkET"),
  list(repo = "chuiqin/irGSEA"),
  list(repo = "jinworks/CellChat"),
  list(repo = "jokergoo/circlize"),
  list(repo = "navinlabcode/copykat"),
  list(repo = "chris-mcginnis-ucsf/DoubletFinder"),
  list(repo = "immunogenomics/presto"),
  list(repo = "campbio/decontX"),
  list(repo = "powellgenomicslab/scPred", ref = "9f407b7436f40d44224a5976a94cc6815c6e837f"),
  list(repo = "genecell/COSGR"),
  list(repo = "JerryZhang-1222/starTracer"),
  list(repo = "cole-trapnell-lab/monocle3"),
  list(repo = "smorabit/hdWGCNA", ref = "dev"),
  list(repo = "xuranw/MuSiC"),
  list(repo = "jackbibby1/SCPA")
)

version_specs <- list(
  list(package = "crossmatch", version = "1.3.1"),
  list(package = "multicross", version = "2.1.0")
)

local_archives <- c(
  "openai.tar.gz",
  "hdWGCNA-0.3.00.tar.gz"
)

if (!skip_cran) install_if_missing(cran_packages, "CRAN")
if (!skip_bioc) install_if_missing(bioc_packages, "Bioconductor")

if (!skip_github) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  for (spec in github_specs) {
    label <- spec$repo
    message("Installing GitHub package: ", label)
    tryCatch({
      args <- list(repo = spec$repo, upgrade = "never", dependencies = TRUE, force = force)
      if (!is.null(spec$ref)) args$ref <- spec$ref
      if (!is.null(spec$subdir)) args$subdir <- spec$subdir
      do.call(remotes::install_github, args)
    }, error = function(e) {
      message("FAILED GitHub package ", label, ": ", conditionMessage(e))
    })
  }

  for (spec in version_specs) {
    message("Installing version-pinned package: ", spec$package, " ", spec$version)
    tryCatch({
      remotes::install_version(spec$package, version = spec$version, upgrade = "never", force = force)
    }, error = function(e) {
      message("FAILED version package ", spec$package, ": ", conditionMessage(e))
    })
  }
}

if (!skip_local && nzchar(course_root)) {
  for (archive_name in local_archives) {
    matches <- list.files(course_root, pattern = paste0("^", archive_name, "$"), recursive = TRUE, full.names = TRUE)
    archive <- if (length(matches)) matches[[1]] else file.path(course_root, archive_name)
    if (file.exists(archive)) {
      message("Installing local archive: ", archive)
      tryCatch({
        remotes::install_local(archive, upgrade = "never", dependencies = TRUE, force = force)
      }, error = function(e) {
        message("FAILED local archive ", archive, ": ", conditionMessage(e))
      })
    } else {
      message("Missing local archive: ", archive_name)
    }
  }
}

required_imports <- c(
  "AnnoProbe","ape","assertthat","BiocParallel","car","CCA","CellChat","celldex",
  "circlize","clusterProfiler","clustree","ComplexHeatmap","copykat","COSG",
  "cowplot","data.table","decontX","dendextend","doParallel","DOSE",
  "DoubletFinder","dplyr","e1071","fastSave","foreach","future","future.apply",
  "ggalluvial","ggExtra","ggforce","ggplot2","ggpubr","ggsci","ggthemes",
  "glmGamPoi","GSEABase","GSVA","harmony","hdf5r","hdWGCNA","igraph",
  "infercnv","irGSEA","limma","Matrix","miscTools","monocle","monocle3",
  "MuSiC","NMF","openai","org.Hs.eg.db","patchwork","phylogram","plyr",
  "preprocessCore","presto","qs","R.utils","randomcoloR","remotes","reshape2",
  "ridge","scales","scater","scattermore","scCustomize","scDblFinder","SCINA",
  "scPred","SCpubr","Seurat","SeuratObject","SingleCellExperiment","SingleR",
  "starTracer","stringr","tidyverse","tricycle","UCell","vioplot","WGCNA"
)

check <- data.frame(
  package = required_imports,
  installed = vapply(required_imports, requireNamespace, logical(1), quietly = TRUE)
)
check$version <- vapply(required_imports, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else ""
}, character(1))

out_file <- file.path(getwd(), "seuratv5_r_package_check.tsv")
utils::write.table(check, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("R package check written to: ", out_file)
message("Installed required imports: ", sum(check$installed), "/", nrow(check))
if (any(!check$installed)) {
  message("Missing required imports: ", paste(check$package[!check$installed], collapse = ", "))
}
