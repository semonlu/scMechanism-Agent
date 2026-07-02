# Course-derived reference script
# English filename: 02_install_github_r_packages.R
# Original course path: see source_manifest.csv
# Role: dependency inventory
# Adaptations applied for the skill package:
# - Filename flattened and translated to English.
# - Interactive working-directory selection replaced with SC_WORK_DIR/getwd() when present.
# - Example-specific object names, thresholds, metadata columns, and local filenames still require review.
# - Prefer scripts/course_adapted/ for runnable project workflows.
####################################################
####################################################
#安装GitHub包####
##生成Github API####
usethis::create_github_token()
##设置Github API###
usethis::edit_r_environ()
##GitHub R包下载####
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE,force = TRUE)
devtools::install_github("Moonerss/CIBERSORT",force = TRUE)
devtools::install_github('dviraran/xCell',force = TRUE)
devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source",force = TRUE)
devtools::install_github("cansysbio/ConsensusTME",force = TRUE)
devtools::install_github("cit-bioinfo/mMCP-counter",force = TRUE)
devtools::install_github("omnideconv/immunedeconv",force = TRUE)
devtools::install_github("Hy4m/linkET", force = TRUE)
devtools::install_github("chuiqin/irGSEA", force =T)
devtools::install_github("jinworks/CellChat", force =T)
devtools::install_github("jokergoo/circlize", force =T)
devtools::install_github("navinlabcode/copykat", force =T)
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', force =T)
devtools::install_github('immunogenomics/presto', force =T)
devtools::install_github("campbio/decontX",force=T)
devtools::install_github(repo="powellgenomicslab/scPred", ref="9f407b7436f40d44224a5976a94cc6815c6e837f", force = TRUE)
devtools::install_github(repo = 'genecell/COSGR',force=T)
devtools::install_github("JerryZhang-1222/starTracer",force=T)
devtools::install_github('cole-trapnell-lab/monocle3',force=T)
devtools::install_github('smorabit/hdWGCNA', ref='dev',force=T)
devtools::install_github('xuranw/MuSiC',force=T)
devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org",force=T)
devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org",force=T)
devtools::install_github("jackbibby1/SCPA",force=T)
####################################################
####################################################


