# =============================================================================
# Immune Deconvolution & Gene-Immune Correlation Explorer
# =============================================================================
# Loads gene expression data (real TCGA RNA-seq via TCGAbiolinks, a synthetic
# demo dataset, or an uploaded matrix), runs immune-cell deconvolution
# (quanTIseq, TIMER, CIBERSORT, MCP-counter, xCell, EPIC, ABIS, ConsensusTME,
# ESTIMATE), compares cancer types/methods, and correlates gene expression
# with immune cell fractions.
#
# SETUP
#   install.packages(c("shiny","DT","dplyr","tidyr","tibble","stringr",
#                       "purrr","ggplot2","openxlsx","shinyjs"))
#   remotes::install_github("icbi-lab/immunedeconv")
#   BiocManager::install(c("TCGAbiolinks","SummarizedExperiment"))
#
#   CIBERSORT: place CIBERSORT.R and LM22.txt in cibersort_files/ next to
#   app.R (obtained via registration at cibersortx.stanford.edu).
#
#   Run: shiny::runApp("app.R")
#
# CITATIONS
#   CIBERSORT: Newman AM, Liu CL, Green MR, et al. Robust enumeration of
#     cell subsets from tissue expression profiles. Nat Methods.
#     2015;12(5):453-457.
#   TCGA breast cancer data: Cancer Genome Atlas Network. Comprehensive
#     molecular portraits of human breast tumours. Nature.
#     2012;490(7418):61-70. Data from the TCGA Research Network:
#     https://www.cancer.gov/tcga
#   immunedeconv: Sturm G, et al. Comprehensive evaluation of transcriptome-
#     based cell-type quantification methods for immuno-oncology.
#     Bioinformatics. 2019;35(14):i436-i445.
# =============================================================================

## 0. Packages

required_pkgs <- c("shiny", "DT", "dplyr", "tidyr", "tibble", "stringr",
                   "purrr", "ggplot2", "openxlsx", "shinyjs")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall with: install.packages(c(", paste(sprintf("'%s'", missing_pkgs), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(openxlsx)
  library(shinyjs)
})
options(shiny.maxRequestSize = 500 * 1024^2)  # allow uploads up to 500 MB
HAS_IMMUNEDECONV <- requireNamespace("immunedeconv", quietly = TRUE)
HAS_TCGABIOLINKS <- requireNamespace("TCGAbiolinks", quietly = TRUE) &&
  requireNamespace("SummarizedExperiment", quietly = TRUE)

if (HAS_IMMUNEDECONV) suppressPackageStartupMessages(library(immunedeconv))


## 0b. CIBERSORT setup
# Case-insensitive search inside cibersort_files/ so naming like
# "Cibersort.R" or "lm22.TXT" still gets picked up.

CIBERSORT_DIR <- "cibersort_files"

find_cibersort_file <- function(pattern) {
  if (!dir.exists(CIBERSORT_DIR)) return(NA_character_)
  hits <- list.files(CIBERSORT_DIR, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (length(hits) == 0) NA_character_ else hits[1]
}

CIBERSORT_R_PATH   <- find_cibersort_file("^cibersort\\.r$")
CIBERSORT_MAT_PATH <- find_cibersort_file("^lm22\\.txt$")
HAS_CIBERSORT_FILES <- !is.na(CIBERSORT_R_PATH) && !is.na(CIBERSORT_MAT_PATH)

message("Working directory: ", getwd())
message("Looking for CIBERSORT files in: ", normalizePath(CIBERSORT_DIR, mustWork = FALSE))
if (dir.exists(CIBERSORT_DIR)) {
  message("Files found in cibersort_files/: ", paste(list.files(CIBERSORT_DIR), collapse = ", "))
} else {
  message("cibersort_files/ directory does not exist at that path.")
}
message("CIBERSORT ready: ", HAS_CIBERSORT_FILES)

if (HAS_IMMUNEDECONV && HAS_CIBERSORT_FILES) {
  immunedeconv::set_cibersort_binary(CIBERSORT_R_PATH)
  immunedeconv::set_cibersort_mat(CIBERSORT_MAT_PATH)
}

## 1. Configuration

CACHE_DIR <- "tcga_cache"
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, showWarnings = FALSE)

# Maximum number of deconvolution methods a user can run in a single pass.
# Running many methods at once (especially CIBERSORT/xCell/ABIS) is slow and
# memory-heavy, so the UI enforces this cap directly (see checkboxGroupInput
# "methods" + the observeEvent that trims selections below).
MAX_METHODS <- 3

TCGA_PROJECTS <- c(
  "BRCA - Breast invasive carcinoma"               = "BRCA",
  "LUAD - Lung adenocarcinoma"                     = "LUAD",
  "LUSC - Lung squamous cell carcinoma"            = "LUSC",
  "COAD - Colon adenocarcinoma"                    = "COAD",
  "SKCM - Skin cutaneous melanoma"                 = "SKCM",
  "STAD - Stomach adenocarcinoma"                  = "STAD",
  "KIRC - Kidney renal clear cell carcinoma"       = "KIRC",
  "PRAD - Prostate adenocarcinoma"                 = "PRAD",
  "HNSC - Head & neck squamous cell carcinoma"     = "HNSC",
  "OV - Ovarian serous cystadenocarcinoma"         = "OV"
)

# Methods that require a per-sample "indications" (cancer type) argument.
# TIMER and ConsensusTME are both calibrated per tumor type.
INDICATION_METHODS <- c("timer", "consensus_tme")

# TIMER / ConsensusTME only recognize a fixed set of TCGA-style lowercase
# codes. Anything outside this list will fail for those specific methods
# (other methods are unaffected).
TIMER_SUPPORTED_INDICATIONS <- c(
  "kich","blca","brca","cesc","gbm","hnsc","kirp","lgg","lihc","luad",
  "lusc","prad","sarc","pcpg","paad","tgct","ucec","ov","skcm","stad",
  "coad","read","esca","thca","ucs","kirc","dlbc","meso","uvm","acc",
  "thym","laml","chol"
)

ALL_DECONV_METHODS <- c(
  "quanTIseq"        = "quantiseq",
  "TIMER"            = "timer",
  "CIBERSORT"        = "cibersort",
  "CIBERSORT (abs.)" = "cibersort_abs",
  "MCP-counter"      = "mcp_counter",
  "xCell"            = "xcell",
  "EPIC"             = "epic",
  "ABIS"             = "abis",
  "ConsensusTME"     = "consensus_tme",
  "ESTIMATE"         = "estimate"
)

DECONV_METHODS <- if (HAS_IMMUNEDECONV) {
  ALL_DECONV_METHODS
} else {
  c("Simple NNLS (demo)" = "demo_nnls")
}

METHOD_NOTES <- tibble::tribble(
  ~method,          ~approach,                              ~output,
  "quantiseq",      "Constrained least squares regression",  "Cell fractions (sum to 1, incl. 'other')",
  "timer",          "Batch-corrected linear regression",     "Relative abundance score per cancer type (requires per-sample cancer-type label)",
  "cibersort",      "Support Vector Regression (nu-SVR)",    "Cell fractions (sum to 1)",
  "cibersort_abs",  "SVR, no sum-to-one constraint",          "Absolute (arbitrary-unit) score",
  "mcp_counter",     "Marker-gene transcriptomic score",      "Arbitrary abundance score (not a fraction)",
  "xcell",          "ssGSEA-based enrichment + calibration", "Enrichment score (not a fraction) - needs broad real gene coverage",
  "epic",           "Constrained least squares (mRNA-weighted)", "Cell fractions (sum to 1, incl. 'other') - needs real gene coverage",
  "abis",           "Linear least-squares against sorted-cell reference", "Cell fractions - can fail with few samples/genes (singular fit)",
  "consensus_tme",  "ssGSEA against tissue-specific gene sets", "Enrichment score (not a fraction) - requires per-sample cancer-type label",
  "estimate",       "ssGSEA (stromal/immune signatures)",     "Stromal / Immune / ESTIMATE score",
  "demo_nnls",      "Non-negative least squares (toy signature)", "Cell fractions (sum to 1) - illustrative only"
)

## 2. Helper functions

#' A small, hand-built marker-gene signature used ONLY when `immunedeconv` is
#' not installed. This lets the app demonstrate the deconvolution math
#' (bulk = signature %*% fractions) without any external dependency.
#' It is NOT intended for real biological inference.
toy_signature <- function() {
  cell_types <- c("B_cell", "CD4_T_cell", "CD8_T_cell", "NK_cell", "Monocyte", "Neutrophil")
  genes <- c("MS4A1", "CD79A", "CD19", "CD3D", "CD3E", "CD4", "IL7R",
             "CD8A", "CD8B", "GZMK", "NCAM1", "NKG7", "GNLY", "KLRD1",
             "CD14", "CD68", "LYZ", "FCN1", "FCGR3B", "CSF3R", "S100A8",
             "S100A9", "ITGAM")
  mat <- matrix(0.2, nrow = length(genes), ncol = length(cell_types),
                dimnames = list(genes, cell_types))
  high <- list(
    B_cell     = c("MS4A1", "CD79A", "CD19"),
    CD4_T_cell = c("CD3D", "CD3E", "CD4", "IL7R"),
    CD8_T_cell = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMK"),
    NK_cell    = c("NCAM1", "NKG7", "GNLY", "KLRD1"),
    Monocyte   = c("CD14", "CD68", "LYZ", "FCN1"),
    Neutrophil = c("FCGR3B", "CSF3R", "S100A8", "S100A9", "ITGAM")
  )
  set.seed(1)
  for (ct in names(high)) {
    mat[high[[ct]], ct] <- 5 + rnorm(length(high[[ct]]), 0, 0.3)
  }
  list(matrix = mat, cell_types = cell_types, genes = genes)
}

#' A pool of ~300 real human gene symbols spanning housekeeping genes,
#' immune-related genes, and commonly-referenced markers across many
#' pathways. Used as "background" filler genes in the synthetic demo
#' expression matrix so that real deconvolution methods (EPIC, xCell, ABIS,
#' MCP-counter, etc.), which match against real gene-symbol signatures, have
#' meaningful gene overlap instead of near-zero overlap. This is still a
#' synthetic/illustrative dataset - not real biological expression data -
#' but it lets those methods actually run instead of failing outright.
real_gene_pool <- function() {
  c(
    "GAPDH","ACTB","TUBB","B2M","RPL13A","RPLP0","YWHAZ","SDHA","HPRT1","TBP",
    "PPIA","POLR2A","UBC","GUSB","HMBS","TFRC","PGK1","EEF1A1","RPL19","RPS18",
    "STAT1","STAT3","STAT5A","JAK1","JAK2","IRF1","IRF4","IRF7","NFKB1","RELA",
    "TNF","IL6","IL10","IL2","IL4","IL5","IL12A","IL13","IL15","IL17A",
    "IL18","IL21","IL23A","IFNG","IFNA1","IFNB1","TGFB1","CSF1","CSF2","CSF3",
    "CCL3","CCL4","CCL5","CCL19","CCL20","CCL21","CXCL10","CXCL12","CXCL13","CXCR3",
    "CXCR4","CXCR5","CCR2","CCR5","CCR7","PTPRC","CD2","CD3G","CD5","CD6",
    "CD7","CD9","CD24","CD27","CD28","CD33","CD34","CD36","CD38","CD40",
    "CD40LG","CD44","CD45RO","CD52","CD53","CD55","CD58","CD59","CD63","CD69",
    "CD70","CD72","CD74","CD80","CD81","CD83","CD84","CD86","CD93","CD96",
    "CD160","CD163","CD200","CD226","CD244","CD247","CD274","CTLA4","PDCD1","PDCD1LG2",
    "HAVCR2","LAG3","TIGIT","ICOS","ICOSLG","TNFRSF4","TNFRSF9","TNFRSF18","TNFRSF14","BTLA",
    "FOXP3","GATA3","TBX21","RORC","EOMES","PRDM1","BCL6","BCL2","BCL2L1","BAX",
    "MKI67","PCNA","TOP2A","CCNB1","CCND1","CDK1","CDK4","CDKN1A","CDKN2A","TP53",
    "MYC","EGFR","ERBB2","KRAS","PIK3CA","PTEN","BRAF","MTOR","AKT1","VEGFA",
    "VEGFB","FLT1","KDR","PDGFRB","FGFR1","MET","ALK","ROS1","RET","NTRK1",
    "COL1A1","COL1A2","COL3A1","FN1","ACTA2","VIM","FAP","PDGFRA","THY1","DCN",
    "LUM","POSTN","SPARC","MMP2","MMP9","TIMP1","PECAM1","VWF","CDH5","ENG",
    "KIT","MPO","ELANE","PRTN3","CTSG","DEFA1","CAMP","LTF","LCN2","OLFM4",
    "ARG1","NOS2","MRC1","MSR1","CD163L1","MARCO","SIGLEC1","ITGAX","ITGAM","ITGB2",
    "FCGR1A","FCGR2A","FCGR2B","FCGR3A","FCER1A","FCER1G","TLR2","TLR4","TLR7","TLR9",
    "NLRP3","MYD88","TICAM1","CGAS","STING1","OAS1","MX1","ISG15","IFIT1","IFIT3",
    "IRF3","DDX58","IFIH1","TRIM25","APOE","APOB","ALB","TTR","SERPINA1","HP",
    "TF","CP","FGA","FGB","FGG","PLG","F2","F3","VTN","AHSG",
    "GC","ORM1","AMBP","AGT","REN","ACE","NR3C1","ESR1","PGR","AR",
    "INS","GCG","SST","PPY","GHRL","LEP","ADIPOQ","PPARG","FABP4","UCP1",
    "MYH7","MYH6","TNNT2","TNNI3","ACTC1","NPPA","NPPB","RYR2","PLN","ATP2A2",
    "SCN5A","KCNQ1","KCNH2","GJA1","DES","LMNA","TTN","MYBPC3","CASQ2","CALM1",
    "HBB","HBA1","HBA2","GYPA","SPTA1","ANK1","EPB42","SLC4A1","ALAS2","EPOR",
    "GATA1","KLF1","TAL1","RUNX1","CEBPA","CEBPE","SPI1","GFI1","MPL","THPO",
    "PF4","ITGA2B","GP9","VWF","SELP","THBS1","F13A1","PLEK","TUBB1","GP1BA"
  )
}

#' Simulate a bulk RNA-seq-like expression matrix for one or more "cancer
#' types" purely for demo purposes when no internet / TCGAbiolinks is
#' available. Cancer types get different (made-up but directionally
#' plausible) immune-composition profiles so cross-cancer comparisons in the
#' app are not flat. CCL2 (a real monocyte chemoattractant, MCP-1) is wired
#' to correlate with the simulated Monocyte fraction so the Gene-Correlation
#' tab has an immediately interesting example to explore, matching the
#' CCL2-vs-monocyte example from the assignment.
#'
#' Background genes are drawn from a ~300-gene pool of REAL human gene
#' symbols (real_gene_pool()) rather than fake placeholder names, so that
#' real deconvolution methods (EPIC/xCell/ABIS/MCP-counter/etc.) find
#' meaningful overlap with their reference signatures instead of failing.
simulate_bulk_expression <- function(cancer_types, n_per_type = 15, seed = 123) {
  set.seed(seed)
  sig <- toy_signature()
  pool <- real_gene_pool()
  
  ct_alpha <- list(
    BRCA = c(3, 3, 1, 1, 2, 1), LUAD = c(1, 2, 3, 2, 2, 2),
    LUSC = c(1, 2, 3, 1, 3, 2), COAD = c(1, 1, 1, 1, 4, 3),
    SKCM = c(1, 2, 4, 3, 1, 1), STAD = c(1, 1, 2, 1, 3, 3),
    KIRC = c(1, 2, 2, 3, 2, 1), PRAD = c(2, 2, 1, 1, 3, 1),
    HNSC = c(1, 2, 3, 1, 3, 2), OV   = c(2, 2, 2, 2, 2, 2)
  )
  
  # Background genes: real gene symbols (deduplicated, excluding anything
  # already used as a marker gene) plus the three genes wired to specific
  # immune fractions below.
  bg_genes <- setdiff(unique(pool), sig$genes)
  bg_genes <- unique(c("CCL2", "CXCL9", "IL2RA", bg_genes))
  n_bg_genes <- length(bg_genes)
  
  expr_list <- list(); frac_list <- list(); meta_list <- list()
  
  for (ct in cancer_types) {
    alpha <- ct_alpha[[ct]]
    if (is.null(alpha)) alpha <- rep(1, length(sig$cell_types))
    
    fracs <- t(replicate(n_per_type, {
      x <- rgamma(length(alpha), shape = alpha, rate = 1)
      x / sum(x)
    }))
    colnames(fracs) <- sig$cell_types
    samp_ids <- paste0(ct, "_S", seq_len(n_per_type))
    rownames(fracs) <- samp_ids
    
    bulk_sig <- sig$matrix %*% t(fracs)                       # marker genes x samples
    noise <- matrix(rlnorm(length(bulk_sig), 0, 0.15), nrow = nrow(bulk_sig))
    bulk_sig <- bulk_sig * noise
    
    bg <- matrix(rlnorm(n_bg_genes * n_per_type, meanlog = 2, sdlog = 1),
                 nrow = n_bg_genes, ncol = n_per_type,
                 dimnames = list(bg_genes, samp_ids))
    bg["CCL2", ]  <- 40 * fracs[, "Monocyte"]    + rlnorm(n_per_type, 0, 0.25)
    bg["CXCL9", ] <- 35 * fracs[, "CD8_T_cell"]  + rlnorm(n_per_type, 0, 0.25)
    bg["IL2RA", ] <- 30 * fracs[, "CD4_T_cell"]  + rlnorm(n_per_type, 0, 0.25)
    
    mat <- rbind(bulk_sig, bg)
    expr_list[[ct]] <- mat
    frac_list[[ct]] <- fracs
    meta_list[[ct]] <- tibble(sample = samp_ids, cancer_type = ct)
  }
  
  list(
    expr      = do.call(cbind, expr_list),
    fractions = as_tibble(do.call(rbind, frac_list), rownames = "sample"),
    meta      = bind_rows(meta_list)
  )
}

#' Fetch (and cache) a TPM expression matrix for one TCGA project via
#' TCGAbiolinks. Only used when the "Real TCGA data" source is selected.
fetch_tcga_expression <- function(project_code, n_samples = 20) {
  cache_file <- file.path(CACHE_DIR, paste0(project_code, "_n", n_samples, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  if (!HAS_TCGABIOLINKS) {
    stop("TCGAbiolinks / SummarizedExperiment are not installed.\n",
         "Install with: BiocManager::install(c('TCGAbiolinks','SummarizedExperiment'))\n",
         "...or switch 'Data source' to Demo data.")
  }
  query <- TCGAbiolinks::GDCquery(
    project = paste0("TCGA-", project_code),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  TCGAbiolinks::GDCdownload(query, method = "api", files.per.chunk = 10)
  se <- TCGAbiolinks::GDCprepare(query)
  tpm <- SummarizedExperiment::assay(se, "tpm_unstrand")
  gene_symbols <- SummarizedExperiment::rowData(se)$gene_name
  rownames(tpm) <- make.unique(as.character(gene_symbols))
  tpm <- tpm[!is.na(rownames(tpm)) & rownames(tpm) != "", , drop = FALSE]
  
  set.seed(42)
  keep <- sample(seq_len(ncol(tpm)), min(n_samples, ncol(tpm)))
  tpm <- tpm[, keep, drop = FALSE]
  colnames(tpm) <- paste0(project_code, "_", colnames(tpm))
  
  saveRDS(tpm, cache_file)
  tpm
}

#' Parse a user-uploaded expression matrix. Handles both simple CSVs (gene
#' column + sample columns) and cBioPortal-style RSEM files, which are
#' tab-delimited and include TWO leading ID columns (Hugo_Symbol,
#' Entrez_Gene_Id) before the sample columns.
#'
#' Steps:
#'   1. Auto-detect delimiter (tab vs comma) from the file's first line.
#'   2. Identify the gene-symbol column by name (Hugo_Symbol / Gene / etc.),
#'      falling back to the first column if none match.
#'   3. Drop other known non-expression ID columns (e.g. Entrez_Gene_Id).
#'   4. Drop rows with blank/NA gene symbols.
#'   5. Collapse duplicate gene symbols by keeping the row with the highest
#'      mean expression (standard practice for multi-mapped symbols).
#'   6. Optionally log2(x + 1) transform (for raw-scale RSEM/counts data;
#'      several methods such as MCP-counter and xCell expect log space).
read_uploaded_matrix <- function(path, log_transform = FALSE) {
  first_line <- readLines(path, n = 1)
  sep <- if (grepl("\t", first_line)) "\t" else ","
  
  df <- utils::read.delim(path, sep = sep, check.names = FALSE,
                          stringsAsFactors = FALSE, quote = "")
  
  gene_col_candidates <- c("Hugo_Symbol", "Gene", "gene", "SYMBOL",
                           "Gene_Symbol", "GeneSymbol", "gene_symbol")
  gene_col <- intersect(gene_col_candidates, names(df))
  gene_col <- if (length(gene_col) > 0) gene_col[1] else names(df)[1]
  
  drop_col_names <- c("Entrez_Gene_Id", "Entrez_ID", "EntrezID", "entrez_id")
  drop_cols <- intersect(drop_col_names, names(df))
  
  genes <- as.character(df[[gene_col]])
  keep_cols <- setdiff(names(df), c(gene_col, drop_cols))
  mat <- as.matrix(df[, keep_cols, drop = FALSE])
  mode(mat) <- "numeric"
  
  valid <- !is.na(genes) & genes != ""
  genes <- genes[valid]
  mat <- mat[valid, , drop = FALSE]
  
  if (any(duplicated(genes))) {
    ord <- order(genes, -rowMeans(mat, na.rm = TRUE))
    genes <- genes[ord]
    mat <- mat[ord, , drop = FALSE]
    keep <- !duplicated(genes)
    genes <- genes[keep]
    mat <- mat[keep, , drop = FALSE]
  }
  rownames(mat) <- genes
  
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  
  if (log_transform) mat <- log2(mat + 1)
  
  mat
}

#' Build a lowercase "indications" vector (one entry per column of
#' expr_mat, in matching order) from sample_meta$cancer_type. Used by
#' TIMER and ConsensusTME, which are calibrated per cancer type and require
#' this argument. Returns NULL (with a warning) if any sample is missing a
#' recognized indication.
build_indications <- function(expr_mat, sample_meta, method_name) {
  ind_lookup <- setNames(tolower(sample_meta$cancer_type), sample_meta$sample)
  indications <- unname(ind_lookup[colnames(expr_mat)])
  
  if (any(is.na(indications)) || any(indications == "")) {
    warning(sprintf(
      "Method '%s' skipped: missing cancer_type for one or more samples (needed for every sample).",
      method_name
    ))
    return(NULL)
  }
  unsupported <- setdiff(unique(indications), TIMER_SUPPORTED_INDICATIONS)
  if (length(unsupported) > 0) {
    warning(sprintf(
      "Method '%s' skipped: unsupported cancer_type code(s) %s (must be one of: %s).",
      method_name, paste(unsupported, collapse = ", "),
      paste(TIMER_SUPPORTED_INDICATIONS, collapse = ", ")
    ))
    return(NULL)
  }
  indications
}

#' Run one or more deconvolution methods on an expression matrix and return
#' one tidy ("long") data frame: sample | cancer_type | method | cell_type | fraction
run_deconvolution <- function(expr_mat, methods, sample_meta) {
  results <- list()
  for (m in methods) {
    res <- tryCatch({
      if (m == "demo_nnls") {
        sig <- toy_signature()
        common_genes <- intersect(rownames(expr_mat), rownames(sig$matrix))
        if (length(common_genes) < 3) {
          stop("Fewer than 3 of the demo marker genes were found in this expression matrix.")
        }
        S <- sig$matrix[common_genes, , drop = FALSE]
        Y <- expr_mat[common_genes, , drop = FALSE]
        k <- ncol(S)
        frac <- vapply(colnames(Y), function(s) {
          y <- Y[, s]
          fit <- optim(rep(1 / k, k),
                       function(b) sum((y - S %*% b)^2),
                       method = "L-BFGS-B", lower = rep(0, k), upper = rep(Inf, k))
          b <- fit$par
          if (sum(b) == 0) b <- rep(1 / k, k)
          b / sum(b)
        }, FUN.VALUE = numeric(k))
        rownames(frac) <- colnames(S)
        as_tibble(frac, rownames = "cell_type")
        
      } else if (m %in% INDICATION_METHODS) {
        indications <- build_indications(expr_mat, sample_meta, m)
        if (is.null(indications)) {
          NULL
        } else {
          immunedeconv::deconvolute(expr_mat, method = m, indications = indications)
        }
        
      } else if (m %in% c("cibersort", "cibersort_abs") && !HAS_CIBERSORT_FILES) {
        warning(sprintf("Method '%s' skipped: CIBERSORT.R/LM22.txt not found in '%s/'.", m, CIBERSORT_DIR))
        NULL
        
      } else {
        immunedeconv::deconvolute(expr_mat, method = m)
      }
    }, error = function(e) {
      warning(sprintf("Method '%s' failed: %s", m, conditionMessage(e)))
      NULL
    })
    
    if (!is.null(res)) {
      results[[m]] <- res %>%
        tidyr::pivot_longer(-cell_type, names_to = "sample", values_to = "fraction") %>%
        mutate(method = m)
    }
  }
  if (length(results) == 0) return(tibble())
  bind_rows(results) %>% left_join(sample_meta, by = "sample")
}

#' Compute Spearman correlations + build one ggplot per (gene, cell_type,
#' method) combination between gene expression and immune cell fraction.
compute_gene_correlations <- function(expr_mat, deconv_long, genes, cell_types, methods) {
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) == 0 || nrow(deconv_long) == 0) return(list(table = tibble(), plots = list()))
  
  expr_long <- as_tibble(t(expr_mat[genes, , drop = FALSE]), rownames = "sample") %>%
    tidyr::pivot_longer(-sample, names_to = "gene", values_to = "expression")
  
  merged <- deconv_long %>%
    filter(cell_type %in% cell_types, method %in% methods) %>%
    inner_join(expr_long, by = "sample", relationship = "many-to-many")
  
  if (nrow(merged) == 0) return(list(table = tibble(), plots = list()))
  
  combos <- merged %>% distinct(gene, cell_type, method, cancer_type)
  
  # Base-R subsetting on purpose (avoids tidy-eval data-masking ambiguity
  # between the loop variables and the `merged` columns of the same name).
  rows <- vector("list", nrow(combos))
  for (i in seq_len(nrow(combos))) {
    g <- combos$gene[i]; ct <- combos$cell_type[i]; me <- combos$method[i]; can <- combos$cancer_type[i]
    d <- merged[merged$gene == g & merged$cell_type == ct & merged$method == me & merged$cancer_type == can, ]
    n <- nrow(d)
    rho <- NA_real_; p_value <- NA_real_
    if (n >= 3 && length(unique(d$fraction)) > 1 && length(unique(d$expression)) > 1) {
      ct_test <- suppressWarnings(cor.test(d$expression, d$fraction, method = "spearman"))
      rho <- unname(ct_test$estimate)
      p_value <- ct_test$p.value
    }
    rows[[i]] <- tibble(gene = g, cell_type = ct, method = me, cancer_type = can,
                        n = n, rho = rho, p_value = p_value)
  }
  stat_table <- bind_rows(rows) %>% arrange(desc(abs(rho)))
  
  plots <- list()
  for (i in seq_len(nrow(combos))) {
    g <- combos$gene[i]; ct <- combos$cell_type[i]; me <- combos$method[i]; can <- combos$cancer_type[i]
    d <- merged[merged$gene == g & merged$cell_type == ct & merged$method == me & merged$cancer_type == can, ]
    st <- stat_table[stat_table$gene == g & stat_table$cell_type == ct &
                       stat_table$method == me & stat_table$cancer_type == can, ]
    rho_txt <- if (nrow(st) == 1 && !is.na(st$rho)) sprintf("%.3f", st$rho) else "NA"
    p_txt   <- if (nrow(st) == 1 && !is.na(st$p_value)) sprintf("%.3g", st$p_value) else "NA"
    
    key <- paste(g, ct, me, can, sep = "__")
    plots[[key]] <- ggplot(d, aes(x = fraction, y = expression)) +
      geom_point(alpha = 0.75, size = 2.4, color = "#2C6E9B") +
      geom_smooth(method = "lm", se = FALSE, color = "#D1495B", linewidth = 0.8, na.rm = TRUE) +
      labs(
        title = paste0("Cell Type: ", ct, "   |   Gene Name: ", g,
                       "   |   Spearman's correlation coefficient: ", rho_txt),
        subtitle = paste0("Cancer type: ", can, "   |   Method: ", me,
                          "   |   p = ", p_txt, "   |   n = ", ifelse(nrow(st) == 1, st$n, nrow(d))),
        x = paste0(ct, " fraction / score"),
        y = paste0(g, " expression")
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9, color = "grey30"))
  }
  
  list(table = stat_table %>% arrange(desc(abs(rho))), plots = plots)
}

## 3. UI
ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML("
    .well { background-color: #f7f9fb; }
    h4.section-title { border-bottom: 2px solid #2C6E9B; padding-bottom: 4px; margin-top: 4px; }
    .status-msg { color: #555; font-style: italic; }
  "))),
  
  titlePanel("Immune Deconvolution & Gene-Immune Correlation Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("1. Data source", class = "section-title"),
      radioButtons("data_source", NULL,
                   choices = c("Demo data (synthetic, instant, no setup)" = "demo",
                               "Real TCGA data (needs TCGAbiolinks + internet)" = "tcga",
                               "Upload my own expression matrix" = "upload"),
                   selected = "demo"),
      
      conditionalPanel(
        condition = "input.data_source != 'upload'",
        selectizeInput("cancer_types", "Cancer type(s) (TCGA project codes)",
                       choices = TCGA_PROJECTS, selected = c("BRCA", "SKCM"),
                       multiple = TRUE,
                       options = list(placeholder = "Select one or more cancer types")),
        numericInput("n_samples", "Samples per cancer type", value = 15, min = 4, max = 200, step = 1)
      ),
      
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("upload_file", "Expression matrix (CSV/TSV: gene column + sample columns; cBioPortal RSEM format supported)",
                  accept = c(".csv", ".tsv", ".txt")),
        fileInput("upload_meta", "Optional metadata CSV (columns: sample, cancer_type)",
                  accept = c(".csv", ".tsv", ".txt")),
        checkboxInput("upload_log_transform", "Log2(x + 1) transform expression values", value = TRUE),
        helpText(class = "status-msg",
                 "Leave log2 transform on for raw-scale RSEM/count data. ",
                 "For TIMER / ConsensusTME, cancer_type must be a supported lowercase TCGA code ",
                 "(e.g. brca, luad, skcm) - see Method notes tab."),
        helpText(class = "status-msg",
                 "If you don't upload a metadata file, all samples are grouped under ",
                 "'Not specified' as a placeholder - this is NOT a real cancer type. ",
                 "TIMER and ConsensusTME will be skipped in that case, since both require ",
                 "a genuine per-sample TCGA cancer-type code.")
      ),
      
      if (!HAS_TCGABIOLINKS) helpText(class = "status-msg",
                                      "Note: TCGAbiolinks not detected \u2014 'Real TCGA data' will show an install hint if selected."),
      
      h4("2. Deconvolution method(s)", class = "section-title"),
      checkboxGroupInput("methods", NULL, choices = DECONV_METHODS,
                         selected = DECONV_METHODS[1]),
      helpText(class = "status-msg",
               sprintf("You can select up to %d methods per run \u2014 each method (especially CIBERSORT, xCell, and ABIS) can take a while to compute.", MAX_METHODS)),
      if (!HAS_IMMUNEDECONV) helpText(class = "status-msg",
                                      "Note: 'immunedeconv' not detected \u2014 using a small built-in NNLS demo method. ",
                                      "Install immunedeconv for quanTIseq / CIBERSORT / xCell / MCP-counter / EPIC / etc."),
      if (HAS_IMMUNEDECONV && !HAS_CIBERSORT_FILES) helpText(class = "status-msg",
                                                             "Note: CIBERSORT.R/LM22.txt not found in 'cibersort_files/'."),
      
      actionButton("run_btn", "Run deconvolution", class = "btn-primary", width = "100%"),
      br(), br(),
      downloadButton("download_deconv", "Download deconvolution results (.xlsx)", width = "100%")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel("Overview",
                 br(),
                 uiOutput("overview_ui")
        ),
        
        tabPanel("Deconvolution results",
                 br(),
                 h4("Cell-type fractions / scores by sample"),
                 DTOutput("results_table"),
                 br(),
                 h4("Average composition by cancer type & method"),
                 plotOutput("barplot_avg", height = 420)
        ),
        
        tabPanel("Compare cancer types & methods",
                 br(),
                 fluidRow(
                   column(6, selectInput("compare_method", "Method", choices = NULL)),
                   column(6, selectInput("compare_celltype", "Cell type", choices = NULL))
                 ),
                 plotOutput("boxplot_compare", height = 420),
                 br(),
                 h4("Method agreement (matching cell-type names only)"),
                 plotOutput("method_agreement", height = 380)
        ),
        
        tabPanel("Gene <-> Immune correlations",
                 br(),
                 fluidRow(
                   column(4, textAreaInput("genes_input", "Gene(s), comma-separated",
                                           value = "CCL2, CXCL9, IL2RA", rows = 2)),
                   column(4, selectizeInput("corr_celltypes", "Cell type(s)", choices = NULL, multiple = TRUE)),
                   column(4, selectizeInput("corr_methods", "Method(s)", choices = NULL, multiple = TRUE))
                 ),
                 actionButton("run_corr_btn", "Compute correlations", class = "btn-primary"),
                 downloadButton("download_corr", "Download correlation table (.xlsx)"),
                 br(), br(),
                 DTOutput("corr_table"),
                 br(),
                 h4("Scatter plots"),
                 uiOutput("corr_plots_ui")
        ),
        
        tabPanel("Method notes",
                 br(),
                 DTOutput("method_notes_table")
        )
      )
    )
  )
)

## 4. Server                                   
server <- function(input, output, session) {
  
  ## Enforce the MAX_METHODS cap on the method checkboxes. If a user checks
  ## more than MAX_METHODS boxes, revert to their first MAX_METHODS choices
  ## and let them know why.
  observeEvent(input$methods, {
    if (length(input$methods) > MAX_METHODS) {
      showNotification(
        sprintf("You can only run up to %d deconvolution methods at a time. Keeping your first %d selections.",
                MAX_METHODS, MAX_METHODS),
        type = "warning"
      )
      updateCheckboxGroupInput(session, "methods", selected = input$methods[seq_len(MAX_METHODS)])
    }
  }, ignoreNULL = FALSE)
  
  ##    Overview tab   
  output$overview_ui <- renderUI({
    tagList(
      h4("What is immune deconvolution?"),
      p("Bulk RNA-seq measures the average gene expression across every cell in a tissue sample - ",
        "tumor cells, T cells, B cells, macrophages, fibroblasts, and more, all mixed together. ",
        "Immune deconvolution is the computational 'unmixing' of that average signal back into ",
        "estimates of how much of each immune cell type is present, without ever physically sorting ",
        "the cells (e.g., with flow cytometry)."),
      p("Mathematically, most methods start from the same mixture model:"),
      tags$pre("bulk_expression (genes x samples)  =  signature (genes x cell types)  x  fractions (cell types x samples)"),
      p("Given a reference 'signature matrix' of marker genes for each pure cell type, the methods solve ",
        "the inverse problem for the unknown fractions - but they differ in exactly how:"),
      tags$ul(
        tags$li(strong("CIBERSORT:"), " support vector regression (nu-SVR) against a leukocyte signature matrix (LM22)."),
        tags$li(strong("quanTIseq:"), " constrained least-squares regression with an empirically-calibrated signature."),
        tags$li(strong("EPIC:"), " constrained least squares that also corrects for how much mRNA each cell type produces."),
        tags$li(strong("TIMER:"), " linear regression on immune gene sets, calibrated per cancer type - requires a cancer-type label for every sample."),
        tags$li(strong("MCP-counter / xCell / ConsensusTME:"), " marker-gene enrichment scores (rank-based, like ssGSEA) - relative scores, not fractions. ConsensusTME also requires a per-sample cancer-type label."),
        tags$li(strong("ESTIMATE:"), " ssGSEA-based stromal and immune scores, mainly used to estimate tumor purity.")
      ),
      p("See the comparison in the 'Method notes' tab, and the benchmarking paper: ",
        a("Avila Cobos et al., Nature Communications 2023",
          href = "https://www.nature.com/articles/s41467-024-50618-0", target = "_blank"), "."),
      hr(),
      h4("How to use this app"),
      tags$ol(
        tags$li("Pick a data source on the left (Demo data works instantly; Real TCGA data needs TCGAbiolinks + internet; or upload your own matrix)."),
        tags$li("Pick one or more cancer types and up to ", MAX_METHODS, " deconvolution methods, then click 'Run deconvolution'."),
        tags$li("Explore results in 'Deconvolution results' and 'Compare cancer types & methods'."),
        tags$li("Go to 'Gene <-> Immune correlations', type in genes (e.g. CCL2), pick cell types, and click 'Compute correlations' to get scatter plots with Spearman's rho.")
      ),
      hr(),
      h4("Citations"),
      tags$ul(
        tags$li("CIBERSORT: Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. ",
                em("Nat Methods."), " 2015;12(5):453-457."),
        tags$li("quanTIseq: Finotello F, Mayer C, Plattner C, et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. ",
                em("Genome Med."), " 2019;11(1):34."),
        tags$li("TIMER: Li T, Fan J, Wang B, et al. TIMER: A Web Server for Comprehensive Analysis of Tumor-Infiltrating Immune Cells. ",
                em("Cancer Res."), " 2017;77(21):e108-e110."),
        tags$li("MCP-counter: Becht E, Giraldo NA, Lacroix L, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. ",
                em("Genome Biol."), " 2016;17(1):218."),
        tags$li("xCell: Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. ",
                em("Genome Biol."), " 2017;18(1):220."),
        tags$li("EPIC: Racle J, de Jonge K, Baumgaertner P, Speiser DE, Gfeller D. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ",
                em("eLife."), " 2017;6:e26476."),
        tags$li("ABIS: Monaco G, Lee B, Xu W, et al. RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. ",
                em("Cell Rep."), " 2019;26(6):1627-1640.e7."),
        tags$li("ConsensusTME: Jim\u00e9nez-S\u00e1nchez A, Cast O, Miller ML. Comprehensive Benchmarking and Integration of Tumor Microenvironment Cell Estimation Methods. ",
                em("Cancer Res."), " 2019;79(24):6238-6246."),
        tags$li("ESTIMATE: Yoshihara K, Shahmoradgoli M, Mart\u00ednez E, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. ",
                em("Nat Commun."), " 2013;4:2612."),
        tags$li("immunedeconv (wrapper package): Sturm G, et al. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. ",
                em("Bioinformatics."), " 2019;35(14):i436-i445."),
        tags$li("TCGA breast cancer data: Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. ",
                em("Nature."), " 2012;490(7418):61-70. Data from the TCGA Research Network: ",
                a("cancer.gov/tcga", href = "https://www.cancer.gov/tcga", target = "_blank"), ".")
      )
    )
  })
  
  output$method_notes_table <- renderDT({
    datatable(METHOD_NOTES, rownames = FALSE, options = list(pageLength = 15))
  })
  
  ##Load expression data (reactive)   
  loaded_data <- eventReactive(input$run_btn, {
    withProgress(message = "Loading expression data...", value = 0.2, {
      if (input$data_source == "demo") {
        req(length(input$cancer_types) > 0)
        sim <- simulate_bulk_expression(input$cancer_types, n_per_type = input$n_samples)
        list(expr = sim$expr, meta = sim$meta, truth = sim$fractions)
      } else if (input$data_source == "tcga") {
        req(length(input$cancer_types) > 0)
        mats <- list(); metas <- list()
        n_ct <- length(input$cancer_types)
        for (i in seq_along(input$cancer_types)) {
          ct <- input$cancer_types[i]
          incProgress(0.6 / n_ct, detail = paste("Fetching", ct))
          m <- fetch_tcga_expression(ct, n_samples = input$n_samples)
          mats[[ct]] <- m
          metas[[ct]] <- tibble(sample = colnames(m), cancer_type = ct)
        }
        common_genes <- Reduce(intersect, lapply(mats, rownames))
        expr <- do.call(cbind, lapply(mats, function(m) m[common_genes, , drop = FALSE]))
        list(expr = expr, meta = bind_rows(metas), truth = NULL)
      } else {
        req(input$upload_file)
        expr <- read_uploaded_matrix(input$upload_file$datapath,
                                     log_transform = isTRUE(input$upload_log_transform))
        if (!is.null(input$upload_meta)) {
          meta <- utils::read.csv(input$upload_meta$datapath, stringsAsFactors = FALSE)
        } else {
          # Placeholder only - NOT a real TCGA cancer type. Shown in tables/
          # plots as "Not specified" so it can't be mistaken for one, and
          # build_indications() will correctly skip TIMER/ConsensusTME since
          # "not specified" isn't in TIMER_SUPPORTED_INDICATIONS.
          meta <- tibble(sample = colnames(expr), cancer_type = "Not specified")
        }
        list(expr = expr, meta = as_tibble(meta), truth = NULL)
      }
    })
  })
  
  ##    Run deconvolution (reactive)   
  deconv_results <- eventReactive(input$run_btn, {
    dat <- loaded_data()
    req(length(input$methods) > 0)
    validate(need(length(input$methods) <= MAX_METHODS,
                  sprintf("Please select at most %d deconvolution methods at a time.", MAX_METHODS)))
    withProgress(message = "Running deconvolution...", value = 0.5, {
      run_deconvolution(dat$expr, input$methods, dat$meta)
    })
  })
  
  observeEvent(deconv_results(), {
    res <- deconv_results()
    if (nrow(res) == 0) {
      showNotification("No method produced results - check the R console for warnings.", type = "error")
    } else {
      showNotification(sprintf("Deconvolution complete: %d samples x %d methods x %d cell types.",
                               dplyr::n_distinct(res$sample), dplyr::n_distinct(res$method),
                               dplyr::n_distinct(res$cell_type)), type = "message")
      updateSelectInput(session, "compare_method", choices = unique(res$method))
      updateSelectInput(session, "compare_celltype", choices = unique(res$cell_type))
      updateSelectizeInput(session, "corr_celltypes", choices = unique(res$cell_type),
                           selected = unique(res$cell_type)[seq_len(min(3, length(unique(res$cell_type))))])
      updateSelectizeInput(session, "corr_methods", choices = unique(res$method), selected = unique(res$method))
    }
  })
  
  ##Results table   
  output$results_table <- renderDT({
    res <- deconv_results()
    validate(need(nrow(res) > 0, "Click 'Run deconvolution' on the left to see results here."))
    wide <- res %>%
      mutate(fraction = round(fraction, 4)) %>%
      pivot_wider(names_from = cell_type, values_from = fraction)
    datatable(wide, rownames = FALSE, filter = "top",
              options = list(scrollX = TRUE, pageLength = 12))
  })
  
  ##Barplot of average composition   
  output$barplot_avg <- renderPlot({
    res <- deconv_results()
    validate(need(nrow(res) > 0, ""))
    avg <- res %>% group_by(method, cancer_type, cell_type) %>%
      summarise(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop")
    ggplot(avg, aes(x = cancer_type, y = mean_fraction, fill = cell_type)) +
      geom_col(position = "stack") +
      facet_wrap(~method, scales = "free_y") +
      labs(x = "Cancer type", y = "Mean fraction / score", fill = "Cell type",
           title = "Average immune composition by cancer type and method") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })
  
  ##Compare tab: boxplot across cancer types for one method/cell type   
  output$boxplot_compare <- renderPlot({
    res <- deconv_results()
    validate(need(nrow(res) > 0, "Run deconvolution first."))
    req(input$compare_method, input$compare_celltype)
    d <- res %>% filter(method == input$compare_method, cell_type == input$compare_celltype)
    validate(need(nrow(d) > 0, "No data for this method/cell type combination."))
    ggplot(d, aes(x = cancer_type, y = fraction, fill = cancer_type)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
      labs(title = paste0(input$compare_celltype, " across cancer types (", input$compare_method, ")"),
           x = NULL, y = "Fraction / score") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
  })
  
  ##    Compare tab: method agreement scatter (matching cell-type names)   
  output$method_agreement <- renderPlot({
    res <- deconv_results()
    validate(need(length(unique(res$method)) >= 2,
                  "Select two or more methods and re-run to compare method agreement."))
    wide <- res %>% select(sample, method, cell_type, fraction) %>%
      pivot_wider(names_from = method, values_from = fraction)
    methods_present <- setdiff(names(wide), c("sample", "cell_type"))
    validate(need(length(methods_present) >= 2, "Need at least two methods with overlapping cell-type names."))
    m1 <- methods_present[1]; m2 <- methods_present[2]
    d <- wide %>% filter(!is.na(.data[[m1]]), !is.na(.data[[m2]]))
    validate(need(nrow(d) > 0, "No overlapping cell-type names between the first two selected methods."))
    ggplot(d, aes(x = .data[[m1]], y = .data[[m2]], color = cell_type)) +
      geom_point(size = 2, alpha = 0.8) +
      labs(title = paste0("Method agreement: ", m1, " vs ", m2),
           x = m1, y = m2) +
      theme_minimal(base_size = 12)
  })
  
  ##Gene <-> Immune correlations   
  corr_results <- eventReactive(input$run_corr_btn, {
    dat <- loaded_data()
    res <- deconv_results()
    validate(need(nrow(res) > 0, "Run deconvolution first."))
    genes <- str_trim(str_split(input$genes_input, ",")[[1]])
    genes <- genes[genes != ""]
    validate(need(length(genes) > 0, "Enter at least one gene."))
    compute_gene_correlations(dat$expr, res, genes, input$corr_celltypes, input$corr_methods)
  })
  
  output$corr_table <- renderDT({
    cr <- corr_results()
    validate(need(nrow(cr$table) > 0,
                  "No valid gene/cell-type combinations found (check gene spelling and selections)."))
    datatable(cr$table %>% mutate(rho = round(rho, 3), p_value = signif(p_value, 3)),
              rownames = FALSE, options = list(pageLength = 10))
  })
  
  output$corr_plots_ui <- renderUI({
    cr <- corr_results()
    if (length(cr$plots) == 0) return(NULL)
    plot_outputs <- lapply(seq_along(cr$plots), function(i) {
      plotname <- paste0("corr_plot_", i)
      local({
        ii <- i
        output[[plotname]] <- renderPlot(cr$plots[[ii]])
      })
      plotOutput(plotname, height = 320)
    })
    do.call(tagList, lapply(plot_outputs, function(p) tagList(p, br())))
  })
  
  ##Downloads   
  output$download_deconv <- downloadHandler(
    filename = function() paste0("immune_deconvolution_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      res <- deconv_results()
      wide <- res %>% pivot_wider(names_from = cell_type, values_from = fraction)
      wb <- createWorkbook()
      addWorksheet(wb, "long_format"); writeData(wb, "long_format", res)
      addWorksheet(wb, "wide_format"); writeData(wb, "wide_format", wide)
      addWorksheet(wb, "method_notes"); writeData(wb, "method_notes", METHOD_NOTES)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$download_corr <- downloadHandler(
    filename = function() paste0("gene_immune_correlations_", Sys.Date(), ".xlsx"),
    content = function(file) {
      cr <- corr_results()
      wb <- createWorkbook()
      addWorksheet(wb, "correlations"); writeData(wb, "correlations", cr$table)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

## 5. Launch                                   -

shinyApp(ui = ui, server = server)