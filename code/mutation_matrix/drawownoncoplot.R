# ============================================================
#  publication_oncoplot()
# ============================================================
library(maftools)
library(ComplexHeatmap)
library(grid)

publication_oncoplot <- function(
    maf,
    genes        = NULL,
    sample_order = NULL,
    show_column_names =FALSE,
    title        = "Oncoplot",
    output       = "oncoplot_final.pdf",
    width_mm     = 297,
    height_mm    = 210
) {

  # ----------------------------------------------------------
  # 1.  VARIANT COLOURS
  # ----------------------------------------------------------
  variant_colors <- c(
    "Missense_Mutation"       = "#4E79A7",
    "Nonsense_Mutation"       = "#E15759",
    "Frame_Shift_Del"         = "#F28E2B",
    "Frame_Shift_Ins"         = "#EDC948",
    "Splice_Site"             = "#B07AA1",
    "In_Frame_Del"            = "#76B7B2",
    "In_Frame_Ins"            = "#59A14F",
    "Translation_Start_Site"  = "#FF9DA7",
    "Nonstop_Mutation"        = "#9C755F",
    "Multi_Hit"               = "#333333"
  )

  # ----------------------------------------------------------
  # 2.  ALTER FUNCTIONS
  # ----------------------------------------------------------
  alter_fun <- list(
    background = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.95,
                gp = gpar(fill = "#F7F7F7", col = NA)),
    Missense_Mutation = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#4E79A7", col = NA)),
    Nonsense_Mutation = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#E15759", col = NA)),
    Frame_Shift_Del = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#F28E2B", col = NA)),
    Frame_Shift_Ins = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#EDC948", col = NA)),
    In_Frame_Del = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#76B7B2", col = NA)),
    In_Frame_Ins = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#59A14F", col = NA)),
    Splice_Site = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#B07AA1", col = NA)),
    Translation_Start_Site = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#FF9DA7", col = NA)),
    Nonstop_Mutation = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#9C755F", col = NA)),
    Multi_Hit = function(x, y, w, h)
      grid.rect(x, y, w * 0.95, h * 0.65,
                gp = gpar(fill = "#BAB0AC", col = NA))
  )

  # ----------------------------------------------------------
  # 3.  CLINICAL ANNOTATION COLOURS  (use display names throughout)
  # ----------------------------------------------------------
  annot_cols <- list(
    MPR = c(
      "MPR"     = "#99cc66",
      "non-MPR" = "#FF6699"
    ),
    `IFN group` = c(
      "LL"   = "#35398e",
      "LH"  = "#dea562",
      "HL"  = "#79a8da",
      "HH" = "#8b1812"
    ),
    Trial = c(
      "OpACINneo" = "#FDAE61",
      "PRADO"     = "#D7191C",
      "DONIMI"    = "#7B3294"
    ),
    `EFS status` = c(          
      "0" = "#D9D9D9",
      "1" = "#E15759"
    ),
    Gender = c(
      "Male"   = "#4E79A7",
      "Female" = "#F28E2B"
    ),
    `Prim | Ln met` = c(          
      "primary" = "#333333",
      "ln_met"  = "#999999"
    ),
    Type = c(
      "SSM"      = "#66C2A5",
      "nodulair" = "#FC8D62",
      "nevoid"   = "#8DA0CB",
      "Acral"    = "#E78AC3"
    ),
    `T stage` = c(
      "T1a" = "#EDE7F6",
      "T1b" = "#B39DDB",
      "T2a" = "#7E57C2",
      "T2b" = "#5C35A8",
      "T3a" = "#4A235A",
      "T3b" = "#6D1B4E",
      "T4a" = "#8B0000",
      "T4b" = "#4A0010"
    ),
    `Primary site` = c(
      "head"       = "#4E79A7",
      "lower limb" = "#F28E2B",
      "neck"       = "#E15759",
      "torso"      = "#76B7B2",
      "upper limb" = "#59A14F"
    )
  )

  # ----------------------------------------------------------
  # 4.  GENE LIST
  # ----------------------------------------------------------
  if (is.null(genes)) {
    genes <- getGeneSummary(maf)$Hugo_Symbol[1:20]
    message("No genes supplied — using top 20 by mutation frequency.")
  }

  # ----------------------------------------------------------
  # 5.  BUILD MUTATION MATRIX
  # ----------------------------------------------------------
  mat_raw    <- maftools:::createOncoMatrix(maf, g = genes)$oncoMatrix
  all_samples <- as.character(maf@clinical.data$Tumor_Sample_Barcode)
  mat <- matrix("", nrow = length(genes), ncol = length(all_samples),
                dimnames = list(genes, all_samples))
  common_genes   <- intersect(rownames(mat_raw), rownames(mat))
  common_samples <- intersect(colnames(mat_raw), colnames(mat))
  mat[common_genes, common_samples] <- mat_raw[common_genes, common_samples]

  # ----------------------------------------------------------
  # 6.  CLINICAL DATA — filter, sort, then rename
  # ----------------------------------------------------------
  clin <- as.data.frame(maf@clinical.data)
  rownames(clin) <- clin$Tumor_Sample_Barcode

  common_samples <- intersect(colnames(mat), rownames(clin))
  clin <- clin[common_samples, , drop = FALSE]

  clin$source <- trimws(as.character(clin$source))

  # Sort by IFN group, then case_id, then source (primary before LN)
source_priority <- ifelse(clin$source == "primary", 1,
                          ifelse(clin$source == "ln_met", 2, 3))
ifn_priority <- as.integer(factor(clin$IFN, levels = c(
  "LL", "LH", "HL", "HH"
)))
target_order <- order(ifn_priority, clin$case_id, source_priority)
  clin <- clin[target_order, , drop = FALSE]

  # Sync matrix columns to clinical order
  mat <- mat[, rownames(clin), drop = FALSE]

  # ----------------------------------------------------------
  # 7.  RENAME columns BEFORE feature detection
  # ----------------------------------------------------------
  clin <- clin |>
    dplyr::rename(
      `Prim | Ln met`         = source,
      `IFN group`    = IFN,
      Type           = type,
      `Primary site` = primary_site,
      `T stage`      = Tstage,
      `EFS status`   = EFSstat
    )

  top_features    <- c("MPR", "EFS status", "IFN group")
  bottom_features <- c("Prim | Ln met", "Trial", "Gender", "Type", "T stage", "Primary site")

  available_top    <- intersect(top_features,    colnames(clin))
  available_bottom <- intersect(bottom_features, colnames(clin))

  annot_cols_top    <- annot_cols[intersect(names(annot_cols), available_top)]
  annot_cols_bottom <- annot_cols[intersect(names(annot_cols), available_bottom)]

  message("Matrix and clinical data synchronized.")

  # ----------------------------------------------------------
  # 8.  TMB
  # ----------------------------------------------------------
  sample_summary  <- getSampleSummary(maf)
  tmb_raw         <- setNames(sample_summary$total, sample_summary$Tumor_Sample_Barcode)
  tmb             <- tmb_raw[colnames(mat)]
  tmb[is.na(tmb)] <- 0

  # ----------------------------------------------------------
  # 9.  TOP ANNOTATION
  # ----------------------------------------------------------
  top_args <- list(
    TMB = anno_barplot(
      tmb,
      gp         = gpar(fill = "#4E79A7", col = NA),
      bar_width  = 0.7,
      axis_param = list(side = "left", gp = gpar(fontsize = 8)),
      height     = unit(20, "mm")
    )
  )

  if (length(available_top) > 0) {
    top_args <- c(
      top_args,
      lapply(available_top, function(f) clin[[f]]) |> setNames(available_top)
    )
  }

  top_args <- c(top_args, list(
    col                  = annot_cols_top,
    which                = "column",
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 11, fontface = "bold"),
    annotation_name_rot  = 0,

    simple_anno_size     = unit(5, "mm"),
    gap                  = unit(1, "mm"),
    border               = FALSE
  ))

  top_annot <- do.call(HeatmapAnnotation, top_args)

  # ----------------------------------------------------------
  # 10. BOTTOM ANNOTATION
  # ----------------------------------------------------------
  bottom_annot <- NULL
  if (length(available_bottom) > 0) {
    bottom_args <- c(
      lapply(available_bottom, function(f) clin[[f]]) |> setNames(available_bottom),
      list(
        col                  = annot_cols_bottom,
        which                = "column",
        annotation_name_side = "left",
        annotation_name_gp   = gpar(fontsize = 11, fontface = "bold"),
        simple_anno_size     = unit(5, "mm"),
        gap                  = unit(1, "mm"),
        border               = FALSE
      )
    )
    bottom_annot <- do.call(HeatmapAnnotation, bottom_args)
  }

  # ----------------------------------------------------------
  # 11. COLUMN SPLIT  (by IFN group)
  # ----------------------------------------------------------
  if ("case_id" %in% colnames(clin)) {
    column_split <- factor(clin$case_id, levels = unique(as.character(clin$case_id)))
  } else {
    column_split <- NULL
  }

  custom_labels <- rep("", ncol(mat))
  names(custom_labels) <- colnames(mat)
  
  # Loop through each patient pair/group
  for (cid in unique(clin$case_id)) {
    # Find all columns belonging to this patient
    idx <- which(clin$case_id == cid)
    
    # Pick the "middle" index (if there are 2 samples, this grabs the 1st)
    mid_idx <- idx[ceiling(length(idx) / 2)]
    
    # Assign the annonID to just that one column
    custom_labels[mid_idx] <- as.character(clin$annonID[mid_idx])
  }

  # ----------------------------------------------------------
  # 12. BUILD ONCOPRINT
  # ----------------------------------------------------------
  ht <- oncoPrint(
    mat,
    alter_fun               = alter_fun,
    alter_fun_is_vectorized = FALSE,
    col                     = variant_colors,
    top_annotation          = top_annot,
    bottom_annotation       = bottom_annot,
    right_annotation        = NULL,
    show_column_names       = TRUE,                 
    column_labels           = custom_labels,
    column_names_side       = "bottom",
    column_names_rot        = 45,                  
    column_order            = colnames(mat),
    column_split            = column_split,
    column_gap              = unit(0.8, "mm"),
    row_gap                 = unit(0.4, "mm"),
    row_names_gp            = gpar(fontsize = 11),
    column_names_gp         = gpar(fontsize =  11),
    pct_gp                  = gpar(fontsize = 10),
    column_title_gp         = gpar(fontsize = 13, fontface = "bold"),
    column_title            = title,
    show_pct                = TRUE,
    remove_empty_columns    = FALSE
  )

  # ----------------------------------------------------------
  # 13. DRAW HELPER
  # ----------------------------------------------------------
  .draw <- function() {
    draw(
      ht,
      heatmap_legend_side    = "bottom",
      annotation_legend_side = "bottom",
      legend_grouping        = "original",
      padding                = unit(c(3, 3, 3, 3), "mm")
    )
  }

  # ----------------------------------------------------------
  # 14. EXPORT OR DISPLAY
  # ----------------------------------------------------------
  if (!is.null(output)) {
    ext <- tolower(tools::file_ext(output))
    if (ext == "pdf") {
      pdf(output, width = width_mm / 25.4, height = height_mm / 18)
      .draw()
      dev.off()
      message("Saved: ", output)
    } else if (ext %in% c("png", "tiff", "jpg", "jpeg")) {
      grDevices::png(output, width = width_mm, height = height_mm,
                     units = "mm", res = 300)
      .draw()
      dev.off()
      message("Saved: ", output)
    } else {
      stop("Unsupported file extension: ", ext, ". Use .pdf, .png, .tiff, or .jpg")
    }
  }

  .draw()
  invisible(ht)
}
