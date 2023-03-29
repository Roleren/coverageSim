#' Plot transcript coordinates window
#' @export
ORF_frame_plot_window <- function(tx_window, sim_reads, uORFs,
                                  cds_region = startRegion(cds_RPL, upstream = 0, downstream = 100),
                                  cov_real = coveragePerTiling(tx_window, merged_gr, is.sorted = T, as.data.table = T, withFrames = T, fraction = "Real"),
                                  cov_sim = coveragePerTiling(tx_window, sim_reads, is.sorted = T, as.data.table = T, withFrames = T, fraction = "Sim"),
                                  break_string = "    ", y_breaks = 50, prop = c(5, 2), x_lab = "Position") {

  cov_combined <- rbindlist(list(cov_real, cov_sim))
  cov_combined[, frame := factor(frame)]
  cov_plot <- ggplot(cov_combined, aes(y = count, x = position, fill = frame)) +
    geom_col() +
    scale_y_continuous(breaks=y_breaks) +
    facet_wrap(~ fraction, ncol = 1, strip.position = "right") +
    theme_classic() + theme(legend.position = "none",
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            plot.margin = margin(0, 0, 0, 0, "cm")); cov_plot

  display_genes <- c(uORFs, cds_region); names(display_genes) <- c(paste0("U", seq(length(uORFs))), "CDS")
  ORF_Panel <- gene_model_panel(display_range = tx_window, display_genes,
                                custom_regions = NULL, viewMode = "any")[[1]]
  ORF_Panel <- ORF_Panel +
    ylab("") +
    xlab(x_lab) +
    xlim(0, widthPerGroup(tx_window, F)) +
    scale_y_continuous(labels = function(breaks) {rlang::rep_along(breaks, break_string)}) +
    theme_classic() +
    facet_wrap(~ "", ncol = 1, strip.position = "right") +
    theme(strip.text.y = element_text(),
          strip.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  plot <- gridExtra::arrangeGrob(cov_plot, ORF_Panel,
                                 layout_matrix = t(matrix(c(rep(1, prop[1]), rep(2, prop[2])), ncol = sum(prop))),
                                 padding = unit(0.0, "line"))
  return(plot)
}

gene_model_panel <- function (display_range, annotation, frame = 1, custom_regions,
                              viewMode) {
  if (!is.null(custom_regions)) {
    same_names <- names(custom_regions) %in% names(annotation)
    names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names],
                                               "_1", sep = "")
    annotation <- c(annotation, custom_regions)
  }
  overlaps <- subsetByOverlaps(annotation, display_range,
                               type = ifelse(viewMode == "tx", "within", "any"))
  if (length(overlaps) > 0) {
    plot_width <- widthPerGroup(display_range)
    onames <- rep(names(overlaps), numExonsPerGroup(overlaps,
                                                    FALSE))
    overlaps <- unlistGrl(overlaps)
    names(overlaps) <- onames
    overlaps$rel_frame <- getRelativeFrames(overlaps)
    rel_frame <- getRelativeFrames(overlaps)
    overlaps <- subsetByOverlaps(overlaps, display_range)
    intersections <- trimOverlaps(overlaps, display_range)
    intersections <- groupGRangesBy(intersections)
    locations <- pmapToTranscriptF(intersections, display_range)
    layers <- geneTrackLayer(locations)
    locations <- unlistGrl(locations)
    rel_frame <- getRelativeFrames(overlaps)
    names(rel_frame) <- names(overlaps)
    if (length(rel_frame) != length(locations))
      rel_frame <- selectFrames(rel_frame, locations)
    locations$rel_frame <- rel_frame
    cols <- colour_bars(locations, overlaps, display_range)
    locations <- ranges(locations)
    blocks <- c(start(locations), end(locations))
    names(blocks) <- rep(names(locations), 2)
    blocks <- sort(blocks)
    lines_locations <- blocks[!(blocks %in% c(1, plot_width))]
    rect_locations <- locations
    locations <- resize(locations, width = 1, fix = "center")
    labels_locations <- start(locations)
    too_close <- labels_locations < 0.02 * plot_width
    too_far <- labels_locations > 0.98 * plot_width
    labels_locations[too_close] <- 1
    labels_locations[too_far] <- plot_width
    hjusts <- rep("center", length(labels_locations))
    hjusts[too_close] <- "left"
    hjusts[too_far] <- "right"
    gene_names <- names(locations)
    custom_region_names <- which(names(lines_locations) %in%
                                   names(custom_regions))
    names(lines_locations) <- rep("black", length(lines_locations))
    names(lines_locations)[custom_region_names] <- "orange4"
    result_dt <- data.table(layers = layers, rect_starts = start(rect_locations),
                            rect_ends = end(rect_locations), cols = cols, labels_locations = labels_locations,
                            hjusts = hjusts, gene_names = gene_names)
  }
  else {
    result_dt <- data.table()
    lines_locations <- NULL
  }
  return(list(result_dt, lines_locations))
}
