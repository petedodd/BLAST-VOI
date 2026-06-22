## this is a set of utility functions for flux calculations
## NOTE various things must be loaded in the global environment for these to work

## flux indices/names
note_flux_idx <- grep("cum_note_flux", BLASTtbmod::get_cols)
(nmz <- grep("cum_note_flux", BLASTtbmod::get_cols, value = TRUE)) # check
inf_flux_idx <- grep("cum_inf_flux", BLASTtbmod::get_cols)
(nmz <- grep("cum_inf_flux", BLASTtbmod::get_cols, value = TRUE)) # check

## --- extract fluxes from model runs
get_fluxes <- function(output) {
  ## get infection flux
  tmp <- t(apply(
    output[inf_flux_idx, , ],
    MARGIN = c(1, 3), FUN = mean
  )) # matrix over time, mean over particles
  tmp_sd <- t(apply(
    output[inf_flux_idx, , ],
    MARGIN = c(1, 3), FUN = sd
  )) # matrix over time, SD over particles
  inf_flux <- matrix(
    tail(tmp, n = 1),
    ncol = args$patch_dims, nrow = args$patch_dims
  )
  inf_flux_sd <- matrix(
    tail(tmp_sd, n = 1),
    ncol = args$patch_dims, nrow = args$patch_dims
  )
  inf_flux_sd <- inf_flux_sd / sum(inf_flux) # normalize
  inf_flux <- inf_flux / sum(inf_flux) # normalize

  ## get notification flux
  tmp <- t(apply(
    output[note_flux_idx, , ],
    MARGIN = c(1, 3), FUN = mean
  )) # matrix over time, mean over particles
  note_flux <- matrix(
    tail(tmp, n = 1), # last for cumulative
    ncol = args$patch_dims, nrow = args$patch_dims
  )
  tmp_sd <- t(apply(
    output[note_flux_idx, , ],
    MARGIN = c(1, 3), FUN = sd
  )) # matrix over time, SD over particles
  note_flux_sd <- matrix(
    tail(tmp_sd, n = 1), # last for cumulative
    ncol = args$patch_dims,
    nrow = args$patch_dims
  )
  note_flux_sd <- note_flux_sd / sum(note_flux) #normalize
  note_flux <- note_flux / sum(note_flux) # normalize

  return(list(
    note_flux = note_flux,
    note_flux_sd = note_flux_sd,
    inf_flux = inf_flux,
    inf_flux_sd = inf_flux_sd
  ))
}

## --- comparison plotter
make_flux_compare_graph <- function(CF) {
  ## FOI vs Note
  w <- 1e-3
  GP1 <- ggplot(
    CF,
    aes(
      x = output.FOI.flux,
      xmin = output.FOI.flux - 1.96 * output.FOI.flux.sd,
      xmax = output.FOI.flux + 1.96 * output.FOI.flux.sd,
      y = output.note.flux,
      ymin = output.note.flux - 1.96 * output.note.flux.sd,
      ymax = output.note.flux + 1.96 * output.note.flux.sd,
      label = variable, shape = `source zone`, col = `source zone`
    )
  ) +
    geom_abline(intercept = 0, slope = 1, col = 2, lty = 2) +
    scale_shape_manual(values = 0:6) +
    geom_errorbar(width = w, alpha = 0.4) +
    geom_errorbarh(width = w, alpha = 0.4) +
    geom_point(size = 2) +
    geom_text_repel(show.legend = FALSE) +
    theme_classic() +
    ggpubr::grids() +
    xlab("Model infection flux") +
    ylab("Model notification flux") +
    theme(legend.position = "top") +
    paletteer::scale_colour_paletteer_d("ggthemes::calc") +
    guides(shape = guide_legend(nrow = 1, byrow = TRUE))
   ## GP1


  ## data vs note
  GP2 <- ggplot(
    CF,
    aes(
      x = genomic.flux,
      y = output.note.flux,
      xmin = genomic.flux - 1.96 * genomic.flux.sd,
      xmax = genomic.flux + 1.96 * genomic.flux.sd,
      ymin = output.note.flux - 1.96 * output.note.flux.sd,
      ymax = output.note.flux + 1.96 * output.note.flux.sd,
      label = variable, shape = `source zone`, col = `source zone`
    )
  ) +
    geom_abline(intercept = 0, slope = 1, col = 2, lty = 2) +
    scale_shape_manual(values = 0:6) +
    geom_errorbar(width = w, alpha = 0.4) +
    geom_errorbarh(width = w, alpha = 0.4) +
    geom_point(size = 2) +
    geom_text_repel(show.legend = FALSE) +
    theme_classic() +
    ggpubr::grids() +
    paletteer::scale_colour_paletteer_d("ggthemes::calc") +
    xlab("Genomic-derived flux") +
    ylab("Model notification flux") +
    theme(legend.position = "top") +
    guides(shape = guide_legend(nrow = 1, byrow = TRUE))
   ## GP2

  ## combine plot
  GP <- ggarrange(
    GP1, GP2,
    nrow = 1, ncol = 2,
    common.legend = TRUE,
    labels = c("A", "B")
  )
  GP
}

get_fluxes_from_MM <- function(MM, n_particles = 10) {
  ## modify arguments
  args$MM <- MM

  ## run model
  output <- run.model(
    args,
    0:end,
    n.particles = n_particles,
    convert = FALSE
  )

  ## get fluxes
  flx <- get_fluxes(output)
  return(flx)
}

get_single_results <- function(n, MM = NULL) {
  if (is.null(MM)) {
    ## create random mixing matrix
    MM <- get_random_MM()
  }

  ## get fluxes
  flx <- get_fluxes_from_MM(MM)

  ## comparison
  CF <- data.table(
    iter = n,
    genomic.flux = c(t(FM) / sum(FM)), # transpose
    output.FOI.flux = c(flx$inf_flux),
    output.note.flux = c(flx$note_flux)
  )
  CF[, variable := gsub("cum_inf_flux", "", nmz)]
  CF[, `recipient zone` := gsub("\\[", "", gsub(",.\\]", "", variable))]
  CF[, `source zone` := gsub("\\]", "", gsub("\\[.,", "", variable))]

  return(CF)
}



make_MM_from_x <- function(x) {
  ## create random mixing matrix
  ## MM <- matrix(exp(c(0, x)), nrow = 7, ncol = 7)
  MM <- matrix(x, nrow = 7, ncol = 7)
  return(MM)
}

get_error_from_x <- function(x) {
  MM <- make_MM_from_x(x)

  ## get fluxes
  flx <- get_fluxes_from_MM(MM, 20)

  ## comparison
  genomic.flux <- c(t(FM) / sum(FM)) # transpose
  offset <- mean(genomic.flux)
  genomic.flux <- genomic.flux + offset
  output.note.flux <- c(flx$inf_flux) + offset

  ## calculate error
  err <- sum((output.note.flux / genomic.flux - 1)^2)
  return(err)
}

get_ratio_from_x <- function(x, n = 20) {
  MM <- make_MM_from_x(x)

  ## get fluxes
  flx <- get_fluxes_from_MM(MM, n)

  ## comparison
  genomic.flux <- c(t(FM) / sum(FM)) # transpose
  offset <- mean(genomic.flux)
  genomic.flux <- genomic.flux + offset
  output.note.flux <- c(flx$inf_flux) + offset

  ## calculate ratio
  ratio <- output.note.flux / genomic.flux
  return(ratio)
}

get_random_MM <- function() {
  ## create random mixing matrix
  MM <- FM <- matrix(runif(49), nrow = 7)

  ## sample as assortativity + random
  FM <- FM / sum(FM) # normalize
  for (i in 1:nrow(MM)) { # loops for transparency
    for (j in 1:nrow(MM)) {
      MM[i, j] <- FM[i, j] / (pops[i] * prev[j])
    }
  }
  MM <- MM / max(Re(eigen(MM)$values)) # renormalize
  return(MM)
}
