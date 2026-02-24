### ggseg atlas prep ###

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("cv_sample", "pheno_cat", "sig_fdr"))
}

# dk for cortical
dk_clean <- ggseg::dk
dk_clean$data <- dk_clean$data |>
  dplyr::filter(!is.na(region) & region != "corpus callosum")

dk_region_gap <- 0.4
dk_border_inset <- 0.4

dk_fill <- dk_clean
dk_fill$data <- dk_clean$data |>
  sf::st_as_sf() |>
  sf::st_buffer(dist = -dk_region_gap) |>
  tibble::as_tibble() |>
  sf::st_as_sf()

dk_border <- dk_clean
dk_border$data <- dk_clean$data |>
  sf::st_as_sf() |>
  sf::st_buffer(dist = -(dk_region_gap + dk_border_inset)) |>
  tibble::as_tibble() |>
  sf::st_as_sf()

# aseg for subcortical
aseg_region_gap <- 0
aseg_border_inset <- 1e-5

aseg_fill <- ggseg::aseg
aseg_fill$data <- ggseg::aseg$data |>
  sf::st_as_sf() |>
  sf::st_buffer(
    dist = 0,
    joinStyle = "MITRE",
    nQuadSegs = 500
  ) |>
  sf::st_buffer(0) |>
  tibble::as_tibble() |>
  sf::st_as_sf()

aseg_border <- ggseg::aseg
aseg_border$data <- ggseg::aseg$data |>
  sf::st_as_sf() |>
  sf::st_buffer(dist = -(aseg_region_gap + aseg_border_inset), nQuadSegs = 500) |>
  tibble::as_tibble() |>
  sf::st_as_sf()

# aseg for global
aseg_all_clean <- ggsegTissue::aseg_all
aseg_all_clean$data <- aseg_all_clean$data |>
  dplyr::filter(label != "Brain-Stem") |>
  dplyr::filter(side != "axial")

glob_region_gap <- 1e-10
glob_border_inset <- 1e-10

glob_fill <- aseg_all_clean
glob_fill$data <- aseg_all_clean$data |>
  sf::st_as_sf() |>
  sf::st_buffer(dist = -glob_region_gap) |>
  tibble::as_tibble() |>
  sf::st_as_sf()

glob_border <- aseg_all_clean
glob_border$data <- aseg_all_clean$data |>
  sf::st_as_sf() |>
  sf::st_buffer(dist = -(glob_region_gap + glob_border_inset)) |>
  tibble::as_tibble() |>
  sf::st_as_sf()

### plotting helpers ###

plot_cortex <- function(df, fill_var, show_cv_facet_labels = FALSE, margin_top=5) {
  p <- df |>
    dplyr::filter(
      pheno_cat == "Regional Vol" |
        pheno_cat == "Regional SA" |
        pheno_cat == "Regional CT"
    ) |>
    dplyr::group_by(cv_sample, pheno_cat) |>
    ggplot2::ggplot() +
    ggplot2::theme_linedraw(base_size = 12) +
    ggseg::geom_brain(
      atlas = dk_fill,
      position = ggseg::position_brain(hemi ~ side),
      ggplot2::aes(fill = {{ fill_var }}),
      color = NA,
      size = 0
    ) +
    ggseg::geom_brain(
      atlas = dk_border,
      position = ggseg::position_brain(hemi ~ side),
      ggplot2::aes(color = sig_fdr, size = sig_fdr),
      fill = NA
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.x = if (show_cv_facet_labels) ggplot2::element_text() else ggplot2::element_blank(),
      plot.margin = ggplot2::margin(margin_top, 0, 5, 0),
      legend.position = "none"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "gray", "TRUE" = "black"),
      na.value = NA,
      guide = "none"
    ) +
    ggplot2::scale_size_manual(
      values = c("FALSE" = .2, "TRUE" = .2),
      na.value = 0,
      guide = "none"
    )

  if (show_cv_facet_labels) {
    p <- p + ggplot2::facet_grid(pheno_cat ~ cv_sample, labeller = ggplot2::labeller(cv_sample = c("A" = "Split A", "B" = "Split B")))
  } else {
    p <- p + ggplot2::facet_grid(pheno_cat ~ cv_sample)
  }
  p
}

plot_subcortex <- function(df, fill_var) {
  p <- df |>
    dplyr::filter(pheno_cat == "Subcortical Vol") |>
    dplyr::group_by(cv_sample, pheno_cat) |>
    ggplot2::ggplot() +
    ggplot2::theme_linedraw(base_size = 12) +
    ggseg::geom_brain(
      atlas = aseg_fill,
      ggplot2::aes(fill = {{ fill_var }}),
      color = NA,
      size = 0
    ) +
    ggseg::geom_brain(
      atlas = aseg_border,
      ggplot2::aes(color = sig_fdr, size = sig_fdr),
      fill = NA
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "none"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "gray", "TRUE" = "black"),
      na.value = NA,
      guide = "none"
    ) +
    ggplot2::scale_size_manual(
      values = c("FALSE" = .2, "TRUE" = .2),
      na.value = 0,
      guide = "none"
    ) +
    ggplot2::facet_grid(pheno_cat ~ cv_sample)

  p
}

plot_global <- function(df, fill_var) {
  df |>
    dplyr::select(-dplyr::any_of(c("geometry", "tissue_class"))) |>
    dplyr::group_by(cv_sample, pheno_cat) |>
    dplyr::filter(pheno_cat == "Global Vol") |>
    ggplot2::ggplot() +
    ggplot2::theme_linedraw(base_size = 12) +
    ggseg::geom_brain(
      atlas = aseg_all_clean,
      ggplot2::aes(fill = {{ fill_var }}),
      color = NA,
      size = 0
    ) +
    ggseg::geom_brain(
      atlas = aseg_all_clean,
      ggplot2::aes(color = sig_fdr, size = sig_fdr),
      fill = NA
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "right",
      legend.box.margin = ggplot2::margin(0, 0, 0, 5)
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "gray", "TRUE" = "black"),
      na.value = NA,
      guide = "none"
    ) +
    ggplot2::scale_size_manual(
      values = c("FALSE" = .2, "TRUE" = .2),
      na.value = 0,
      guide = "none"
    ) +
    ggplot2::facet_grid(
      pheno_cat ~ cv_sample,
      labeller = ggplot2::labeller(cv_sample = c("A" = "Split A", "B" = "Split B"))
    )
}

format_fill_var <- function(df, fill_var, plt, fill_name, fill_color, fill_limits) {
  fill_var_name <- rlang::as_name(rlang::enquo(fill_var))
  is_categorical <- is.factor(df[[fill_var_name]]) || is.character(df[[fill_var_name]])

  if (is_categorical) {
    plt <- plt +
      paletteer::scale_fill_paletteer_d(
        fill_color,
        name = fill_name,
        na.value = "gray90",
        limits = fill_limits,
        drop = FALSE
      )
  } else {
    plt <- plt +
      paletteer::scale_fill_paletteer_c(
        fill_color,
        name = fill_name,
        limits = fill_limits,
        direction = -1,
        na.value = "gray90"
      )
  }
  return(plt)
}

plot_cv_brain <- function(df, fill_var, fill_name, fill_color, fill_limits, height_list=c(1.03, 4.5, .93), return_plt = TRUE, include_global = TRUE) {
  show_cv_labels <- !include_global

  cortex_plt <- plot_cortex(
    df, {{ fill_var }},
    show_cv_facet_labels = show_cv_labels,
    margin_top = if (include_global) 5 else 0
  )
  cortex_plt <- format_fill_var(df, {{ fill_var }}, cortex_plt, fill_name, fill_color, fill_limits)

  subcort_plt <- plot_subcortex(df, {{ fill_var }})
  subcort_plt <- format_fill_var(df, {{ fill_var }}, subcort_plt, fill_name, fill_color, fill_limits)

  if (include_global) {
    glob_plt <- plot_global(df, {{ fill_var }})
    glob_plt <- format_fill_var(df, {{ fill_var }}, glob_plt, fill_name, fill_color, fill_limits)
    legend <- cowplot::get_legend(glob_plt)
    glob_plt <- glob_plt + ggplot2::theme(legend.position = "none")
  } else {
    subcort_with_legend <- subcort_plt + ggplot2::theme(legend.position = "right", legend.box.margin = ggplot2::margin(0, 0, 0, 5))
    legend <- cowplot::get_legend(subcort_with_legend)
  }

  if (return_plt == TRUE) {
    if (include_global) {
      plts <- cowplot::plot_grid(
        glob_plt,
        cortex_plt,
        subcort_plt,
        ncol = 1,
        nrow = 3,
        rel_heights = height_list,
        align = "v",
        axis = "lr",
        greedy = TRUE
      )
    } else {
      plts <- cowplot::plot_grid(
        cortex_plt,
        subcort_plt,
        ncol = 1,
        nrow = 2,
        rel_heights = height_list,
        align = "v",
        axis = "lr",
        greedy = TRUE
      )
    }

    full_plt_obj <- cowplot::plot_grid(plts, legend, ncol = 2, rel_widths = c(3, .6))
    return(full_plt_obj)
  } else {
    if (include_global) {
      plt_list <- list(glob_plt, cortex_plt, subcort_plt, legend)
    } else {
      plt_list <- list(cortex_plt, subcort_plt, legend)
    }
    return(plt_list)
  }
}
