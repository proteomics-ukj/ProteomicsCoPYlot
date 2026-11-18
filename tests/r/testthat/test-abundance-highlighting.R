test_that("abundance selection is finite, exact and combines both tails", {
  features <- c("AKT1_HUMAN;MTOR_HUMAN", "P~S1", "P~S2", "LOW", "BAD")
  values <- c(4, 3, 2, 1, Inf)
  expect_equal(which(.abundance_selection(features, values, c("mtor", "P~S1"), 1, 1)), c(1L, 2L, 4L))
  expect_false(any(.abundance_selection(features, values, "P")))
  expect_equal(which(.abundance_selection(features, values, "P~S1")), 2L)
  expect_equal(sum(.abundance_selection(features, values, NULL, 100, 100)), 4)
})

test_that("interactive abundance labels retain conditions and use GUI defaults", {
  skip_if_not_installed("plotly")
  data <- data.frame(ProteinNames=c("HIGH", "MID", "LOW"), a=c(100,10,1), b=c(1,10,100))
  meta <- data.frame(sample=c("a", "b"), condition=c("A", "B"))
  default <- plotly::plotly_build(interactive_abundance_plot(data, meta, "A", search="MID", data_scale="linear"))
  selected <- Filter(function(trace) grepl("Highlighted", trace$name) && trace$legendgroup != "_copylot_highlight_legend", default$x$data)
  expect_length(selected, 1)
  expect_equal(selected[[1]]$mode, "markers")
  expect_equal(as.numeric(selected[[1]]$x), 2)
  expect_equal(selected[[1]]$marker$symbol, "circle-open")
  expect_equal(selected[[1]]$marker$color, "#333333")
  expect_false(identical(selected[[1]]$marker$color, default$x$data[[1]]$marker$color))
  expect_equal(selected[[1]]$marker$size, 6)
  expect_true(all(vapply(default$x$data, function(trace) !grepl("text", trace$mode), logical(1))))
  for (outline in c(TRUE, FALSE)) for (labels in c(TRUE, FALSE)) {
    fig <- plotly::plotly_build(interactive_abundance_plot(data, meta, "All Conditions", data_scale="linear",
      top_high=1, top_low=1, search="MID", highlight_outline=outline, highlight_labels=labels))
    overlays <- Filter(function(trace) grepl("Highlighted", trace$name) && trace$legendgroup != "_copylot_highlight_legend", fig$x$data)
    expect_length(overlays, if (outline) 2 else 0)
    for (trace in overlays) {
      expect_equal(as.numeric(trace$x), c(1,2,3))
      expect_equal(as.numeric(trace$y), c(2,1,0))
      expect_equal(grepl("markers", trace$mode), outline)
      expect_false(grepl("text", trace$mode))
      expect_null(trace[["text"]])
    }
  }
})

test_that("abundance highlighting works for every notebook workflow and scale", {
  skip_if_not_installed("plotly")
  meta <- data.frame(sample=c("a", "b"), condition=c("A", "B"))
  for (workflow in c("Protein", "Phosphosite", "Phosphoprotein")) {
    data <- data.frame(a=c(100,10,1), b=c(1,10,100))
    data[[.workflow_col(workflow)]] <- c("HIGH", "MID", "LOW")
    for (scale in c("linear", "log2")) {
      input <- data
      if (scale == "log2") input[c("a", "b")] <- log2(input[c("a", "b")])
      fig <- plotly::plotly_build(interactive_abundance_plot(input, meta, "All Conditions", workflow,
        data_scale=scale, top_high=1, top_low=1, highlight_labels=TRUE))
      overlays <- Filter(function(trace) grepl("Highlighted", trace$name) && trace$legendgroup != "_copylot_highlight_legend", fig$x$data)
      expect_length(overlays, 2)
      expect_equal(as.numeric(overlays[[1]]$y), c(2, 0))
      expect_equal(as.numeric(overlays[[2]]$y), c(2, 0))
      path <- tempfile(fileext=".pdf")
      grDevices::pdf(path)
      expect_no_error(abundance_plot(input, meta, workflow, data_scale=scale,
                                     top_high=1, top_low=1, highlight_labels=TRUE))
      grDevices::dev.off()
      unlink(path)
    }
  }
})

test_that("mouse matching, gene labels and point-following rings agree", {
  expect_equal(.abundance_selection(c("GFRP_MOUSE", "GFRP2_MOUSE", "GFRP_MOUSE~S1"), c(3,2,1), "GFRP"), c(TRUE,FALSE,FALSE))
  d <- data.frame(ProteinNames=c("GFRP_MOUSE", "OTHER_MOUSE"), GeneNames=c("Gfrp", "Other"), a=c(100,10))
  m <- data.frame(sample="a", condition="A")
  for (size in c(4, 12)) {
    fig <- plotly::plotly_build(interactive_abundance_plot(d,m,"A",data_scale="linear",label_column="GeneNames",
      search="Gfrp",point_size=size,highlight_width=2.5,highlight_labels=TRUE))
    keys <- Filter(function(t) isTRUE(t$showlegend) && grepl("Highlighted",t$name), fig$x$data)
    expect_length(keys, 1)
    expect_equal(keys[[1]]$name, "Highlighted")
    expect_equal(fig$x$data[[1]]$marker$size, size)
    ring <- Filter(function(t) grepl("Highlighted",t$name) && t$legendgroup != "_copylot_highlight_legend", fig$x$data)[[1]]
    expect_equal(ring$marker$size,size)
    expect_equal(ring$marker$line$width,2.5)
    expect_equal(ring$mode,"markers")
    expect_match(ring$hovertext,"Gfrp")
  }
  v <- data.frame(ProteinNames=c("GFRP_MOUSE","OTHER_MOUSE"),log2FC=c(2,-2),neg_log10_selected_pval=c(3,2),significance=c("Upregulated","Downregulated"))
  fig <- plotly::plotly_build(volcano_plot_h(v,"A","B",highlight_terms="GFRP",top_down=1,point_size=11,highlight_width=2))
  keys <- Filter(function(t) isTRUE(t$showlegend) && grepl("Highlighted",t$name), fig$x$data)
  expect_length(keys, 1)
  expect_equal(keys[[1]]$name, "Highlighted")
  expect_equal(sum(vapply(Filter(function(t) grepl("highlighted", t$name),fig$x$data),function(t) length(t$x), integer(1))),2L)
  expect_true(all(vapply(Filter(function(t) t$legendgroup != "_copylot_highlight_legend", fig$x$data),function(t) t$marker$size==11,logical(1))))
})
