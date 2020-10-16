## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/results/data_290917/Figures/Triplet_upset/triplet_overlaps/TRIPLET_upset.R

## Upset plots
## Overlap numbers are from Intervene tool: https://github.com/asntech/intervene

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(UpSetR)

expressionInput <- c('HDAC2'=5026,'HDAC2&H3K27ac'=12617,'H3K27ac'=29728,'MYC&HDAC2&H3K27ac'=16966,'MYC&HDAC2'=1691,'MYC'=1747,'MYC&H3K27ac'=1025)

## metadata
met = data.frame(PROTEIN = c("H3K27ac", "MYC", "HDAC2"),
                 COL = c("h3k27ac", "myc", "hdac2"))

pdf("../figures/PRIMARY_UPSET_PLOT_CUSTOM.pdf", height = 7, width = 10)

upset(fromExpression(expressionInput),
      nsets=4,
      nintersects=30,
      show.numbers="yes",
      #main.bar.color="#ea5d4e",
      #sets.bar.color="#317eab",
      empty.intersections=NULL,
      order.by = "freq",
      number.angles = 0,
      mainbar.y.label ="number of uniqe/overlaping binding sites",
      sets.x.label ="Total binding sites",
      matrix.dot.alpha = 0,
      sets = c("H3K27ac", "HDAC2", "MYC"),
      keep.order = TRUE,
      point.size = 8,
      line.size = 2,
      text.scale = 2,
      set_size.angles = 0,
      matrix.color = "gray50",
      main.bar.color = "gray50",
      sets.bar.color = "gray50",
      queries = list(list(query = intersects, params = list("MYC", "HDAC2", "H3K27ac"), active = T, color = "black")),
      set.metadata = list(data = met,
                          plots = list(list(type = "matrix_rows",
                                            column = "COL",
                                            colors = c(h3k27ac = "#F22C1E", hdac2 = "#AD0AFD", myc="#78CA20"),
                                            alpha = 0.9))))


upset(fromExpression(expressionInput),
      nsets=4,
      nintersects=30,
      show.numbers="yes",
      #main.bar.color="#ea5d4e",
      #sets.bar.color="#317eab",
      empty.intersections=NULL,
      order.by = "freq",
      number.angles = 0,
      mainbar.y.label ="uniqe/overlaping binding sites",
      sets.x.label ="Total binding sites",
      matrix.dot.alpha = 0,
      sets = c("H3K27ac", "HDAC2", "MYC"),
      keep.order = TRUE,
      point.size = 8,
      line.size = 3,
      text.scale = 2,
      set_size.angles = 90,
      matrix.color = "black",
      main.bar.color = "black",
      sets.bar.color = "black",
      #text.scale = c(2, 0.5, 0, 0.5, 6, 6),
      #queries = list(list(query = intersects, params = list("MYC", "HDAC2", "H3K27ac"), active = T)),
      set.metadata = list(data = met,
                          plots = list(list(type = "matrix_rows",
                                            column = "COL",
                                            colors = c(h3k27ac = "#F22C1E", hdac2 = "#AD0AFD", myc="#78CA20"),
                                            alpha = 0.9))))

dev.off()


upset(fromExpression(expressionInput),
      nsets=4,
      nintersects=30,
      show.numbers="yes",
      #main.bar.color="#ea5d4e",
      #sets.bar.color="#317eab",
      empty.intersections=NULL,
      order.by = "freq",
      number.angles = 0,
      mainbar.y.label ="uniqe/overlaping binding sites",
      sets.x.label ="Total binding sites",
      matrix.dot.alpha = 0,
      sets = c("H3K27ac", "HDAC2", "MYC"),
      keep.order = TRUE,
      point.size = 8,
      line.size = 3,
      text.scale = 2,
      set_size.angles = 90,
      matrix.color = "black",
      main.bar.color = "black",
      sets.bar.color = "black",
      #text.scale = c(2, 0.5, 0, 0.5, 6, 6),
      #queries = list(list(query = intersects, params = list("MYC", "HDAC2", "H3K27ac"), active = T)),
      set.metadata = list(data = met,
                          plots = list(list(type = "matrix_rows",
                                            column = "COL",
                                            colors = c(h3k27ac = "#F22C1E", hdac2 = "#AD0AFD", myc="#78CA20"),
                                            alpha = 0.9))))


