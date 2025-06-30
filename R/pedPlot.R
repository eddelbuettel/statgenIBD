#' Helper function for plotting the pedigree
#'
#' @importFrom utils tail head
#' @keywords internal
pedPlot <- function(pedigree,
                    offSpring,
                    popType,
                    genoCross,
                    multiCross = NULL,
                    title) {
  pedDatTot <- pedigree
  isDH <- endsWith(popType, "DH")
  isS <- grepl(pattern = "S|F", x = popType)
  isC3C4 <- popType %in% c("C3", "C4")
  isBC <- startsWith(popType, "BC")
  if (isBC) {
    nBC <- as.numeric(substring(text = popType, first = 3, last = 3))
  }
  ## Restrict to parents and first progeny.
  pedDatPar <- pedDatTot[!pedDatTot[["ID"]] %in% offSpring, ]
  pedDatOff <- pedDatTot[pedDatTot[["ID"]] %in% offSpring, ]
  pedDatOff <- pedDatOff[!duplicated(pedDatOff[c("par1", "par2")]), ]
  pedDatOff[["ID"]] <- if (!isC3C4) "F1" else popType
  if (!is.null(genoCross)) {
    pedDatTot[["cross"]] <- genoCross[["cross"]][match(pedDatTot[["ID"]],
                                                       genoCross[["geno"]])]
  } else {
    pedDatTot[pedDatTot[["ID"]] %in% offSpring, "cross"] <- "cross1"
  }
  pedDat <- rbind(pedDatPar, pedDatOff)
  ## Put most used parent central.
  parTab <- table(c(pedDat[["par1"]], pedDat[["par2"]]))
  ## First order by ID in pedigree to prevent crossing arrows as much as possible.
  parTab <- parTab[unique(pedigree[["ID"]])]
  parTab <- parTab[names(parTab) %in%
                     pedDat[pedDat[["type"]] == "INBPAR", "ID"]]
  if (length(parTab) > 1 && length(unique(parTab)) > 1) {
    parTab <- sort(parTab)
    parTab <- c(head(parTab, length(parTab) / 2), tail(parTab, 1),
                rev(head(rev(parTab)[-1], (length(parTab) - 1) / 2)))
    pedDat <- pedDat[c(match(names(parTab), table = pedDat[["ID"]], nomatch = 0),
                       (length(parTab) + 1):nrow(pedDat)), ]
  }
  generation <- as.numeric(factor(pedDat[["type"]],
                                  levels = unique(pedDat[["type"]]))) - 1
  ## Determine the row and column numbers in the plot.
  plotCols <- max(table(generation))
  plotRows <- length(unique(generation))
  if (!isC3C4) plotRows <- plotRows + 1
  ## Determine x and y positions of all parents and individuals in the plot
  xPos <- yPos <- NULL
  for (g in unique(generation)) {
    nGen <- sum(generation == g)
    if (nGen == plotCols) {
      xPosGen <- seq(1, plotCols, length.out = nGen)
    } else {
      xPosGen <- seq(1, plotCols, by = (plotCols - 1) / (nGen + 1))
      xPosGen <- xPosGen[-c(1, length(xPosGen))]
    }
    yPosGen <- rep(plotRows - g, nGen)
    xPos <- c(xPos, xPosGen)
    yPos <- c(yPos, yPosGen)
  }
  pedDat[["xPos"]] <- xPos
  pedDat[["yPos"]] <- yPos
  ## Construct data for plotting arrows.
  arrowDat <- rbind(merge(pedDat,
                          pedDat[, !colnames(pedDat) %in% c("par1", "par2")],
                          by.x = "par1", by.y = "ID", sort = FALSE),
                    merge(pedDat,
                          pedDat[, !colnames(pedDat) %in% c("par1", "par2")],
                          by.x = "par2", by.y = "ID", sort = FALSE))
  arrowDat <- arrowDat[order(arrowDat[["yPos.x"]], decreasing = TRUE), ]
  arrowDat[["linetype"]] <- "solid"
  ## Add extra arrow to bottom of plot.
  if (!isC3C4) {
    extArrow <- tail(arrowDat, sum(generation == max(generation)))
    extArrow[["linetype"]] <- "dotted"
    extArrow[["yPos.y"]] <- extArrow[["yPos.x"]]
    extArrow[["yPos.x"]] <- extArrow[["yPos.x"]] - 1
    extArrow[["xPos.y"]] <- extArrow[["xPos.x"]]
    arrowDat <- rbind(arrowDat, extArrow)
  } else {
    extArrow <- data.frame()
  }
  ## Construct data for labels.
  labDat <- arrowDat[colnames(arrowDat) != "ID"]
  labDat <- merge(labDat, pedDat[c("ID", "xPos", "yPos")],
                  by.x = c(c("xPos.y", "yPos.y")), by.y = c("xPos", "yPos"))
  extLab <- tail(arrowDat,
                 if (!is.null(genoCross)) length(unique(genoCross[["cross"]])) else 1)
  extLab[["yPos.y"]] <- extLab[["yPos.x"]]
  extLab[["xPos.y"]] <- extLab[["xPos.x"]]
  extLab[["ID"]] <- popType
  labDat <- rbind(labDat, extLab)
  ## Construct texts.
  ## Get number of individuals per cross.
  pedDatTot[["cross"]] <- factor(pedDatTot[["cross"]],
                                 levels = unique(pedDatTot[["cross"]]))
  crossSizes <- table(pedDatTot[["cross"]])
  textDat <- data.frame(xPos.y = c(0, 0, 0, extLab[["xPos.y"]]),
                        yPos.y = c(plotRows, 1, 0, rep(0, nrow(extLab))),
                        text = c("Parent:", "Population type:", "size:",
                                 crossSizes))
  if (nrow(extArrow) > 0) {
    extText <- tail(arrowDat, nrow(extArrow))
    ## Construct text label for last - dashed - line in plot.
    ## Exact text depends on population type.
    if (popType == "DH") {
      lastText <- "haploid"
    } else {
      lastText <- NULL
      if (isBC) {
        BCtxt <- if (nBC == 1) "F1" else paste0("BC", nBC - 1)
        lastText <- paste(BCtxt , "x", pedDat[2, "ID"])
      }
      if (isS) {
        lastText <- c(lastText, "selfing")
      }
      if (isDH) {
        lastText <- c(lastText, "DH")
      }
      lastText <- paste(lastText, collapse = " + ")
    }
    extText[["text"]] <- lastText
    extText[["xPos.y"]] <- (extText[["xPos.x"]] + extText[["xPos.y"]]) / 2
    extText[["yPos.y"]] <- (extText[["yPos.x"]] + extText[["yPos.y"]]) / 2
    textDat <- rbind(textDat, extText[c("xPos.y", "yPos.y", "text")])
  }
  ## For back crosses add an extra dotted line from parent 2 to bottom.
  ## Should be done after other computations are done to keep correct
  ## label positions.
  if (isBC) {
    extArrow2 <- extArrow
    extArrow2[["yPos.y"]] <- plotRows
    extArrow2[["xPos.y"]] <- max(arrowDat[["xPos.y"]])
    arrowDat <- rbind(arrowDat, extArrow2)
  }
  arrowDat[["linetype"]] <- factor(arrowDat[["linetype"]],
                                   levels = unique(arrowDat[["linetype"]]))
  ## Construct title.
  if (is.null(title)) {
    title <- "pedigree"
  }
  ## Plot pedigree.
  ## Segments, arrows, labels and text separately and an additional
  ## segment for the arrow in the explanation part.
  p <- ggplot2::ggplot(arrowDat,
                       ggplot2::aes(x = .data[["xPos.y"]],
                                    y = .data[["yPos.y"]])) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data[["xPos.x"]],
                                       yend = .data[["yPos.x"]],
                                       linetype = .data[["linetype"]]),
                          linewidth = 1, color = "blue") +
    ggplot2::geom_segment(
      ggplot2::aes(xend = (.data[["xPos.x"]] + .data[["xPos.y"]]) / 2,
                   yend = (.data[["yPos.x"]] + .data[["yPos.y"]]) / 2),
      linewidth = 1, color = "blue",
      data = arrowDat[arrowDat[["linetype"]] == "solid", ],
      arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"),
                             type = "closed")) +
    ggplot2::geom_label(ggplot2::aes(label = .data[["ID"]]),
                        data = labDat, fill = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data[["text"]]),
                       data = textDat) +
    ggplot2::geom_segment(x = 0, y = plotRows - 0.2, xend = 0, yend  = 1.2,
                          linewidth = 1, color = "blue",
                          arrow = ggplot2::arrow()) +
    ggplot2::xlim(-0.5, plotCols + 0.5) +
    ggplot2::ylim(-0.5, plotRows + 0.5) +
    ggplot2::labs(x = "", y = "", title = title) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none")
  return(p)
}
