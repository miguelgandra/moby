#######################################################################################################
## Internal helpers: colour palettes #################################################################
#######################################################################################################


##################################################################################################
## Okabe-Ito qualitative palette #################################################################
## Colourblind-safe qualitative palette (Okabe & Ito 2008), the de-facto standard for categorical
## scientific figures. Reordered to place high-contrast hues first (yellow/black deferred), which
## suits dense timeline/abacus plots drawn on light backgrounds.

#' Okabe-Ito colourblind-safe qualitative palette
#'
#' @description Returns `n` colours from the Okabe-Ito qualitative palette (Okabe & Ito 2008), the
#' de-facto colourblind-safe standard for categorical scientific figures. Colourblind-safety is a
#' property of the 8 base colours: for `n` greater than 8 the palette falls back to a spaced HCL
#' qualitative palette (`"Dark 3"`), which maximises separation but — like any categorical scheme —
#' cannot keep colours distinguishable under colour-vision deficiency beyond ~8 categories. (The older
#' behaviour of interpolating the 8 base colours was worse on both counts and was dropped.)
#' @param n Number of colours.
#' @return A character vector of hex colours.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.okabe_ito_pal <- function(n) {
  oi <- c("#E69F00", "#0072B2", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#F0E442", "#000000")
  if (n <= length(oi)) oi[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}


##################################################################################################
## Economist color palette #######################################################################
## Based in ggthemes::economist_pal()

#' Economist color palette
#'
#' @description This function generates a color palette using the Economist color scheme,
#' inspired by the color palette provided in the ggthemes::economist_pal() function.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.economist_pal <- function(n){
  economist_colors <- c("#6794a7", "#014d64", "#01a2d9", "#7ad2f6","#00887d",
                        "#76c0c1", "#7c260b", "#ee8f71", "#a18376", "#adadad")
  if(n==1) economist_colors[2]
  else if(n==2) return(economist_colors[c(3,2)])
  else if(n==3) return(economist_colors[c(1,2,3)])
  else if(n==4) return(economist_colors[c(1,2,3,9)])
  else if(n==5) return(economist_colors[c(1,2,4,3,6)])
  else if(n==6) return(economist_colors[c(1,2,4,3,6,5)])
  else if(n==7) return(economist_colors[c(1:6,9)])
  else if(n==8) return(economist_colors[c(1:6,8:9)])
  else if(n==9) return(economist_colors[1:9])
  else if(n==10) return(economist_colors[1:10])
  else stop("Invalid number of colors requested for the economist palette. Please choose a number between 1 and 10.")
}


##################################################################################################
## Viridis color palette #########################################################################
## Sourced from viridis::viridis(100)

#' Viridis color palette
#'
#' @description This function generates a color palette using the Viridis color scheme,
#' known for its perceptually uniform properties. The palette is adapted from the viridis package.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.viridis_pal <- function(n) {
  viridis_colors <- c(
    "#440154FF", "#450558FF", "#46085CFF", "#470D60FF", "#471063FF",
    "#481467FF", "#481769FF", "#481B6DFF", "#481E70FF", "#482173FF",
    "#482576FF", "#482878FF", "#472C7AFF", "#472F7CFF", "#46327EFF",
    "#453581FF", "#453882FF", "#443B84FF", "#433E85FF", "#424186FF",
    "#404587FF", "#3F4788FF", "#3E4A89FF", "#3D4D8AFF", "#3C508BFF",
    "#3B528BFF", "#39558CFF", "#38598CFF", "#375B8DFF", "#355E8DFF",
    "#34608DFF", "#33638DFF", "#32658EFF", "#31688EFF", "#2F6B8EFF",
    "#2E6D8EFF", "#2D708EFF", "#2C718EFF", "#2B748EFF", "#2A768EFF",
    "#29798EFF", "#287C8EFF", "#277E8EFF", "#26818EFF", "#26828EFF",
    "#25858EFF", "#24878EFF", "#238A8DFF", "#228D8DFF", "#218F8DFF",
    "#20928CFF", "#20938CFF", "#1F968BFF", "#1F998AFF", "#1E9B8AFF",
    "#1F9E89FF", "#1FA088FF", "#1FA287FF", "#20A486FF", "#22A785FF",
    "#24AA83FF", "#25AC82FF", "#28AE80FF", "#2BB07FFF", "#2EB37CFF",
    "#31B67BFF", "#35B779FF", "#39BA76FF", "#3DBC74FF", "#41BE71FF",
    "#47C06FFF", "#4CC26CFF", "#51C56AFF", "#56C667FF", "#5BC863FF",
    "#61CA60FF", "#67CC5CFF", "#6DCD59FF", "#73D056FF", "#78D152FF",
    "#7FD34EFF", "#85D54AFF", "#8CD646FF", "#92D741FF", "#99D83DFF",
    "#A0DA39FF", "#A7DB35FF", "#ADDC30FF", "#B4DE2CFF", "#BBDE28FF",
    "#C2DF23FF", "#C9E020FF", "#D0E11CFF", "#D7E219FF", "#DDE318FF",
    "#E4E419FF", "#EBE51AFF", "#F1E51DFF", "#F7E620FF", "#FDE725FF"
  )
  return(colorRampPalette(viridis_colors)(n))
}



##################################################################################################
## Palr bathy deep color palette #########################################################################
## Sourced from palr::bathy_deep_pal

#' Bathy deep color palette
#'
#' @description This function generates a color palette using the Bathy Deep color scheme,
#' sourced from palr::bathy_deep_pal function.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.bathy_deep_pal <- function (x, palette = FALSE, alpha = 1) {
  breaks <- c(-5500, seq(-5000, -1000, by = 1000), -500, 0)
  breaks <- seq(-5500, 0, length = 255)
  cols <- colorRampPalette(rgb(c(0,18,60,103,141,194,255), c(0,18,60,111,163,216,255),
                               c(0,26,85,135,173,216,255), maxColorValue=255))(256)
  hexalpha <- as.hexmode(round(255 * alpha))
  if (nchar(hexalpha) == 1L) hexalpha <- paste(rep(hexalpha, 2L), collapse = "")
  cols <- paste0(cols, hexalpha)
  if (palette)  return(list(breaks = breaks, cols = cols))
  if (missing(x)) return(colorRampPalette(cols))
  if (length(x) == 1L) {return(paste0(colorRampPalette(cols)(x), hexalpha))}
  else {return(cols[findInterval(x, breaks)])}
}
