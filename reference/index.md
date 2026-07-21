# Package index

## Data model

- [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
  : Create a moby telemetry dataset
- [`mobyMeta()`](https://miguelgandra.github.io/moby/reference/mobyMeta.md)
  : Retrieve mobyData metadata
- [`is_moby()`](https://miguelgandra.github.io/moby/reference/is_moby.md)
  : Check whether an object is a mobyData dataset
- [`rays`](https://miguelgandra.github.io/moby/reference/rays.md)
  [`rays_detections`](https://miguelgandra.github.io/moby/reference/rays.md)
  [`rays_tags`](https://miguelgandra.github.io/moby/reference/rays.md)
  [`rays_deployments`](https://miguelgandra.github.io/moby/reference/rays.md)
  : Toy acoustic-telemetry dataset (two rays in an MPA)

## Import, quality control & curation

- [`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md)
  : Import and harmonise acoustic detection data
- [`importDeployments()`](https://miguelgandra.github.io/moby/reference/importDeployments.md)
  : Import and harmonise receiver deployment metadata
- [`importTags()`](https://miguelgandra.github.io/moby/reference/importTags.md)
  : Import and harmonise tag / animal metadata
- [`checkDeployments()`](https://miguelgandra.github.io/moby/reference/checkDeployments.md)
  : Check receiver deployment metadata (quality control)
- [`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md)
  : Match detections to receiver deployments and back-fill metadata
- [`plotArray()`](https://miguelgandra.github.io/moby/reference/plotArray.md)
  : Map the receiver array
- [`assignAnimalIDs()`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
  : Assign animal IDs (and tagging dates) to detections from tag
  metadata
- [`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
  : Filter and quality-control acoustic detections
- [`plotMinLag()`](https://miguelgandra.github.io/moby/reference/plotMinLag.md)
  : Diagnose the min_lag false-detection threshold

## Temporal & spatial classification

- [`getTimeBins()`](https://miguelgandra.github.io/moby/reference/getTimeBins.md)
  : Assign time bins
- [`getDielPhase()`](https://miguelgandra.github.io/moby/reference/getDielPhase.md)
  : Estimate diel phase
- [`getSeason()`](https://miguelgandra.github.io/moby/reference/getSeason.md)
  : Estimate annual season
- [`getReprodPeriod()`](https://miguelgandra.github.io/moby/reference/getReprodPeriod.md)
  : Assign reproductive status
- [`getSunTimes()`](https://miguelgandra.github.io/moby/reference/getSunTimes.md)
  : Retrieve diel phase boundary times (hours)
- [`calculateCOAs()`](https://miguelgandra.github.io/moby/reference/calculateCOAs.md)
  : Calculate Centers of Activity (COAs)
- [`calculateStepDistances()`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
  : Estimate step distances and reconstruct movement tracks
- [`getTrajectories()`](https://miguelgandra.github.io/moby/reference/getTrajectories.md)
  : Extract movement trajectories from a calculateStepDistances() result
- [`interpolateDistances()`](https://miguelgandra.github.io/moby/reference/interpolateDistances.md)
  : Distance interpolation
- [`calculateLandDists()`](https://miguelgandra.github.io/moby/reference/calculateLandDists.md)
  : Calculate distance to nearest shore/land feature and related
  movement metrics.
- [`correctPositions()`](https://miguelgandra.github.io/moby/reference/correctPositions.md)
  : Relocate points on land to the nearest marine cell.
- [`createWideTable()`](https://miguelgandra.github.io/moby/reference/createWideTable.md)
  : Create detections table in wide format

## Summaries, residency & home range

- [`calculateResidency()`](https://miguelgandra.github.io/moby/reference/calculateResidency.md)
  : Calculate residency indices
- [`calculateVisits()`](https://miguelgandra.github.io/moby/reference/calculateVisits.md)
  : Extract residence events (visits)
- [`summaryTable()`](https://miguelgandra.github.io/moby/reference/summaryTable.md)
  : Generate summary table for tagged animals
- [`calculateROM()`](https://miguelgandra.github.io/moby/reference/calculateROM.md)
  : Calculate rates of movement (and total distance travelled)
- [`calculateLinearityIndex()`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md)
  : Calculate the linearity index of individual movements
- [`movementTable()`](https://miguelgandra.github.io/moby/reference/movementTable.md)
  : Create movement stats table
- [`calculateUDs()`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
  : Calculate animals' kernel utilization areas
- [`calculateUDOverlap()`](https://miguelgandra.github.io/moby/reference/calculateUDOverlap.md)
  : Pairwise overlap between utilization distributions
- [`plotMaps()`](https://miguelgandra.github.io/moby/reference/plotMaps.md)
  : Plot home-range maps with movement trajectories

## Network analysis

- [`is_mobyNetwork()`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
  [`networkEdges()`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
  [`networkNodes()`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
  [`networkType()`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
  : Inspect a moby network object
- [`plot(`*`<mobyNetwork>`*`)`](https://miguelgandra.github.io/moby/reference/plot.mobyNetwork.md)
  : Plot a moby network
- [`calculateAssociations()`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md)
  : Calculate pairwise spatiotemporal associations (co-occurrence
  network)
- [`randomizeAssociations()`](https://miguelgandra.github.io/moby/reference/randomizeAssociations.md)
  : Test association-network co-occurrences against a null model
- [`plotAssociations()`](https://miguelgandra.github.io/moby/reference/plotAssociations.md)
  : Network representation of pairwise overlaps.
- [`plotAssociationMatrix()`](https://miguelgandra.github.io/moby/reference/plotAssociationMatrix.md)
  : Plot an association significance / overlap matrix
- [`plotGroupSizeDistribution()`](https://miguelgandra.github.io/moby/reference/plotGroupSizeDistribution.md)
  : Plot the distribution of co-occurring group sizes
- [`calculateTransitions()`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md)
  : Build a movement (transition) network between locations
- [`randomizeTransitions()`](https://miguelgandra.github.io/moby/reference/randomizeTransitions.md)
  : Test movement-network transitions against a null model
- [`plotMovements()`](https://miguelgandra.github.io/moby/reference/plotMovements.md)
  : Plot a movement (transition) network on a map
- [`transitionsTable()`](https://miguelgandra.github.io/moby/reference/transitionsTable.md)
  : Summarise a movement network as a transitions table
- [`networkMetrics()`](https://miguelgandra.github.io/moby/reference/networkMetrics.md)
  : Compute node- and network-level metrics for a moby network

## Visualization

- [`plotAbacus()`](https://miguelgandra.github.io/moby/reference/plotAbacus.md)
  : Abacus plot
- [`plotDeployments()`](https://miguelgandra.github.io/moby/reference/plotDeployments.md)
  : Plot receiver operating periods (deployment timeline)
- [`plotActograms()`](https://miguelgandra.github.io/moby/reference/plotActograms.md)
  : Actogram plot (per-individual diel activity)
- [`plotChronogram()`](https://miguelgandra.github.io/moby/reference/plotChronogram.md)
  : Chronogram (hour x date heatmap)
- [`plotContours()`](https://miguelgandra.github.io/moby/reference/plotContours.md)
  : Contour plot (hour x date heatmap of a continuous variable)
- [`plotMetricComparison()`](https://miguelgandra.github.io/moby/reference/plotMetricComparison.md)
  : Compare metrics across the levels of a grouping factor
- [`plotStationStats()`](https://miguelgandra.github.io/moby/reference/plotStationStats.md)
  : Plot receiver-based statistics
- [`plotPeriodogram()`](https://miguelgandra.github.io/moby/reference/plotPeriodogram.md)
  : Plot detection periodograms (dominant rhythms)
- [`plotScalogram()`](https://miguelgandra.github.io/moby/reference/plotScalogram.md)
  : Plot wavelet scalograms (time-resolved rhythms)
- [`shadeSeasons()`](https://miguelgandra.github.io/moby/reference/shadeSeasons.md)
  : Return season boundaries for background shading.
