### 2021/01/01  (FIRST R CODE OF 2021!)
### resizes the decisiontable icon raster data

library(raster)
library(tidyverse)

icons <- c(
"Icon_Commercial_Fisheries_Herring_Mackerel_Lobster",
"Icon_Environmental_Considerations",
"Icon_Groundfish_Fishery",
"Icon_Groundfish_Species",
"Icon_Herring_Fishery_Option1",
"Icon_Herring_Fishery_Option2",
"Icon_Lobster",
"Icon_Lobster_Fishery",
"Icon_Predator_Fisheries_Tuna_Haddock_Flatfish",
"Icon_Predator_Species_Tuna_Haddock_Flatfish",
"Icon_Primary_Production",
"Icon_Protected_Species",
"Icon_Protected_Species_and_Tourism",
"Icon_Tourism",
"Icon_Tuna",
"Icon_Tuna_Fishery")


icon <- icons[1]
rescale_icon <- function(icon) {
  load(paste0("data/",icon,".RData"))
  xx <- aggregate(raster(get(icon)), fact=10)
  get(icon)@grid@cells.dim
  assign(icon, as(xx, "SpatialPixelsDataFrame"))
  get(icon)@grid@cells.dim
  save(list=icon, file = paste0("data/",icon,".RData"), compress = "xz")
}
purrr::map(icons, rescale_icon)

