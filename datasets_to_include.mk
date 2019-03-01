null  :=
space := $(null) #
comma := ,

## Define the data sets to include in the plots
datasetsms := DCS108 HEK293RNA RNA001 pilot
datasetsmsc := $(subst $(space),$(comma),$(datasetsms))

conditionsms := wt
conditionsmsc := $(subst $(space),$(comma),$(conditionsms))