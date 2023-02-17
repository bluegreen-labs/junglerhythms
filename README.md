# The Jungle Rhythms workflow

This repository provides the basis for the Jungle Rhythms project's context, pre- and post-processing. However, this repository can also serve as a template for other data recovery project. You can cite the code in this repository like this:

> Hufkens K. and Kearsley (2023). The Jungle Rhythms workflow: recovering historical tropical tree phenology data - a citizen science project https://doi.org/10.5281/zenodo.XYZ

REFERENCE

## Introduction

Jungle Rhythms transcribed old observations of tree life cycle events (flowering, leaf shedding, fruit dispersion), which are key to understanding a tree's functioning. The African rainforest, the second largest on Earth, covers ~630 million ha and stores up to 66 Pg of carbon. The forest is presently a persistent carbon sink, offsetting large amounts of human CO2 emissions. Drought events in tropical rainforests have the potential to alter forest structure. However, due to data scarcity, little is known on how droughts affect the structure and function of African rainforests. The Jungle Rhythms data provides us with key information on how sensitive tree species are to drought, and how this sensitivity might alter the structure and function of the forest as drought regimes change.

### Historical observations

From 1937-1958 scientists stationed at the Yangambi research station in the Democratic Republic of the Congo made an effort to observe more than 2,000 trees on a weekly basis, writing down key life cycle events, e.g. fruit development, flowering and leaf shedding. All of these weekly observations were jotted down in little notebooks and finally summarized in large hand-drawn tables. In an effort to recover the key parts of this knowledge, which currently only exists on paper, and preserve the original copy, Jungle Rhythms will transcribe the summary tables. 

![]("https://raw.githubusercontent.com/khufkens/junglerhythms/master/vignettes/images/sheet.jpg")

In addition, these unique data will tell us how tropical forests respond to changing patterns in temperature and rain. As such, the data will allow us to predict the future state of the forest using historical data. A good summary of the scope and value of the data is given by a write up on research activities of among others [the Jungle Rhythms project in The Guardian](https://www.theguardian.com/environment/2017/sep/22/long-lost-congo-notebooks-shed-light-how-trees-react-to-climate-change).

### The recovered data

We present a unique analysis of historical tropical tree phenology data collected in the central Congo Basin (1937-1956), before large-scale impacts of human-induced climate change. A first analysis of the leaf phenological patterns of 129 species (4706 individual observation years, 91.2% basal area, 94 evergreen, 35 deciduous) illustrates the divergent behaviour within and across species, the variability of climate-phenology relationships, and its influence on landscape scale phenology. The full dataset will be published together with an upcoming data paper.

## Use

This R package contains all processing code and the raw data (in due time) of the Jungle Rhythms project. You can use the code in reprocessing. We also document parts of the workflow followed to provide people with inspiration to start their own projects. Note that this code comes AS IS. I refer to [our consulting policy](https://bluegreenlabs.org/labs/#support) if you need custom advice on your project.

To install the package run the following commands:

``` r
if(!require(devtools)){install.packages("devtools")}
devtools::install_github("bluegreen-labs/junglerhythms")
library("junglerhythms")
```

For general information we refer to the vignettes on the [documenation site](https://khufkens.github.io/junglerhythms/articles/).

# Acknowledgements

This project was started as a personal by Koen Hufkens and BlueGreen Labs, but over the years has been supported by the National Science Foundation’s Macro-system Biology Program (awards EF-1065029 and EF-1702697), the Belgian Science Policy office COBECORE project (BELSPO; grant BR/175/A3/COBECORE) and the Marie Skłodowska-Curie Action (H2020 grant 797668).
