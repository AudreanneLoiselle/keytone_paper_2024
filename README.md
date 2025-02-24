# keytone_paper_2024
Data repository for manuscript : Adapting the concept of functionally dominant species for observational data

# Description
This is the data repository for the files linked to the manuscript titles _Data repository for manuscript : Adapting the concept of functionally dominant species for observational data

## Sampling description

### Plants
Plants were sampled twice during the summer of 2018 using 30 quadrats of 0.5x0.5m in total. In addition, species identified during a search time of 5 to 20 minutes, depenting on site size, were also added to the inventory data.

### Bird
Birds were sampled once during the summer of 2019 using digical recorders installed near the shore of each wetland. Species were identified by listening to a random sample of 3 minute audio clips from 5AM to 9AM.

### Fish
Fish species were sampled trice during the summer of 2019 using minnow traps installed in hydrological structures in each site. The number of minnow traps was adjusted according to the amound of hydrological structure in a site.

## Dataset description

### plant_DCi.CSV
This file contains the plant species relative frequencies (columns) in each wetland site (lines).

### bird_DCi.CSV
This file contains the bird species relative frequencies (columns) in each wetland site (lines).

### plant_DCi.CSV
This file contains the fish species relative frequencies (columns) in each wetland site (lines).

## site_data.CSV
This file contains for each sampled wetland (lines) the wetland type (wet_type), as well as the plant species richness (plant_rich), bird species richness (bird_rich) and fish species richness (fish_rich), calculated using species number. NA values in the bird_rich column indicate sites were we couldn't get data due to equipment failure. Type value "alder" refers to alder swamps, which were dominated by _Alnus incana_ ssp. _rugosa_. Type value "peat" refers to site that were poor fens dominated by ericaceous vegetation and _Sphagnum_ mosses. Type value "ash" refers to ash swamps that were dominated by _Fraxinus nigra_ with tree cover reaching at least 30%.

## keystone_code.R
This file contains and example of the method we used to compute the CIi values for the plant dataset.
