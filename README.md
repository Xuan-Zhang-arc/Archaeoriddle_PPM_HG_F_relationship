# Part 1: Basic information:

DOI: https://doi.org/10.5281/zenodo.12803445

This is part of the Archaeo-riddle project (https://theia.arch.cam.ac.uk/archaeoriddle/) answering "RQ1: What was the relationship between the two groups? Was it peaceful or hostile?".

Title: Using Point Process Modelling to detect cooperation vs competition

Author: Xuan Zhang (张璇), PhD student in School of Architecture, Tianjin University, China; visiting PhD student in Department of Archaeology, University of Cambridge, UK (2022-2023) 

The work was one of the proposals at the Archaeoriddle Workshop, EAA, 2023.


# Part 2: Contensts of the repository and the map of the folder/file structure

1.`XZ-FHG.R`: the R script;

2.`Layers`: the maps provided by the organizers;

3.`Sites`: the sites provided by the organizers and the sites gain from the extra 5 grids;

4.`Output`: the output figures;

5.`Presentation_and_proposal`: the methods and results of the proposal, and the EAA presentation slides.

**Map of the folder/file structure**:

```
Xuan-Zhang-arc-Archaeoriddle_PPM_HG_F_relationship
├── here
├── Layers/
│   ├── east_narnia4x.tif
│   ├── east_narnia4x.tif.aux.xml
│   ├── resources.tiff
│   └── resources.tiff.aux.xml
├── Output/
│   ├── 1_select_5_more/
│   │   ├── predict(gam_FHG_map).jpg
│   │   ├── predict(gam_F_map).jpg
│   │   └── predict(gam_HG_map).jpg
│   └── 2_analysis/
│       ├── Kcross_early_FHG.jpg
│       ├── Kcross_late_FHG.jpg
│       ├── Kcross_middle_early_FHG.jpg
│       ├── Kcross_middle_late_FHG.jpg
│       ├── early_pm.jpg
│       ├── late_pm.jpg
│       ├── late_pm_Al.jpg
│       ├── middle_early_pm.jpg
│       └── middle_late_pm.jpg
├── Presentation_and_proposal/
│   ├── Methods_and_results_Xuan_Zhang.doc
│   └── Xuan_Zhang_EAA_presentation_slides.pdf
├── README.md
└── Sites/
├── Biblio_data.csv
├── all_sites_for_R.csv
└── XZ-FHG.R
 ```


# Part 3: How to use the files/folders

If users wish to run this proposal on their own devices, they must download all files except for the `README.md` and the `Presentation_and_proposal` folder.  

**Step 1**: Open `XZ-FHG.R` in Neovim or RStudio.  

**Step 2**: Run `XZ-FHG.R`.

(Note: 
1) The file paths have been standardized using the `here` package, so users do not need to adjust the paths.
2) A file named `all_sites_for_R.csv` is provided, which includes all data obtained after selecting five additional grids and processing the radiocarbon dates via Oxcal Online. Users may choose to skip these steps.
3) The maptools pacakge has retired and removed from the CRAN repository, so please visit 'https://cran.r-project.org/src/contrib/Archive/maptools/' to obtain former available versions and install the package. The users can also use other packages as suggested on 'https://cran.r-project.org/web/packages/maptools/index.html' as a substitution. However, switching to other packages might have an impact on the results.) 

**Step 3**: View the results in the `Output` folder.


# Part 4: Dependencies

1.I used `R` statistical computing language (R Core Team, 2023) through `Neovim` (Neovim, 2023).

2.I used `OxCal Online` to process the radiocarbon dates of each sites.

3.The versions of packages can be found at the beginning of `XZ-FHG.R` (the R script).


# Part 5: Answering the question: methods and results

This proposal was part of the Archaeo-riddle project and used Point Process Modeling (PPM) through “spatstat (Baddeley & Turner, 2005)” R package to analyze the relationship between the farmers and hunter-gatherers (“RQ1. What was the relationship between the two groups? Was it peaceful or hostile?”). 

The aims of the proposal were to answer RQ1 and test the performance of PPM. 

This proposal contains 4 parts: Preparations (selecting 5 more grids and defining the time unit), Methods (defining the proxy for the two groups’ relationship and fitting the models), Results and Conclusions (the answer to RQ1 and PPM’s performance).

**1.Preparations**

Each participants of the Archaeo-riddle project (https://theia.arch.cam.ac.uk/archaeoriddle/) was given “one map with height values and another one with values with probabilities of settlement according to environmental fitness”, information of some known sites including site names, coordinates, radiocarbon dates and cultural affiliations. Participants were also asked to select 5 more grids to get more sites.

**1) Selecting 5 more grids**

To select 5 more grids, I built first-order point process models (Bevan, 2020) with the given environmental maps for farmer sites, hunter-gatherer sites and all sites and got predictive maps. One grid was selected for farmer sites, one for hunter-gatherer sites and three for all sites. 

**2) Defining the time unit**

After processing the radiocarbon dates of each sites via OxCal Online, the chronology of the sites was divided into the following stages: 1) the early period (before 6000 BC, 48 sites); 2) the middle-early period (from 6000 BC to 5500 BC, 77 sites) and the middle-late period (from 5500 BC to 5000 BC, 47 sites); 3) the late period (after 5000 BC, 7 sites). 

**2.Methods**

**1) Defining the proxy for relationship**

The distance between the two groups was used as a proxy for their relationships. According to evolutionary ecology studies (Field 2003), the distance between farms and hunter-gatherers could be a proxy for their relationship. During each period, if farmer sites and hunter-gatherer sites were crowded with each other, they were likely to be hostile to each other. If they kept reasonable distances from each other, they might have had a peaceful relationship. 

**2) Fitting the models**

For each period, I applied Cross-type K-Function to test whether the two groups were clustered or not (Baddeley, Rubak, & Turner, 2015). If they were clustered with each other, I generated a Multitype Strauss Model combined with Maximum Pseudolikelihood to get the optimum value of irregular parameter R, or rather the interaction distance between the two groups of sites (Baddeley, Rubak, & Turner, 2015). If they were not clustered with each other, I generated an Area-interaction Model to distinguish a regular process from a Poisson process (Baddeley, Rubak, & Turner, 2015). 

**3.Results**

Seen from the results (see the Output file) of the Cross-type K-Function, the two groups were clustered in all the periods except the late period. The interaction distance R in the early period was 1.38. The R of the middle-early period was 1.24 and the R was 1.25 in the middle-late period. Both R values were lower than that in the early period, indicating a possibly more hostile relationship than the beginning. During the late period, the fitted interaction parameter eta was 1.38197e-70 showing a regular process and there were more farmer sites than hunter-gatherer sites left. 

**4.Conclusions**

**1) The answer to RQ1**

Overall, the farmers and the hunter-gatherers might have had some conflicts during the middle period and the farmers were the winners. 

**2) PPM’s performance**

As for PPM’s performance in answering the question, it had advantages in dealing with different-order interactions and analyzing multitype second-order interactions. However, the method was limited by the existing points when sampling. It also suffered from the bias caused by time averaging.


