# Part 1: Basic information:

DOI: https://doi.org/10.5281/zenodo.12803445

This is part of the Archaeo-riddle project (https://theia.arch.cam.ac.uk/archaeoriddle/) answering "RQ1: What was the relationship between the two groups? Was it peaceful or hostile?".

Title: Using Point Process Modelling to detect cooperation vs competition

Author: Xuan Zhang (张璇), PhD student in School of Architecture, Tianjin University, China; visiting PhD student in Department of Archaeology, University of Cambridge, UK (2022-2023) 

The work was one of the proposals at the Archaeoriddle Workshop, EAA, 2023.


# Part 2: Contensts of the repository and the map of the folder/file structure

1.XZ-FHG.R: the R script;

2.Layers: the maps provided by the organizers;

3.Sites: the sites provided by the organizers and the sites gain from the extra 5 grids;

4.Output: the output figures;

5.Presentation_and_proposal: the methods and results of the proposal, and the EAA presentation slides.

Map of the folder/file structure:

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
│   ├── Methods_and_results_Xuan Zhang.doc
│   └── Xuan_Zhang_EAA_presentation_slides.pdf
├── README.md
└── Sites/
├── Biblio_data.csv
├── all_sites_for_R.csv
└── XZ-FHG.R
 ```

# Part 3: How to use the files/folders

If users wish to run this proposal on their own devices, they must download all files except for the README.md and the "Presentation_and_proposal" folder.  
**Step 1**: Open `XZ-FHG.R` in Neovim or RStudio.  
**Step 2**: Run `XZ-FHG.R` (Note: 1) The file paths have been standardized using the `here` package, so users do not need to adjust the paths; 2) A file named `all_sites_for_R.csv` is provided, which includes all data obtained after selecting five additional grids and processing chronological data via Oxcal online. Users may choose to skip these steps).  
**Step 3**: View the results in the `Output` folder.



# Part 4: Dependencies

1.I used R statistical computing language (R Core Team, 2023) through Neovim (Neovim, 2023).

2.I used OxCal Online to process the radiocarbon dates of each sites.

3.The versions of packages can be found at the beginning of XZ-FHG.R (the R script).



# Part 5: Answering the question: methods and results

This can be found in the Presentation_and_proposal\Methods_and_results_Xuan Zhang.doc.
