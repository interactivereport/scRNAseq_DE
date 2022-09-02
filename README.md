<div id="top"></div>
<!--
*** Used a trimmed-down version of the markdown template here: https://github.com/othneildrew/Best-README-Template/blob/master/BLANK_README.md
-->

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

<!-- PROJECT LOGO -->
<br />
<div align="left">

<h3 align="left">scDEpipelineR6</h3>

  <p align="left">
    R code for performing scRNA-seq DE and evaluation of DE method performance using simulated data
    <br />
    <br />
    ·
    <a href="https://github.com/mryals/scDEpipelineR6/issues">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#nebula-output-note">Nebula Output Note</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

This project aims to evaluate performance of differential gene expression methods using simulated single cell RNA-seq data.  

<p align="right">(<a href="#top">back to top</a>)</p>

### Built With

* [R6](https://r6.r-lib.org/)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

Welcome to the source code for performing simulations and DE analysis for scRNA-seq data.  Please see the markdown document for more information.  This repository is a work in progress so check back for more information.  Currently the package is designed so users can have a look at the code used to generate the simulation data for DE methods presented in the manuscript attached to this package.

### Prerequisites

It is recommended to use R/3.5 for this code.

### Installation

```r
library(devtools)
devtools::install_github("interactivereport/scRNAseq_DE")
```

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Usage

This package contains an introductory vignette that shows the usage of the DE method pipeline, and an example of the call to the DE benchmarking simulation function.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- NEBULA OUTPUT NOTE -->
## Nebula Output Note

When using scRNA-seq DE methods in the code pipeline_class.R, the output is an R data.frame for DE methods T-test, ANCOVA, U-test, EdgeR, limma-voom, limma-cell, DESeq2, glmmTMB, and MAST. The data.frame of DE results contains at minimum columns for the tested gene name, log2FC, pvalue, and adjusted p-value (FDR) for these DE methods.  
In contrast, the output of DE results using the DE method NEBULA in pipeline_class.R is a list.  This list contains two elements with names ‘res.tab’ and ‘res.ls’.  The list element named ‘res.tab’ is an R data.frame similar to the output of the other DE methods in pipeline_class.R, and contains columns for the tested gene name, log2FC, pvalue, and FDR. 
For convenience to the user, the ‘res.tab’ data.frame also contains additional information in four additional columns that are specific to the NEBULA DE method.  The first additional column is ‘algorithm’, which will output the NEBULA method used to test a gene.  If the NEBULA method was specified as ‘HL’, then all values in this column will be ‘HL’; otherwise, if the method used was set to ‘LN’, then the output will correspond to the NEBULA method used for a specific tested gene. The ‘convergence’ column corresponds to the convergence of the algorithm. Briefly, a value in the ‘convergence’ column of ‘-20’ or ‘-30’ could indicate a failure of the algorithm for the tested gene.   For more information on convergence, please refer to the NEBULA publication [https://www.nature.com/articles/s42003-021-02146-6](https://www.nature.com/articles/s42003-021-02146-6).  The ‘DE_quality_score’ column is equal to the ratio of the tested gene mean counts per sample divided by the estimated cell-level overdispersion.  The ‘DE_quality_score’ ratio is intended to be a metric that corresponds to the true positive reliability of a DEG.  The ‘DE_quality_indicator’ column is used to flag if the row value of the ‘DE_quality_score’ column is below 20, which may indicate an unreliable DEG result.
The second element in the NEBULA DE result list, named ‘res.ls’, contains the raw output list from running the NEBULA function.  For more information on the ‘res.ls’ list and the definition of the objects in that list, please refer to the NEBULA function source documentation [https://github.com/lhe17/nebula](https://github.com/lhe17/nebula).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

Feel free to modify source package and code as needed for your scRNA-seq simulation and DE studies within the permissions of the license.

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTACT -->
## Contact and acknowledgements

**Contributors:**  

Jake Gagnon, Matthew Ryals, Lira Pi, Liping Hou, Qingwen Wan

Contact: [https://www.biogen.com/en_us/contact.html](https://www.biogen.com/en_us/contact.html)

Project Link: [https://github.com/interactivereport/scRNAseq_DE](https://github.com/interactivereport/scRNAseq_DE)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
