# why-natural-selection

Code for the paper "Human capital mediates natural selection in 
contemporary humans" by Abdel Abdellaoui and David Hugh-Jones.

# The paper

The current version is available at https://github.com/hughjonesd/why-natural-selection/blob/master/why-natural-selection.pdf.

# BGA poster

The BGA poster is available [https://github.com/hughjonesd/why-natural-selection/blob/master/BGA poster.pdf](here).

# To download

From your command line, run:

```
git clone https://github.com/hughjonesd/import-ukbb-data.git
git clone https://github.com/hughjonesd/why-natural-selection.git
```

# The data

Data is stored separately. To get UK Biobank data you'll need to
make an application to UK Biobank (or join our application). Other data is
publicly available, so we can share it.

We also use weights from 

Alten, Sjoerd van, Benjamin W Domingue, Titus Galama, and Andries Marees. 2021. 
“The Effects of Demographic-Based Selection Bias on GWAS Results in the UK Biobank.” 

Thank you to Sjoerd and coauthors for their kindness in providing these!
You'll need their code to recreate the weights. 

# To run


1. Start R within the project directory. This should automatically
   download the `{renv}` package.

2. Run `renv::restore()` to download the appropriate R packages.

3. Edit directory locations in `_drake.R`.

4. Within R, run

```r
r_make()
```

This will take some time to run, depending on your machine. It will rebuild
the PDF `why-natural-selection.pdf`.
