# CrediSite
Bayesian Inference of Amino Acid Level Geno-/Phenotype Correlation in Protein Multiple Sequence Alignments

# TODO:
- [ ] Make .rmd documentation
- [ ] Add descriptions to all elements in package
- [x] Run SigniSite benchmark
- [x] Discuss model configuration
  - [ ] Test model changes: 
    - [ ] Separate BLOSUM and positional prior
    - [ ] Test LKJ mix correlation matrix
    - [ ] BLOSUM across positions
- [ ] Add R-hat, etc. to paper?
- [ ] Redo regression
- [ ] Redo structure plots
- [ ] Redo density plot
- [x] Redo performance plots
- [x] Redo standard performance runs with new data
- [x] Transfer Docs to LaTeX template
- [ ] Reread and -write paper 
- [ ] Include example data


## Introduction

## Installation
This package requires rstan (follow the instructions for installation [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)).

Afterwards, installation of CrediSite can be done using devtools:
``` r
devtools::install_github("FSGade/CrediSite")
```

## Tutorial
``` r
input_data <- credisite::msa_to_model_data("msa/ATV.msa",
                                   consensus = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
                                   phenotype_transform = function(x){max(log10(x), -2)})

genotype <- input_data$genotype
phenotype <- input_data$phenotype

model <- credisite::fit_model(
            genotype = genotype,
            phenotype = phenotype,
            model_type = "regularised",
            cores = 4,
            iter = 3000,
            adapt_delta = 0.995,
            max_treedepth = 15
          )

outfilename <- "ATV_CrediSite_regularised.rds"
saveRDS(model, outfilename)
```

