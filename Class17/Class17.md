# Class17
Yuntian Zhu (PID: A17816597)

- [Section 1. Proportion of G/G in a
  population](#section-1-proportion-of-gg-in-a-population)
- [Section 2. Population scale analysis
  (Homework)](#section-2-population-scale-analysis-homework)

## Section 1. Proportion of G/G in a population

Here we start by reading the downloaded CSV file.

``` r
mxl <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(mxl)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
table(mxl$Genotype..forward.strand.)
```


    A|A A|G G|A G|G 
     22  21  12   9 

``` r
table(mxl$Genotype..forward.strand.)/nrow(mxl)
```


         A|A      A|G      G|A      G|G 
    0.343750 0.328125 0.187500 0.140625 

Now, let’s look at a different population, GBR.

``` r
gbr <- read.csv("373522-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(gbr)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  HG00096 (M)                       A|A ALL, EUR, GBR      -
    2                  HG00097 (F)                       G|A ALL, EUR, GBR      -
    3                  HG00099 (F)                       G|G ALL, EUR, GBR      -
    4                  HG00100 (F)                       A|A ALL, EUR, GBR      -
    5                  HG00101 (M)                       A|A ALL, EUR, GBR      -
    6                  HG00102 (F)                       A|A ALL, EUR, GBR      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
table(gbr$Genotype..forward.strand.)/nrow(gbr)
```


          A|A       A|G       G|A       G|G 
    0.2527473 0.1868132 0.2637363 0.2967033 

This variant that is associated with childhood asthma is more frequent
in the GBR population than the MKL population.

Lets now dig into this further.

## Section 2. Population scale analysis (Homework)

Now, let’s start some population scale analysis

> Q13: Read this file into R and determine the sample size for each
> genotype and their corresponding median expression levels for each of
> these genotypes.

``` r
expr <- read.table("rs8067378_ENSG00000172057.6.txt")
head(expr)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
table(expr$geno)
```


    A/A A/G G/G 
    108 233 121 

Therefore, we have 108 A/A individuals, 233 A/G individuals and 121 G/G
individuals.

``` r
tapply(expr$exp, expr$geno, median)
```

         A/A      A/G      G/G 
    31.24847 25.06486 20.07363 

Therefore, the median expression level for A/A is 31.24847, for A/G is
25.06486 and for G/G is 20.07363.

> Q14: Generate a boxplot with a box per genotype, what could you infer
> from the relative expression value between A/A and G/G displayed in
> this plot? Does the SNP effect the expression of ORMDL3?

``` r
library(ggplot2)

ggplot(expr) + aes(geno, exp, fill = geno) +
  geom_boxplot(notch = TRUE)
```

![](Class17_files/figure-commonmark/unnamed-chunk-9-1.png)

The boxplot shows a clear, graded difference in ORMDL3 expression across
the three genotypes. Individuals with the A/A genotype have the highest
median expression, A/G individuals show an intermediate median and G/G
individuals have the lowest median expression. The notches in the
boxplot do not overlap substantially between A/A and G/G, indicating
that their medians are significantly different.

Therefore, the SNP appears to have a strong effect on ORMDL3 expression:
carrying the G allele is associated with decreased ORMDL3 expression,
while the A allele is associated with higher expression.
