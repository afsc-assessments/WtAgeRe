# Get predicted weight-at-age

Read ADMB output and return predicted weight-at-age by year.

## Usage

``` r
fn_get_pred(file = c("wt"), dat = NULL, source = c("model"), age_range = NULL)
```

## Arguments

- file:

  Character vector of base filenames (no extension). Default is "wt".

- dat:

  Optional data list used to infer age range (indices 8 and 9).

- source:

  Character vector of model labels for each file.

- age_range:

  Optional numeric vector of ages to use as column names.

## Value

A data frame with year, age columns, and a source column.
