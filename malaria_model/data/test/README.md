!! Files in this directory are not excluded from being tracked on the repository. So please don't put any sensitive data in this directory.

# Files

## `regions_*.txt` (Mandatory)

A `tab`-delimited txt-based tabular file, columns of which contains the following information:

1. population size
2. longitude
3. latitude
4. transmission coefficient (mosquito bites)
5. presence/absence of malaria post
6. time malaria post opened
7. proportion of the pop that is migrant
8. proportion of resistant parasites

~~ 9. timing of MDA campaigns ~~
~~ 10. elevation ~~

9. proportion of forest goer population
10. type of location
    - "v", a village
    - "e", forest entry point
    - "f", poi in a forest

## `age_weights_vector_*.txt` (Mandatory)

A txt-based file containing only one column representing the weights of age groups. Used to initialise the population's age property.