# File structure

- src/ contains source code
-- functions.R contains functions used for simulation
-- simulation.R the main script for running the simulation
- docs/ contains documents
- data/ contains reference data
- results/ contains results from simulation

# Functions

## simulate_data(data, nrow, ncols)

### arguments
- data: data matrix to be used for reference distributions
- nrow, ncol: number of rows & columns for the data matrix to be returned

### returns
- matrix of nrow x ncol, containing simulated data

### comments
- use mvrnorm() function for simulating data

## simulate_missingness(data, mcar=0, mar=0, mnar=0, mnar.type="left")

### arguments
- data: data matrix where missingness will be introduced.
- mcar: percentage of MCAR missingness (0..1)
- mar: percentage of MAR missingness (0..1)
- mnar: percentage of MNAR missingness (0..1)
- mnar.type: type of MNAR missingness, either "left" or "right"

### returns
- data with given types & proportions of NA values introduced

## detect_missingness_type(data)

### arguments
- data: data matrix of the data to be analyzed for missingness

### returns
- vector with length = ncol(data), each element describing the type of detected missingness for each column of the data. Possible values "NONE", "MCAR", "MAR" or "MNAR"

## select_imputation_method(types)

### arguments
- type: vector with each element describing the type of missingness. Possible values "NONE", "MCAR", "MAR" or "MNAR"

### returns
- vector with length = length(types), each element describing the imputation method for each element of the types. Example values "min", "mean", "PPCA" or "KNN"

## impute(data, methods)

### arguments

- data: data matrix to be imputed
- methods: vector of imputation methods for each column in the matrix. Possible values "RF", "KNN"

### returns

- data matrix with the missing values imputed
