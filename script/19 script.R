library(usethis)
 
use_git_config(user.name  = "Yaqi-liu-2023",
               user.email = "rnr22dea@uea.ac.uk")

# Required packages
# Load the necessary libraries

library(tidymodels)
library(readxl)
library(tidyverse)
library(here)
library(performance)
library(skimr)
        
# Read the dataset into R
ch4 <- read_xlsx(here("data", "DNA methylation data.xlsm"), sheet = 1)

# Explore the first few rows
head(ch4)

skim(ch4)

names(ch4)
dim(ch4)
class(ch4)
ch5 <- select(
  # the data object
  .data = ch4,
  # the variables you want to select
  'Age', 'CpG 1 TET2', 'CpG 2 TET2', 'CpG 3 TET2', 'CpG 4 TET2', 'CpG GRIA2 1', 'CpG GRIA2 2', 'ASPA 1')

linear_model <- linear_reg() |> 
  set_engine("lm")  # We are using the "lm" engine for linear regression

age_recipe <- recipe(Age ~ ., data = ch5) |> 
  step_center(all_predictors())  # Centering the predictors

workflow_model <- workflow() |> 
  add_model(linear_model) |> 
  add_recipe(age_recipe)

fit_model <- fit(workflow_model, data = ch5)

fit_model_summary <- tidy(fit_model)
fit_model_summary

# Get predictions on the training data

predictions_lm <- augment(fit_model, new_data = ch5)

# Plot observed vs. predicted values
ggplot(data = ch5, aes(x = Age, y = predictions_lm$.pred)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Observed vs Predicted Bat Age")

mse <- mean((ch5$Age - predictions_lm$.pred)^2)

# A cleaner function for calculating MSE
mse_impl <- function(model, data, predictor) {
  augment(model, new_data = data) |> 
    mutate(squared_error = (.pred - {{predictor}})^2) |> 
    summarise(mse = mean(squared_error)) |> 
    pull(mse)
}

mse <-  mse_impl(fit_model, ch5, Age)

rmse <- sqrt(mse)
rmse

# Get more comprehensive model statistics
glance(fit_model)

fit_model |> 
  extract_fit_engine() |> 
  check_model()

# Split the data into 80% training and 20% testing
set.seed(123)  # For reproducibility
split <- initial_split(ch4, prop = 0.8)

# Extract training and testing sets
train_data <- training(split)
test_data <- testing(split)

# Check the dimensions of the splits
glimpse(train_data)
glimpse(test_data)

# Fit the model on the training data
lm_fit <- fit(workflow_model, data = train_data)

# View the model summary
tidy(lm_fit)

# Calculate performance metrics for linear regression
cat("RMSE is:", sqrt(mse_impl(lm_fit, test_data, Age)))

####

# Make predictions on the test set using the linear model
lm_predictions <- predict(lm_fit, new_data = test_data)

# Combine predictions with actual values (truth) into a tibble
results <- tibble(
  truth = test_data$Age,  # Actual values (truth)
  estimate = lm_predictions$.pred  # Predicted values (estimate)
)

# Now use yardstick's rsq function to calculate R-squared
library(yardstick)
rsq_result <- rsq(results, truth = truth, 
                  estimate = estimate)

# Print the R-squared result
rsq_result

# Perform 10-fold cross-validation
#v/k same thing

folds <- vfold_cv(train_data, v = 10)

# Fit models using cross-validation
cv_results <- fit_resamples(workflow_model, 
                            resamples = folds)

# Collect and visualize metrics
cv_metrics <- collect_metrics(cv_results)

ggplot(cv_metrics, aes(x = .metric, y = mean, color = .metric)) +
  geom_boxplot() +
  labs(title = "Cross-Validation Performance")

fold_numbers <- seq(2, 10, 1)

cv_results_by_fold <- map_dfr(fold_numbers, function(k_value) {
  # Create cross-validation with k folds
  k_folds <- vfold_cv(train_data, v = k_value)
  
  # Fit models using cross-validation
  cv_results <- fit_resamples(workflow_model, resamples = k_folds)
  
  # Collect metrics and add k value
  cv_metrics <- collect_metrics(cv_results)
  cv_metrics$.k <- k_value
  
  return(cv_metrics)
})

# Plot performance by number of folds
ggplot(cv_results_by_fold, aes(x = .k, y = mean, color = .metric)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ .metric, scales = "free_y") +
  labs(title = "Performance by Number of Cross-Validation Folds",
       x = "Number of Folds (k)",
       y = "Mean Performance")

# Use 5-fold cross-validation for final assessment
final_k <- vfold_cv(train_data, v = 5)

# Fit with prediction saving
final_cv_results <- fit_resamples(
  workflow_model, 
  resamples = final_k,
  control = control_resamples(save_pred = TRUE)
)

# Examine performance metrics
collect_metrics(final_cv_results)

# Get all predictions across folds for visualization
cv_predictions <- collect_predictions(final_cv_results)

# Visualize predictions vs actual values across all folds
ggplot(cv_predictions, aes(x = Age, y = .pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Cross-Validation Predictions vs Actual Age",
       x = "Actual Age (years)",
       y = "Predicted Age (years)")

# Define a Lasso regression model with a specific penalty
lasso_model <- linear_reg(penalty = 0.1, mixture = 1) |> 
  set_engine("glmnet")

# Create a workflow with the Lasso model
workflow_model_lasso <- workflow() |> 
  add_model(lasso_model) |> 
  add_recipe(age_recipe)

# Fit the Lasso model
fit_lasso <- fit(workflow_model_lasso, data = ch5)

# Check model performance
mse_impl(fit_lasso, ch5, Age)

# Define a KNN model with K=1
knn_model <- nearest_neighbor(neighbors = 1) |> 
  set_engine("kknn") |> 
  set_mode("regression")

# Create a workflow with the KNN model
workflow_model_knn <- workflow() |> 
  add_model(knn_model) |> 
  add_recipe(age_recipe)

# Fit the KNN model
fit_knn <- fit(workflow_model_knn, data = ch5)

# Check model performance
mse_impl(fit_knn, ch4, Age)

# Try different K values
k_values <- c(1, 5, 10, 25, 50, 100)

k_results <- map_dfr(k_values, function(k) {
  # Define KNN model with current K value
  knn_model <- nearest_neighbor(neighbors = k) |>  
    set_engine("kknn") |>  
    set_mode("regression")
  
  # Create workflow
  workflow_model_knn <- workflow() |>  
    add_model(knn_model) |>  
    add_recipe(age_recipe)
  
  # Fit model
  fit_knn <- fit(workflow_model_knn, data = ch4)
  
  # Calculate MSE
  mse_value <- mse_impl(fit_knn, ch5, Age)
  
  # Return results as a data frame row
  tibble(k = k, mse = mse_value, rmse = sqrt(mse_value))
})

# Plot results
ggplot(k_results, aes(x = k, y = rmse)) +
  geom_line() +
  geom_point() +
  labs(title = "KNN Performance by Number of Neighbors (K)",
       x = "Number of Neighbors (K)",
       y = "Root Mean Squared Error (years)") +
  theme_minimal()

# Make sure we've created train/test split first
set.seed(123)
split <- initial_split(ch5, prop = 0.8)
train_data <- training(split)
test_data <- testing(split)

# Train models on training data
fit_lm <- fit(workflow_model, data = train_data)
fit_ridge <- fit(workflow_model_ridge, data = train_data)
fit_lasso <- fit(workflow_model_lasso, data = train_data)
fit_knn <- fit(workflow_model_knn, data = train_data)

# Generate predictions on test data
predictions <- list(
  "Linear" = augment(fit_lm, new_data = test_data)$.pred,
  "Ridge" = augment(fit_ridge, new_data = test_data)$.pred,
  "Lasso" = augment(fit_lasso, new_data = test_data)$.pred,
  "KNN" = augment(fit_knn, new_data = test_data)$.pred
)

# Create visualization of predictions vs actual values
plots_list <- map2(predictions, names(predictions), function(preds, model_name) {
  ggplot(data = test_data, aes(x = Age, y = preds)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = model_name,
         x = "Actual Age (years)",
         y = "Predicted Age (years)") +
    theme_minimal()+
    scale_y_continuous(limits = c(0,9))
})

# Display plots in a grid
library(patchwork)
wrap_plots(plots_list)
