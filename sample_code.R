## Load required packages
install.packages("pak")
pak::pkg_install("YuHuiDeakin/rabc")
library(rabc)
library(tidyverse)
library(readr)
library(remotes)

## Ensure that version 1.7.11.1 of xgboost is installed to avoid errors 
remotes::install_version('xgboost', '1.7.11.1')

## Run function update from package author
max_freq_amp_update <- function (data_vec) 
{
  if(sd(data_vec) == 0) {
    return(c(0,0,0))
  }
  freq_vec <- abs(fft(lm(as.numeric(data_vec) ~ c(1:length(data_vec)))$residuals))
  half_ind <- floor(length(data_vec)/2)
  ind <- min(which(freq_vec[1:half_ind] == max(freq_vec[1:half_ind])))
  results <- c(as.numeric(ind), max(freq_vec[1:half_ind]), 
               entropy::entropy(freq_vec[1:half_ind])^2/half_ind)
  return(results)
}

calculate_feature_freq_update <- function (df_raw = NULL, samp_freq, axis_num = 3) 
{
  if (is.null(df_raw)) {
    stop("Please provide a valid data.frame!")
  }
  if (!exists("samp_freq")) {
    stop("Provide a valid samp_freq")
  }
  if (!is.numeric(samp_freq)) {
    stop("Samp_freq should be a number")
  }
  if (axis_num == 3 & (ncol(df_raw) - 1)%%3 != 0) {
    stop("Number of data column (not including the label column) should be three multiples.")
  }
  if (axis_num == 2 & (ncol(df_raw) - 1)%%2 != 0) {
    stop("Number of data column (not including the label column) should be two multiples.")
  }
  if (is.unsorted(as.character(df_raw[, ncol(df_raw)]))) {
    warning("This function expect the input data sorted by function order_acc.")
  }
  row_num <- dim(df_raw)[[1]]
  col_num <- dim(df_raw)[[2]]
  val_range <- max(as.numeric(as.matrix(df_raw[, -col_num]))) - 
    min(as.numeric(as.matrix(df_raw[, -col_num])))
  df_raw <- as.data.frame(df_raw)
  if (axis_num == 3) {
    sub_x <- df_raw[, seq(from = 1, to = col_num - 1, by = 3)]
    sub_y <- df_raw[, seq(from = 2, to = col_num - 1, by = 3)]
    sub_z <- df_raw[, seq(from = 3, to = col_num - 1, by = 3)]
    matsub_x <- as.matrix(sub_x)
    matsub_y <- as.matrix(sub_y)
    matsub_z <- as.matrix(sub_z)
    freq_index <- samp_freq/length(sub_x)
    freq_x <- apply(matsub_x, 1, max_freq_amp_update)
    freq_y <- apply(matsub_y, 1, max_freq_amp_update)
    freq_z <- apply(matsub_z, 1, max_freq_amp_update)
    x_freqmain <- freq_x[1, ] * freq_index
    y_freqmain <- freq_y[1, ] * freq_index
    z_freqmain <- freq_z[1, ] * freq_index
    x_freqamp <- freq_x[2, ]
    y_freqamp <- freq_y[2, ]
    z_freqamp <- freq_z[2, ]
    x_entropy <- freq_x[3, ]
    y_entropy <- freq_y[3, ]
    z_entropy <- freq_z[3, ]
    df_feature <- data.frame(x_freqmain, y_freqmain, z_freqmain, 
                             x_freqamp, y_freqamp, z_freqamp, x_entropy, y_entropy, 
                             z_entropy)
  }
  else if (axis_num == 2) {
    sub_x <- df_raw[, seq(from = 1, to = col_num - 1, by = 2)]
    sub_y <- df_raw[, seq(from = 2, to = col_num - 1, by = 2)]
    matsub_x <- as.matrix(sub_x)
    matsub_y <- as.matrix(sub_y)
    freq_index <- samp_freq/length(sub_x)
    freq_x <- apply(matsub_x, 1, max_freq_amp_update)
    freq_y <- apply(matsub_y, 1, max_freq_amp_update)
    x_freqmain <- freq_x[1, ] * freq_index
    y_freqmain <- freq_y[1, ] * freq_index
    x_freqamp <- freq_x[2, ]
    y_freqamp <- freq_y[2, ]
    x_entropy <- freq_x[3, ]
    y_entropy <- freq_y[3, ]
    df_feature <- data.frame(x_freqmain, y_freqmain, x_freqamp, 
                             y_freqamp, x_entropy, y_entropy)
  }
  else if (axis_num == 1) {
    sub_x <- df_raw[, seq(from = 1, to = col_num - 1, by = 1)]
    matsub_x <- as.matrix(sub_x)
    freq_index <- samp_freq/length(sub_x)
    freq_x <- apply(matsub_x, 1, max_freq_amp_update)
    x_freqmain <- freq_x[1, ] * freq_index
    x_freqamp <- freq_x[2, ]
    x_entropy <- freq_x[3, ]
    df_feature <- data.frame(x_freqmain, x_freqamp, x_entropy)
  }
  else {
    stop("Please provide valid number of axis from 1, 2, or 3.")
  }
  return(df_feature)
}

## Sort the data in order of activity type
accelo_sample_data <- read_csv("sample_dataset.csv")
accelo_data_sort <- order_acc(df_raw = accelo_sample_data)

## Plot the data, examine the patterns in behaviors 
plot_acc(df_raw = accelo_data_sort, axis_num = 3)

## Calculate the time domain features 
df_time <- calculate_feature_time(df_raw = accelo_data_sort, winlen_dba = 10)

## Calculate the frequency domain features
df_freq <- calculate_feature_freq_update(df_raw = accelo_data_sort, samp_freq = 40)

## Assign behavior labels to own object
labels <- accelo_data_sort$labels

## Examine a subset of time and frequency features for classifiers and view the plotted outcomes to see most relevant classification features
selection <- select_features(df_feature = cbind(df_time, df_freq), filter = FALSE, cutoff = 0.9, vec_label = labels, no_features = 10)

plot_selection_accuracy(results = selection)

## Use the provided feature visualizations to assess how well the selected features classify the behaviors 
# Indicate the feature you want to visualize in quotes

# Line plot
plot_feature(df_feature = df_time[, "y_min", drop = FALSE], vec_label = labels)

# Box plot
plot_grouped_feature(df_feature = df_time[, "y_min", drop = FALSE], vec_label = labels, geom = "boxplot")

# UMAP
plot_UMAP(df_time = df_time, df_freq = df_freq, label_vec = labels)

## Once you have assessed your selected features, create the classification model
class_model <- train_model(df = df_time[, selection$features[1:2]], vec_label = labels, train_ratio = 0.8)

## Assess model metrics using the confusion matrix
predictions <- plot_confusion_matrix(df_feature = df_time[, c(selection$features[1:2])], vec_label = labels)

## Assess incorrect classifications
plot_wrong_classifications(df_raw = accelo_data_sort, df_result = predictions)
