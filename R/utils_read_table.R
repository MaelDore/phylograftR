
##### Internal function to read a table from an external file #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Read a table (.csv, .tsv, .xlsx, ...) as a data frame ####

#' @title Read a table from a file
#'
#' @description Internal function shared by [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] and
#'   [`phylograftR::import_taxonomic_table()`][import_taxonomic_table]. It checks the inputs, chooses the reader from the extension of the file, and
#'   returns the table as a data frame. All the columns are read as text, blank cells and `NA` are missing values, and names of the columns are kept as they are.
#'
#' @param file Character string. Path to the file.
#' @param sep Character string or `NULL`. Field separator of a text file. If `NULL`, chosen from the extension.
#' @param sheet Character string or integer. Sheet of an Excel file.
#' @param reader Function or `NULL`. Optional function to read the file, which takes precedence over the extension.
#' @param ... Further arguments passed to the reader.
#'
#' @importFrom utils read.table
#'
#' @noRd
#'
#' @return A data frame with at least one row.
#'

read_table_file <- function (file, sep = NULL, sheet = 1, reader = NULL, ...)
{
  ## 1.1/ Check the inputs ####

  if (!is.character(file) || length(file) != 1 || is.na(file)) { stop("'file' must be a single path", call. = FALSE) }
  if (!file.exists(file)) { stop("File not found: ", file, call. = FALSE) }
  if (!is.null(reader) && !is.function(reader)) { stop("'reader' must be a function that takes the path of the file and returns a data frame", call. = FALSE) }
  if (!is.null(sep) && (!is.character(sep) || length(sep) != 1 || is.na(sep) || nchar(sep) != 1)) { stop("'sep' must be a single character", call. = FALSE) }

  ## 1.2/ Read the table ####

  extension <- tolower(sub(".*\\.", "", basename(file)))

  if (!is.null(reader))
  {
    data <- reader(file, ...)
  } else if (extension %in% c("xlsx", "xls"))
  {
    if (!requireNamespace("readxl", quietly = TRUE))
    {
      stop("The package 'readxl' is required to read Excel files. Install it with install.packages(\"readxl\"), ",
           "or export the table as a .csv file.", call. = FALSE)
    }
    data <- readxl::read_excel(path = file, sheet = sheet, col_types = "text", ...)
  } else if (extension %in% c("csv", "tsv", "tab", "txt") || !is.null(sep))
  {
    if (is.null(sep)) { sep <- if (extension == "csv") "," else "\t" }
    # Everything is read as text (the cells contain lists of tips, or names of groups), blank cells are missing values, and the names of the columns are kept as they are
    data <- utils::read.table(file = file, header = TRUE, sep = sep, quote = "\"", colClasses = "character",
                              na.strings = c("", "NA"), strip.white = TRUE, check.names = FALSE, comment.char = "",
                              fill = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8", ...)
  } else {
    stop("Unknown file extension '.", extension, "'. Use a '.csv', '.tsv', '.tab', '.txt', '.xlsx' or '.xls' file, ",
         "or provide 'sep' (for a text file) or 'reader' (a function that reads the file).", call. = FALSE)
  }

  if (!is.data.frame(data)) { stop("The table read from '", file, "' is not a data frame", call. = FALSE) }
  if (nrow(data) == 0) { stop("The table read from '", file, "' has no rows", call. = FALSE) }

  return(as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE))
}
