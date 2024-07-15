#!/usr/bin/env Rscript

# A script to setup dependencies using renv. It can be used to install
# dependencies when the renv.lock file is missing or out-of-date. Inspired by
# the python poetry tool chain.
#
# By default if there are discrepencies between the human defined dependencies
# and the renv.lock then an instructional error message witll be printed. But
# this can be bypassed with --force.
#
# Run with --help for usage information.
#
# Currently supports: dependencies.R file format. Does not support: DESCRIPTION
# file format.

#############
# CONSTANTS #
#############

SCRIPT_VERSION <- "0.1.0"
MISSING_VERSION <- "MISSING_VERSION"


#############
# FUNCTIONS #
#############


# Helper function to parse a package string into name and version
parse_renv_package_string <- function(package_string) {
    # Remove leading/trailing whitespace
    package_string <- trimws(package_string)
    
    # Remove quotes if present
    package_string <- gsub('"', '', package_string)
    
    # Split the string by "::" to separate source and package
    # e.g., 'bioc::GenomicRanges@1.50.2' -> 'bioc' and 'GenomicRanges@1.50.2'
    raw_package_string <- tail(strsplit(package_string, "::")[[1]], 1)

    # Split the raw package string by "@" to separate name and version
    package_parts <- strsplit(raw_package_string, "@")[[1]]
    
    # Extract the package name
    name <- package_parts[1]
    
    # Extract the version if available
    if (length(package_parts) > 1) {
        version <- package_parts[2]
    } else {
        version <- MISSING_VERSION
    }
    
    # If the name resembles a path extract the last part as the name
    # e.g. 'VanLoo-lab/ascat/ASCAT' -> 'ASCAT'
    if (grepl("/", name)) {
        name <- sub(".*/(.*)", "\\1", name)
    }
    
    # Return a list with name and version
    return(list(name = name, version = version))
}


# Helper function to get package version from lockfile
get_package_version_from_lock <- function(lockfile, package_name) {
    # Set a flag to show membership of the lockfile
    has_package <- package_name %in% names(lockfile$Packages)

    # If the package is in the lockfile, get the version, defaulting to
    # MISSING_VERSION if its absent. Otherwise, return MISSING_VERSION
    if (has_package) {
        package_info <- lockfile$Packages[[package_name]]
        if (is.null(package_info$Version)) {
            result <- MISSING_VERSION
        } else {
            result <- (package_info$Version)
        }
    } else {
        result <- MISSING_VERSION
    }

    return(result)
}


# Function to parse and compare dependencies
check_configuration_discrepancies <- function(lockfile, required_packages) {
    # Parse the required packages strings into a more convenient structure
    # Going from ["dplyr@1.1.2", "bioc::rtracklayer@1.58.0"]
    # to [[name = "dplyr", version = "1.1.2"], [name = "rtracklayer", version = "1.58.0"]]
    parsed_required_packages <- lapply(required_packages, parse_renv_package_string)
    
    # Check for discrepancies
    discrepancies <- FALSE
    for (package_info in parsed_required_packages) {
        package_name <- package_info$name
        required_version <- package_info$version
        
        # Get version from lockfile using the extracted package name. This can return a semver or MISSING_VERSION
        lock_version <- get_package_version_from_lock(lockfile, package_name)

        is_both_matching <- lock_version == required_version
        is_missing_in_lockfile <- lock_version == MISSING_VERSION

        # Check if package is missing in lockfile (i.e. not up-to-date lock) --> discrepancy!
        if (is_missing_in_lockfile) {
            discrepancies <- TRUE
            break
        }

        # Check if package version match --> not a discrepancy
        if (is_both_matching) {
            next
        }

        # Everything else is a discrepancy (e.g. different versions and so not up-to-date lock)
        discrepancies <- TRUE
        break
    }

    return(discrepancies)
}


# Function to read dependencies from a file, depending on the file format: dependencies.R or DESCRIPTION
get_required_packages <- function(dependencies_file) {
    # Check the file exists otherwise, stop with a warning
    if (!file.exists(dependencies_file)) {
        msg = paste("Dependencies file not found. Please ensure '", dependencies_file, "' is present in the working directory.")
        stop(msg)
    }

    # If the file is a dependencies.R file, source it and return the
    # REQUIRED_PACKAGES value, if its a DESCRIPTION file, raise not implemented
    # error and anythign else raise an error
    if (grepl("dependencies.R", dependencies_file)) {
        source(dependencies_file)
        return(REQUIRED_PACKAGES)
    } else if (grepl("DESCRIPTION", dependencies_file)) {
        stop("Parsing of DESCRIPTION file format not implemented yet.")
    } else {
        msg = paste("Unsupported file format. Got: '", dependencies_file, "'. Supported formats are 'dependencies.R' and 'DESCRIPTION'.")
        stop(msg)
    }
}


handle_discrepancies <- function(required_packages, force) {
    # A multiline message to inform the user about the discrepancies
    msg <- paste("Configuration discrepancies found. You have the following options:\n",
                "1. Install the dependencies by running 'source(\"dependencies.R\"); renv::install(REQUIRED_PACKAGES)' and then 'renv::snapshot()'.\n",
                "2. Re-run this script with the '--force' flag to force installation of dependencies (same as option 1).\n",
                "3. Run 'renv::restore()' to restore the environment from the lockfile.\n",
                "4. Delete 'renv.lock' and run this script again to recreate the lockfile.",
                collapse = "")
    if (force) {
        message("Forcing installation of dependencies.")
        renv::install(required_packages)
        renv::snapshot()
    } else {
        message(msg)
        quit(status = 1)
    }
}


# Function to handle dependency setup
setup_dependencies <- function(lockfile_path, required_packages, force) {

  if (!file.exists(lockfile_path)) {
    message("No lockfile found. Installing dependencies and creating a lockfile.")
    renv::install(required_packages)
    renv::snapshot()
  } else {
    # Parse the lockfile into its native data structure
    lockfile <- renv::lockfile_read(file = lockfile_path)

    # Check for configuration discrepancies, i.e. is the lockfile up-to-date
    # (missing required packages or they have different versions)
    discrepancies_found <- check_configuration_discrepancies(lockfile, required_packages)
    
    if (discrepancies_found) {
        message("There were R package discrepancies found between human defined dependencies and the lockfile")
        # Either exit with user instructions or if forced, install the
        # dependencies and update the lockfile
        handle_discrepancies(required_packages, force)
    } else {
        message("No R package discrepancies found. Restoring environment from lockfile.")
        renv::restore()
    }
  }
}


simple_argparse <- function() {
    # Optparse would be better but we want to avoid a dependency for this
    # bootstrapping script.

    # Get the full list of command-line arguments
    args <- commandArgs(trailingOnly = TRUE)

    # Initialize variables for flags
    force <- FALSE

    # Attempt to retrieve the file name from the call stack
    script_arg <- commandArgs(trailingOnly = FALSE)[4]
    script_name <- basename(sub("^--file=", "", script_arg))

    # Check for and process arguments
    if (length(args) > 0) {
        # Loop through each argument
        for (i in seq_along(args)) {
            switch(args[i],
                "--help" = {
                    cat(paste("Usage: Rscript", script_name, "[--force] [--help]", "\n"))
                    cat("Options:\n")
                    cat("--help  Show this help message\n")
                    cat("--force  Force the installation of dependencies and update the lockfile\n")
                    quit(status = 0)  # Exit after showing help
                },
                "--version" = {
                    cat(paste("Version:", SCRIPT_VERSION, "\n"))
                    quit(status = 0)  # Exit after showing version
                },
                "--force" = {
                    force <- TRUE  # Set force to TRUE if --force is present
                },
                {
                    cat("Unknown argument:", args[i], "\n")
                    cat("Use --help for usage information.\n")
                    quit(status = 1)  # Exit with error
                }
            )
        }
    }

    # Now you can use the 'force' variable in your script logic
    if (force) {
        warning("Force mode enabled.")
    }
    return(list(force = force))
}


main <- function(force = FALSE) {
    # Parse the CLI args
    cli_args <- simple_argparse()
    force <- cli_args$force

    if (!requireNamespace("renv", quietly = TRUE)) {
        message("renv is not installed. Please install renv and try again.")
        quit(status = 1)
    }

    # Define paths to the lockfile and dependencies file
    lockfile_path <- "renv.lock"
    dependencies_file <- "dependencies.R"

    required_packages <- get_required_packages(dependencies_file)

    # Run the setup dependencies function
    setup_dependencies(lockfile_path, required_packages, force)

    print("Dependencies setup complete.")
}

###########
# RUNTIME #
###########


if (interactive()) {
    msg <- paste("This script is not meant to be run interactively. Please use Rscript to run it from the command line.")
    stop(msg)
} else {
    main()
}