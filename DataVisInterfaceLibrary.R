suppressMessages(library(cli))
suppressMessages(library(argparse))

## DATA MENU FUNCTIONS ##

# Run the whole interactive menu
# Only function called outside library
interactive_menu <- function(){

    ## Display initial message
    header <- "Plant breeding simulations Data Visualizer"
    cli_h1(header)
    cli_text()
    cli_h3("Welcome!")
    cli_text("This programme is for visualizing the data generated from our simulation.")
    reminder <- "To make sure the visualization successful, please do not modify the names of the data files."
    cli_text(cli_alert_info(reminder))
    cli_text()

    folder_prompt = "Please enter the name of the folder you want to visualize or 'q' to quit the programme: "

    set <- FALSE
    while (set == FALSE) {
        # Get folder name from the command line
        folder_name <- read_input(folder_prompt)

        set <- read_input_folder(folder_name)
        
        # check if user choose to terminate the programme
        if (folder_name == "q") {
            ans <- read_input_exit("Do you want to quit the programme [Y/n]? ")
            if (ans == TRUE) {
                cli_text("Exiting the program...")
                Sys.sleep(1)
                q(save = "no")
            } else if (ans == FALSE) {
                set <- FALSE
            } else {
                read_input_exit("Do you want to quit the programme [Y/n]? ")
            }
        } 
    }

    
    # Option list for choosing what graph(s) is/are going to generate
    option <- list(c("Correlations (cors)", "Genetic Values (gvs)", "Variances (vars)", "All"))

    chosen <- FALSE
    # while (chosen != TRUE) {
    #     print_options(option)
    #     choice <- read_choice(length(option[[1]]))
    # }
}

# Print options to choose which data type is going to be visualized
print_options <- function(option) {
    cli_ol()
    for (i in seq_along(option[[1]])) {
        cli_li(option[[1]][i])
    }
    cli_end()
    cli_text()
}


## INPUT FUNCTIONS ##


# Prompt a text and read input
# Returns input value.
read_input <- function(prompt=""){
    cat(prompt)
    type.convert(readLines("stdin", 1), as.is=TRUE)
}


# Function to check if a folder exists
# Return boolean for folder existence or q for quitting
read_input_folder <- function(folder_name) {
    if (file.exists(folder_name) && file.info(folder_name)$isdir) {
        return(TRUE)
    } else if (length(folder_name) == 0) {
        cli_alert_warning("Please enter a valid folder name")
        cli_text()
        return(FALSE)
    } else if (folder_name == "q") {
        return("q")
    } else {
        cli_alert_warning("The folder does not exist. Please try again.")
        cli_text()
        return(FALSE)
    }
}

# Read input for exit option (Y/n)
# return boolean 
read_input_exit <- function(prompt="") {
    exit <- read_input(prompt)
    if (exit == "Y") {
        return(TRUE)
    } else if (exit == "n") {
        return(FALSE)
    } else {
        cli_alert_warning("Not a valid option")
        return(read_input_exit(prompt))
    }
}

# Read integer for user's choice of data type to do visualization on
# Return integer
read_choice <- function(max) {
    prompt <- "Please choose which data you would like to visualize by typing the {.emph number}"
    choice <- read_input(prompt)


}



interactive_menu()