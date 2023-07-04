suppressMessages(library(cli))
suppressMessages(library(argparse))

## MENU FUNCTIONS ##

# Run the whole interactive menu
# Only function called outside library
interactive_menu <- function(){

    ## Display initial message
    header <- "Plant breeding simulations"
    cli_h1(header)
    cli_text()
    cli_text("Welcome to plant breeding simulations!")
    cli_text("Before starting, choose the simulation parameters.")
    cli_text()
    inp <- read_input("Press enter to continue...")

    done <- FALSE
    while(done == FALSE){
        ## Prompts parameters to user and allows edition
        params_menu(sim_params, "Simulation")
        params_menu(run_params, "Runtime")

        ## Displays final settings and asks for confirmation
        inp <- confirm_settings(header)

        done <- (inp == 1)
    }

    cli_text("Proceeding to simulation...")
}

# Displays editable parameters and allows user to set them.
params_menu <- function(params_list, setting_name="Other"){
    header <- "{setting_name} settings"
    msg <- "" ## Feedback message when parameter is changed

    while (TRUE){
        # Print header
        cli_h1(header)
        
        ## If there is a message, displays and resets it
        if (msg != "") { 
            cli_alert_success(msg)
            msg <- ""
        }

        ## Prints initial text and available parameters
        cli_text()
        cli_text(paste("Choose an option to edit or", (col_br_cyan("0 to confirm and proceed"))))
        cli_text()
        print_params(params_list, args)
        
        # Get user input
        inp <- read_input_range(lower=0, higher=length(params_list))
        
        if (inp == 0) # Finish current setting
            break
        
        param_name <- names(params_list[inp])
        param <- params_list[[inp]]

        # Print header and text for selected parameter
        cli_h1(header)
        cli_text()
        cli_text("Select {param$help}")
        
        # If parameter has list of choices, display them to user
        # If not, user inputs answer manually
        if (is.null(param$choices) == FALSE){
            cli_text()
            print_list(param$choices, args[[param_name]])
            ind <- read_input_range(length(param$choices))
            
            args[[param_name]] <<- param$choices[ind]
        } else {
            cli_alert_info("Current value: { col_br_cyan(args[[param_name]]) }")
            cli_text()
            args[[param_name]] <<- read_input_type( typeof(args[[param_name]]) )
        }

        # Create feedback message to next menu
        msg <- "{ col_br_cyan(param$help) } has been set to { col_br_cyan(args[[param_name]]) }"
    }
}

# Print list of parameters and their current values between parentheses
print_params <- function(options, args){
    cli_ol()
    for (opt_name in names(options)){
        opt <- options[[opt_name]]
        cli_li("{opt$help} ({ col_br_cyan(args[[opt_name]]) })")
    }
    cli_end()
    cli_text()
}

# Print list of options and highlight current set value
print_list <- function(list, curVal=NULL){
    cli_ol()
    for (item in list){
        if (!is.null(curVal) && item == curVal)
            cli_li("{col_br_cyan(item)}")
        else
            cli_li(item)
    }
        
    cli_end()
    cli_text()
}

# Displays final settings and asks for confirmation
confirm_settings <- function(header){
    cli_h1(header)
    cli_alert_info("Here are your settings.")

    cli_h3("Simulation settings")
    print_params(sim_params, args)
    cli_h3("Runtime settings")
    print_params(run_params, args)
    
    cli_h3("Confirm settings and start simulation?")
    print_list(c("Yes", "No"))
    inp <- read_input_range(2)

    inp
}


## INPUT FUNCTIONS ##

# Prompt a text and read input
# Returns input value.
read_input <- function(prompt=""){
    cat(prompt)
    type.convert(readLines("stdin", 1), as.is=TRUE)
}

# Read an input of a specific type, sending an error message if it's wrong.
# Returns input value.
read_input_type <- function(type){
    inp <- read_input("New value: ")
    while(typeof(inp) != type){
        cli_alert_warning("Wrong type. Variable type should be {col_br_cyan(type)}")
        inp <- read_input("New value: ")
    }
    inp
}

# Read integer number in a specific range (lower <= num <= higher).
# Returns input value.
read_input_range <- function(higher, lower=1){
    prompt <- paste("Choose a value between ", col_br_cyan(lower), " and ", col_br_cyan(higher), ": ", sep="")
    inp <- read_input(prompt)

    while(typeof(inp) != "integer" || inp < lower || inp > higher){
        cli_alert_warning("Invalid choice!")
        inp <- read_input(prompt)
    }
    inp
}
