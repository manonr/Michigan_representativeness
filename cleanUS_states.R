extract_us_state <- function(x) {
  # List of state names and abbreviations
  state_names <- tolower(state.name)
  state_abbs <- tolower(state.abb)
  names(state_names) <- state_names
  names(state_abbs) <- state_names  # to map abbreviation back to full name
  
  # Combine into one named vector for matching
  state_lookup <- c(state_names, state_abbs)
  
  # Function to extract one state from a string
  extract_one <- function(string) {
    string_lower <- tolower(string)
    
    # Match full state name
    for (state in names(state_names)) {
      if (grepl(paste0("\\b", state, "\\b"), string_lower)) return(state.name[which(state_names == state)])
    }
    
    # Match abbreviation
    for (abbr in names(state_abbs)) {
      abbr_code <- state_abbs[abbr]
      if (grepl(paste0("\\b", abbr_code, "\\b"), string_lower)) return(state.name[which(state_abbs == abbr_code)])
    }
    
    return(NA)
  }
  
  # Apply over vector
  sapply(x, extract_one, USE.NAMES = FALSE)
}


#df <- data.frame(location = c(
#  "Lansing, MI", 
#  "California coast", 
#  "random text", 
#  "Austin TX", 
#  "north carolina",
#  "Seattle, WA",
#  "Boston,Massachusetts"
#))



#speicifc issue with Baltimore, make those "MAryland



#setwd("/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/raw-data")
#df <-  read.csv("H1_GISAID_GENBANK_unique.csv")
#df$state <- extract_us_state(df$Location)
