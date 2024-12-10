#sub-divide IDPs into the largest possible complete chunks for ComBat harmonization

#tempdir()
#get_hostname <- function(){
#	    return(as.character(Sys.info()["nodename"]))
#}
#print(get_hostname())
 
#system('du -s /tmp')
#system('df -h /tmp')

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(rlang)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df_orig <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
idp_list <- readRDS(args[2]) #path to .rds obj of list
save_path <- as.character(args[3])
name_prefix <- as.character(args[4])

#DEFINE FUNCTIONS

#function to split out complete dataframes
get_complete_df <- function(var, df, batch=NULL){
  #print(var)
  new_df <- df %>%
    dplyr::filter(!is.na(!!sym(var))) %>% #can't handle var name as string, fix this!
    select_if(~ !any(is.na(.)))
  if (!is.null(batch)){
    new_df <- new_df %>%
      group_by(!!sym(batch)) %>%
      filter(n() >=5) %>% #remove sites with < 5 ppl
      ungroup()
  }
  return(new_df)
}

#re-combine dfs for any phenotypes that are identical - code with help from ChatGPT
merge_identical_dfs <- function(df_list) {
  df_names <- names(df_list)
  processed <- rep(FALSE, length(df_list)) # Track which dataframes have been processed
  new_list <- list() # Initialize an empty list to store merged dataframes
  
  for (i in seq_along(df_list)) {
    print(paste("testing dataframe #", i))
    if (!processed[i]) {
      identical_dfs <- c(df_names[i]) # Start with the current dataframe's name
      
      # Compare the current dataframe with others
      for (j in (i + 1):length(df_list)) {
        if (j > length(df_list)) next  # Avoid index out of bounds
        print(paste("compare to dataframe #", j))
        if (!processed[j] && identical(df_list[[i]], df_list[[j]])) {
          identical_dfs <- c(identical_dfs, df_names[j])
          processed[j] <- TRUE # Mark the duplicate dataframe as processed
        }
      }
      
      # If there are any identical dataframes, rename the first one
      if (length(identical_dfs) > 1) {
        new_name <- paste(identical_dfs, collapse = "_") # Concatenate the names
        new_list[[new_name]] <- df_list[[i]] # Add the merged dataframe with the new name
      } else {
        # If no identicals were found, keep the original dataframe
        new_list[[df_names[i]]] <- df_list[[i]]
      }
      
      processed[i] <- TRUE # Mark the current dataframe as processed
    }
  }
  
  return(new_list)
}

rename_cols <- function(df_list, col_list){
  lapply(names(df_list), function(df_name) {
    # Construct the file path using the dataframe's name
    combat_list <- strsplit(gsub("_", " ", df_name), " ")[[1]]
    
    #find any phenotypes that are NOT the combat targets
    not_target <- setdiff(col_list, combat_list)
    #print(not_target)
    
    #see which ones are in the dataframe
    not_target <- intersect(not_target, names(df_list[[df_name]]))
    #print(not_target)
    
    if (length(not_target)>0){
      df_list[[df_name]] <- df_list[[df_name]] %>%
        rename_with(~ paste(., "X", sep = "_"), any_of(not_target))
    }
    return(df_list[[df_name]])
  })
}

#write out csvs - written with help from ChatGPT
save_csv_list <- function(df_list, path = ".", name_prefix) {
  lapply(seq(length(df_list)), function(idx) {
    # Construct the file path using the dataframe's name
    file_path <- paste0(path,"/", name_prefix, "_", idx, ".csv")
    
    # Use fwrite to save the dataframe
    fwrite(df_list[[idx]], file_path)
  })
}

#PROCESS DATA

#split out data
df_list <- lapply(idp_list, get_complete_df, df_orig, batch="study_site")
#name
names(df_list) <- idp_list

#recombine & rename any identical dfs
df_list_clean <- merge_identical_dfs(df_list)

#identify and mark non-target IDPs (those being used as priors) with _X
df_list_final <- rename_cols(df_list_clean, idp_list)

#save
save_csv_list(df_list_final, save_path, name_prefix)

print(paste("Saved", length(df_list_final), "files :)"))
