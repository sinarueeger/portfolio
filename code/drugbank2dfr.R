###########################################################################
###########################################################################
###                                                                     ###
###                   DRUGBANK XML to DATAFRAME                         ###
###                                                                     ###
###########################################################################
###########################################################################


# packages needed ---------------------------------------------------------

## loaded directly via ::
## library(here)
## library(data.table)
## library(parallel)
## library(import)
## library(dplyr)

library(XML)
import::from(magrittr, "%>%")


# get data ----------------------------------------------------------------

## download from: https://www.drugbank.ca/releases/latest
drugbank.xml <- "drugbank_20190210.xml"


# settings ----------------------------------------------------------------

PATH_IN <-
  here::here("data", "raw", "drugbank.xml") ## writing protected
PATH_OUT <-
  here::here("data", "processed", "drugbank_data_filtered.RData")
NUMBER_OF_CORES <-
  4 ## numbers of cores to use (within the parallel package)




##////////////////////////////////////////////////////////////////
##                         XML DATA                             //
##////////////////////////////////////////////////////////////////

## Parse XML file into a dataframe that contains each drug as a line,
## and features as columns

## Code snippets from here: http://www.informit.com/articles/article.aspx?p=2215520

# load xml file -----------------------------------------------------------
xmlfile <- xmlParse(PATH_IN)

# give content of root ----------------------------------------------------
xmltop <- xmlRoot(xmlfile)


##////////////////////////////////////////////////////////////////
##               FUNCTION TO EXTRACT XML DATA                   //
##////////////////////////////////////////////////////////////////

#' Function to extract drugbank xml data by features (for a single drug!)
#'
#' @param feature character describing the feature to be extracted, e.g. "drugbank-id"
#' @param data XML entry for a single drug processed with xmlRoot(xmlParse(XML-file)).
#'
#' @return scalar
#'
#' @examples
#' ## Extract drugbank-id and name for third drug in xml file
#' data <- xmltop[[3]]
#' extract.xml("drugbank-id", data)
#' extract.xml("name", data)

extract.xml <- function(feature, data)
{
  tmp <- data[[feature]]
  
  ## if empty, return NA
  if (xmlSize(tmp) == 0)
  {
    return(NA)
  }
  
  ## if 1, return character
  if (xmlSize(tmp) == 1 &
      !(feature %in% c("transporters", "targets", "enzymes", "carriers")))
  {
    return(as(tmp[[1]], "character"))
  } else{
    ## if > 1
    if (feature == "categories")
    {
      tmp2 <-
        sapply(1:xmlSize(tmp), function(x)
          as(tmp[[x]][["category"]][[1]], "character"))
      
      tmp.out <- paste0(tmp2, collapse = "/")
      return(tmp.out)
      
    } else{
      if (feature %in% c("transporters", "targets", "enzymes", "carriers"))
      {
        tmp2 <- sapply(1:xmlSize(tmp), function(x)
        {
          gene <- tmp[[x]][[7]][["gene-name"]][[1]]
          if (!is.null(gene))
          {
            as(gene, "character")
          } else{
            NA
          }
        })
        
        ## paste
        tmp.out <- paste0(na.omit(tmp2), collapse = "/")
        return(tmp.out)
        
        ## if larger and not present, then NA
      } else{
        return(NA)
      }
    }
  }
  
}

#' Function to assemble extracted XML entries
#' (wrapper of extract.xml, applying it to a vector of features)
#'
#' @param features character describing the features to be extracted, e.g. "drugbank-id" or "target"
#' @param data XML entry for a single drug processed with xmlRoot(xmlParse(XML-file)).
#'
#' @return vector
#'
#' @examples
#' ## Extracting drugbank-id, name, indication, description and gene targets for drug number 3.
#' data <- xmltop[[3]]
#' features <- c("drugbank-id", "name", "indication", "description", "targets")
#' xml2df(features, data)
xml2df <- function(features, data)
{
  out <-
    data.frame(t(unlist(
      sapply(features, function(x)
        extract.xml(x, data))
    )))
  
  return(out)
  
}


##////////////////////////////////////////////////////////////////
##                      EXTRACT XML DATA                        //
##////////////////////////////////////////////////////////////////

# loop through each drug  ----------------------------------------------

features <-
  c(
    "drugbank-id",
    "name",
    "indication",
    "description",
    "snp-effects",
    "snp-adverse-drug-reactions" ,
    "categories",
    "targets",
    "transporters",
    "enzymes",
    "carriers"
  )

db_list <- parallel::mclapply(1:xmlSize(xmltop), function(k)
{
  ## first, we apply xml2df for each drug k
  ## xml2df will itself loop through each feature
  out <-  xml2df(features, xmltop[[k]])
  
  return(out)
  
}, mc.cores = NUMBER_OF_CORES)


# turn list into dataframe ------------------------------------------------

db <- data.table::rbindlist(db_list)


# turn factors into characters --------------------------------------------

db <- (db %>% dplyr::mutate_if(is.factor,
                               as.character))

# save the data frame -----------------------------------------------------

save(db, file = PATH_OUT)
