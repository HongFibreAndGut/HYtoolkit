
#' @title ClearGenusSilva138
#' @description Make clear and short genus name(s) for taxon based on SILVA 138 (i.e., k__;p__;c__;o__;f__;g__)
#' @param raw_taxon A character vector of taxa, in form of k__xxxx;p__xxxx;c__xxxx;o__xxxx;f__xxxx;g__xxxx
#'
#' @return Returns a character vector of the genus taxa
#' @export
#'
#' @examples raw_taxon <- c("Unassigned;__;__;__;__;__",
#' "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter",
#' "k__Bacteria;p__Actinobacteriota;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriales_Incertae_Sedis;g__uncultured")
#' new_taxon <- ClearGenusSilva138(raw_taxon)

ClearGenusSilva138 <- function(raw_taxon){
  cleargenus <-
    ifelse(grepl(";g__UCG", raw_taxon) & grepl(";f__UCG", raw_taxon) == F, # Add the family information after the UCG
           paste(gsub(".*;f__(.+);g__.*","\\1",raw_taxon), gsub(".*;g__","\\1",raw_taxon)),
           ifelse(grepl(";g__UCG", raw_taxon) & grepl(";f__UCG", raw_taxon),
                  paste(gsub(".*;o__(.+);f__.*","\\1",raw_taxon), gsub(".*;f__(.+);g__.*","\\1",raw_taxon)),
                  ifelse(grepl(";g__uncultured", raw_taxon) &
                           (grepl(";f__uncultured", raw_taxon) == F) &
                           (grepl(";o__uncultured", raw_taxon) == F), # Add the family information after the uncultured
                         paste(gsub(".*;g__","\\1",raw_taxon), gsub(".*;f__(.+);g__.*","\\1",raw_taxon)),

                         ifelse(grepl(";g__uncultured", raw_taxon) &
                                  (grepl(";f__uncultured", raw_taxon) == T) &
                                  (grepl(";o__uncultured", raw_taxon) == F),
                                paste(gsub(".*;f__(.+);g__.*","\\1",raw_taxon),gsub(".*;o__(.+);f__.*","\\1",raw_taxon)),

                                ifelse(grepl(";g__uncultured", raw_taxon) &
                                         (grepl(";f__uncultured", raw_taxon) == T) &
                                         (grepl(";o__uncultured", raw_taxon) == T),
                                       paste(gsub(".*;o__(.+);f__.*","\\1",raw_taxon),gsub(".*;c__(.+);o__.*","\\1",raw_taxon)),

                                       ifelse(grepl(";g__Family", raw_taxon),
                                              paste(gsub(".*;f__(.+);g__.*","\\1",raw_taxon),gsub(".*;g__(.+)","\\1",raw_taxon)),
                                              ifelse(grepl(";g__", raw_taxon),  # Add the f/o/c/p/d information after the ;__
                                                     gsub(".*;g__","\\1",raw_taxon),

                                                     ifelse(grepl(";f__", raw_taxon),
                                                            paste("unclassified",gsub(".*;f__(.+?);.*","\\1",raw_taxon)),
                                                            ifelse(grepl(";o__", raw_taxon),
                                                                   paste("unclassified",gsub(".*;o__(.+?);.*","\\1",raw_taxon)),
                                                                   ifelse(grepl(";c__", raw_taxon),
                                                                          paste("unclassified",gsub(".*;c__(.+?);.*","\\1",raw_taxon)),
                                                                          ifelse(grepl(";p__", raw_taxon),
                                                                                 paste("unclassified",gsub(".*;p__(.+?);.*","\\1",raw_taxon)),
                                                                                 ifelse(grepl("k__", raw_taxon),
                                                                                        paste("unclassified",gsub(".*k__(.+?);.*","\\1",raw_taxon)),
                                                                                        "unassgined taxon"
                                                                                 ))))))))))))
  return(cleargenus)
}



