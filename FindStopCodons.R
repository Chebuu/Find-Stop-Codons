## TODO:
## Support for PAMs not included. Only parameters CDS, stops, inc.codon, inc.position start are supported right now

# BiocManager::install("Biostrings", version = "3.8")
library(Biostrings)

UtoT <- function(string){
  # Convert instances of U to T in a given string or character vector
  if(length(string)==1){
    string <- strsplit(toupper(string), '')[[1]]  
  }
  string[string == 'U'] <- 'T'
  return(string)
}
# UtoT(c('A','U','G')) == UtoT('A''U''G')

permuteSTOPS <- function(stops){
  # Returns a matrix of all single nucleotide polymorphisms of a given codon in a list of codons
  # @ param stops A vector of codons as strings
  bases <- c('a','g','c','t')
  stops <- toupper(stops)
  perm.out <- sapply(stops, function(stop){
    stop <- strsplit(stop, '')[[1]]
    stop.perms <- c()
    for(s in 1:length(stop)){
      for(base in bases){
        stop.new <- stop
        stop.new[s] <- base
        stop.new <- paste(stop.new, collapse='')
        stop.perms <- c(stop.perms, stop.new)
      }
    }
    return(stop.perms)
  })
  rownames(perm.out) <- rep(bases, times=length(stops))
  perm.out <- toupper(perm.out)
  return(perm.out)
}
# permuteSTOPS(c('UAG', 'UgA', 'UAA'))

findSTOPS <- function(CDS, stops=c('UAG', 'UGA', 'UAA'), start=1, PAM=c(), PAM.offset=c(), inc.position=T, inc.codon=T, PAMonly=F){
  # Find codons in the CDS that can be mutated at a single base to yield a stop codon. If the codon is already a stop codon, the program returns 0 as the mutation index.
  # Yields a list where each element represents a codon and the format of each element is canonical mutation format e.g. c("AAA", "A", "1537", "T", "UAA") == c(<original_codon>, <original_base>, <mutation_position>, <resulting_base>, <resulting_codon>)
  # @ param CDS The coding sequence as a DNAStringSet. CDS is a misnomer because the DNSStringSet can include UTRs.
  # @ param stops A vector of stop codons
  # @ param start The +1 position of the coding sequence. The index of A in the AUG start codon.
  # @ param PAM A list of PAM sequences as strings
  # @ param inc.position Whether to include the position of the codon in the output i.e. the position of the amino acid.
  # @ param inc.codon Whether to include the position of each mutation in the output as well as the resulting stop codon. 
  # @ param PAMonly Whether to include only codons that fall within the given PAM.offset from each PAM
  
  if(class(CDS) == 'DNAStringSet' | class(CDS) == 'list'){
    return(findSTOPS.multiCDS( as.list(CDS) ))
  }
  CDS <- strsplit(as.character(CDS), '')[[1]]
  CDS <- CDS[start:length(CDS)]
  stops <- sapply(stops, function(stop) {
    return(paste(UtoT(stop), collapse=''))
  })
  stop.perms <- permuteSTOPS(stops)
  stops.possible <- vector(mode='list', length=ceiling(length(CDS)/3))
  codon.current <- 0
  for(i in seq(3, length(CDS), by=3)){
    codon.current <- codon.current + 1 # the current codon (amino acid position)
    position.current <- codon.current * 3 - 2 # current position of codon
    codon <- paste(CDS[(i-2):i], collapse='') # the actual codon itself
    perms.which <- codon == stop.perms # which stop codons could this codon become with a single mutation
    if(!any(perms.which)){
      stops.possible[[codon.current]] <- FALSE # there is no mutation that will create a stop codon
    }else{
      res <- unique(stop.perms[perms.which]) # the original codon
      app <- unique(colnames(stop.perms)[which(perms.which, arr.ind = T)[,2]]) 
      if(inc.position & !inc.codon){
        res <- c(res, position.current) # the original codon and the position
      }
      if(inc.codon){
        app.position <- sapply(app, function(a){
          a <- UtoT(a)
          codon.split <- strsplit(codon,'')[[1]]
          if(all(a == codon.split)){
            return(0) # the current codon is a stop codon, return 0
          }else{
            mut.position <- which(a != codon.split) # the position of the mutation in the codon, values are 1:3
            app.position <- position.current  + mut.position - 1 # the position of the mutation in the CDS
            mut <- a[mut.position] # the mutated base
            org <- codon.split[mut.position] # the orignal base
            res <- unlist(c(org, as.vector(rbind(app.position, mut)))) # HACKY collapse workaround yields e.g. c("AAA", "A", "1537", "T", "UAA") == c(<original_codon>, <original_base>, <mutation_position>, <resulting_base>, <resulting_codon>)
            return(res)  
          }
        })
        res <- unlist(c(res, as.vector(rbind(app.position, app)))) 
      }
      stops.possible[[codon.current]] <- res
    }
  } 
  return(stops.possible)
}
# findSTOPS(readDNAStringSet('/path/to/dna.fasta'))

findSTOPS.multiCDS <- function(CDS.list, opts.list=NULL, ...){
  # Apply the findSTOPS funtion to a list of DNA sequences
  # @ param CDS.list A list of DNAStringSets or a list of strings
  # @ param opts.list A list of arguments to be passed with each item in CDS.list. Indices should match idices in CDS.list. NULL can be used for default arguments or to use ... arguments instead of list arguments.
  # @ param ... Arguments to be passed to findSTOPS
  if(class(CDS.list)=='DNAStringSet'){
    CDS.list <- as.list(CDS.list)
  }
  res <- lapply(1:length(CDS.list), function(CDS.idx){
    opts <- opts.list[[CDS.idx]]
    CDS <- CDS.list[[CDS.idx]]
    if(!is.null(opts)){
      return(findSTOPS(CDS, opts))
    }else{
      return(findSTOPS(CDS, ...))
    }
  })
  names(res) <- names(CDS.list)
  return(res)
}
# findSTOPS(readDNAStringSet('/path/to/dnas.fasta'))


## TESTS
# DATA_DIR <- './data/'
# fastafiles <- paste0(DATA_DIR, dir(DATA_DIR))
# 
# CDS.test <- readDNAStringSet( './data/KT876022.fasta' )
# 
# findSTOPS(CDS.test)
# length(CDS.test) == length(findSTOPS(CDS.test))
# 
# findSTOPS(CDS.test[[1]])
# length(CDS.test[[1]]) / 3 == length(findSTOPS(CDS.test[[1]]))

