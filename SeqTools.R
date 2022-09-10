library(tidyverse)

gbfile <- '/Users/max/Downloads/GCF_013167035.1_ASM1316703v1_genomic.gbff'
fixedfile <- gsub('([.][^.]*$)', '_fixed\\1', gbfile)

# first gsed marks lines that start with /, next one replaces the preceding \n, next two get rid of < and > for truncated annotations last one converts join(A..B,C..D) to A..D, if B == C
system(paste0('cat ', gbfile, ' | gsed -E \'s/^        *[^ /]/;/\' | gsed -z \'s/\\n;//g\' | gsed -E \'s/[(]</(/\' | gsed \'s/..>/../\' | gsed -E \'s/join[(]([0-9]*)[.][.]([0-9]*),\\2[.][.]([0-9]*)[)]/\\1..\\3/\'> ', fixedfile))

raw <- readLines(fixedfile)

gb_sections <- which(grepl('(^[A-Z]|//)', raw))
start_features <- which(grepl('^FEATURES', raw)) + 1
features <- raw[start_features:(gb_sections[which.max(gb_sections > start_features)] -1)]
features <- grep('^ +[A-z]', features, value = TRUE) %>% gsub('  *', '',.)
features <- data.frame(name = gsub('[0-9]*[.][.].*', '', features), 
                       position = str_extract(features, '[0-9]*[.][.][0-9]*'))

features$start <- str_split(features$position, '[.][.]') %>% sapply("[[", 1)
features$end <- str_split(features$position, '[.][.]') %>% sapply("[[", 2)
features$position <- 'plus_strand'
features[str_detect(features$name, 'complement[(]$'), 'position'] <- 'minus_strand'
features$name <- str_replace(features$name, 'complement[(]$', '')
features <- features[!str_detect(features$name, '[(]$'),]

start_origin <- which(grepl('^ORIGIN', raw))
seq <- raw[start_origin:gb_sections[which.max(gb_sections > start_origin)]]
seq <- paste(gsub('[^a-z]', '', seq), collapse = '')


attrs <- grep('^  */', raw, value =TRUE) %>% str_extract('(?<=/)[^=]*(?==)') %>% unique()
attrs <- c('sequence', attrs[!is.na(attrs)])
features[, attrs] <- NA

features[2:9117, 'sequence'] <- apply(features[2:9117,], 1, function(x) substr(seq, x['start'], x['end']))
