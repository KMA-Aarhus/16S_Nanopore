library(tidyverse)

args <- commandArgs(TRUE)

path = args[1]
out_path = args[2]


files = list.files(path, pattern="_rel-abundance.tsv", all.files=FALSE)
sample_ids = gsub("_rel-abundance.tsv","",files)

drops <- c("superkingdom","species.group","subspecies","species.subgroup","sample_id","sample")

first_file <- files[1][1]

abundance_table <- read.csv(paste0(path,"/",first_file), sep="\t")
sample_id <- gsub("_rel-abundance.tsv","",first_file)
sample <- gsub("\\.._rel-abundance.tsv","",first_file)
abundance_table$sample <- sample
abundance_table$sample_id <- sample_id
names(abundance_table)[names(abundance_table) == "abundance"] <- paste0(sample_id,"_abundance")
names(abundance_table)[names(abundance_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]
abundance_table[is_empty(abundance_table)] = NA
for(f in files[-1]) {
  print("Entering loop")
  print("File")
  print(f)
  temp_table <- read.csv(paste0(path,"/",f), sep="\t")
  print("temp_table1")
  print(temp_table)
  sample_id <- gsub("_rel-abundance.tsv","",f)
  print("sample_id")
  print(sample_id)
  sample <- gsub("_rel-abundance.tsv","",f)
  print("sample")
  print(sample)
  abundance_table$sample <- sample
  print("abundance_table1")
  print(abundance_table)
  temp_table$sample_id <- sample_id
  print("temp_table2")
  print(temp_table)
  names(temp_table)[names(temp_table) == "abundance"] <- paste0(sample_id,"_abundance")
  print("temp_table3")
  print(temp_table)
  names(temp_table)[names(temp_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
  print("temp_table4")
  print(temp_table)
  temp_table <- temp_table[ , !(names(temp_table) %in% drops)]
  print("temp_table5")
  print(temp_table)
  temp_table[is_empty(temp_table)] = NA
  print("temp_table6")
  print(temp_table)
  abundance_table <- abundance_table %>% 
    full_join(temp_table, by = c("tax_id"="tax_id", 
                           "species"="species",
                           "genus"="genus",
                           "family"="family",
                           "order"="order",
                           "class"="class",
                           "phylum"="phylum",
                           "clade"="clade"))
}
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]
print(abundance_table)
abundance_table[abundance_table==""] <- "Unknown"
abundance_table <- abundance_table %>%
  select("tax_id", "species", "genus", "family", "order","class","phylum","clade",everything())
write.csv(abundance_table, file=paste0(out_path,"/full_abundance_table.csv"))


###### Species
abundance_table <- read.csv(paste0(path,"/",first_file), sep="\t")
sample_id <- gsub("_rel-abundance.tsv","",first_file)
sample <- gsub("\\.._rel-abundance.tsv","",first_file)
abundance_table$sample <- sample
abundance_table$sample_id <- sample_id
abundance_table <- abundance_table %>% 
  group_by(species) %>% 
  summarise(estimated.counts = round(sum(estimated.counts),digits=0),
            abundance = sum(abundance))
abundance_table[is_empty(abundance_table)] = NA
names(abundance_table)[names(abundance_table) == "abundance"] <- paste0(sample_id,"_abundance")
names(abundance_table)[names(abundance_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]

for(f in files[-1]) {
  temp_table <- read.csv(paste0(path,"/",f), sep="\t")
  sample_id <- gsub("_rel-abundance.tsv","",f)
  sample <- gsub("\\.._rel-abundance.tsv","",f)
  abundance_table$sample <- sample
  temp_table$sample_id <- sample_id
  temp_table <- temp_table[ , !(names(temp_table) %in% drops)]
  temp_table <- temp_table %>% 
    group_by(species) %>% 
    summarise(estimated.counts = round(sum(estimated.counts),digits =0),
              abundance = sum(abundance))
  names(temp_table)[names(temp_table) == "abundance"] <- paste0(sample_id,"_abundance")
  names(temp_table)[names(temp_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
  temp_table[is_empty(temp_table)] = NA
  abundance_table <- abundance_table %>% 
    full_join(temp_table, by = c("species"="species"))
}
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]
abundance_table[is_empty(abundance_table)] = NA
abundance_table[is.na(abundance_table)] <- 0
abundance_table[abundance_table==""] <- "Unknown"
write.csv(abundance_table, file=paste0(out_path,"/species_abundance_table.csv"))
##### Genus
abundance_table <- read.csv(paste0(path,"/",first_file), sep="\t")
sample_id <- gsub("_rel-abundance.tsv","",first_file)
sample <- gsub("\\.._rel-abundance.tsv","",first_file)
abundance_table$sample <- sample
abundance_table$sample_id <- sample_id
abundance_table <- abundance_table %>% 
  group_by(genus) %>% 
  summarise(estimated.counts = round(sum(estimated.counts),digits =0),
            abundance = sum(abundance))
abundance_table[is_empty(abundance_table)] = NA
names(abundance_table)[names(abundance_table) == "abundance"] <- paste0(sample_id,"_abundance")
names(abundance_table)[names(abundance_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]

for(f in files[-1]) {
  temp_table <- read.csv(paste0(path,"/",f), sep="\t")
  sample_id <- gsub("_rel-abundance.tsv","",f)
  sample <- gsub("\\.._rel-abundance.tsv","",f)
  abundance_table$sample <- sample
  temp_table$sample_id <- sample_id
  temp_table <- temp_table[ , !(names(temp_table) %in% drops)]
  temp_table <- temp_table %>% 
    group_by(genus) %>% 
    summarise(estimated.counts = round(sum(estimated.counts),digits =0),
              abundance = sum(abundance))
  names(temp_table)[names(temp_table) == "abundance"] <- paste0(sample_id,"_abundance")
  names(temp_table)[names(temp_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
  temp_table[is_empty(temp_table)] = NA
  abundance_table <- abundance_table %>% 
    full_join(temp_table, by = c("genus"="genus"))
}
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]
abundance_table[is_empty(abundance_table)] = NA
abundance_table[is.na(abundance_table)] <- 0
abundance_table[abundance_table==""] <- "Unknown"
write.csv(abundance_table, file=paste0(out_path,"/genus_abundance_table.csv"))

##### Family

abundance_table <- read.csv(paste0(path,"/",first_file), sep="\t")
sample_id <- gsub("_rel-abundance.tsv","",first_file)
sample <- gsub("\\.._rel-abundance.tsv","",first_file)
abundance_table$sample <- sample
abundance_table$sample_id <- sample_id
abundance_table <- abundance_table %>% 
  group_by(family) %>% 
  summarise(estimated.counts = round(sum(estimated.counts),digits =0),
            abundance = sum(abundance))
abundance_table[is_empty(abundance_table)] = NA
names(abundance_table)[names(abundance_table) == "abundance"] <- paste0(sample_id,"_abundance")
names(abundance_table)[names(abundance_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]

for(f in files[-1]) {
  temp_table <- read.csv(paste0(path,"/",f), sep="\t")
  sample_id <- gsub("_rel-abundance.tsv","",f)
  sample <- gsub("\\.._rel-abundance.tsv","",f)
  abundance_table$sample <- sample
  temp_table$sample_id <- sample_id
  temp_table <- temp_table[ , !(names(temp_table) %in% drops)]
  temp_table <- temp_table %>% 
    group_by(family) %>% 
    summarise(estimated.counts = round(sum(estimated.counts),digits =0),
              abundance = sum(abundance))
  names(temp_table)[names(temp_table) == "abundance"] <- paste0(sample_id,"_abundance")
  names(temp_table)[names(temp_table) == "estimated.counts"] <- paste0(sample_id,"_estimated.counts")
  temp_table[is_empty(temp_table)] = NA
  abundance_table <- abundance_table %>% 
    full_join(temp_table, by = c("family"="family"))
}
abundance_table <- abundance_table[ , !(names(abundance_table) %in% drops)]
abundance_table[is_empty(abundance_table)] = NA
abundance_table[is.na(abundance_table)] <- 0
abundance_table[abundance_table==""] <- "Unknown"
write.csv(abundance_table, file=paste0(out_path,"/family_abundance_table.csv"))
