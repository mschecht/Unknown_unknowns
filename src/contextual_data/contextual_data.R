### Merging ORF and contextual data sets

library("tidyverse")

# concatenate all headers from the ORF data
# mutate -> adds a column (name of column = "value of column")
# rename -> renames column name to x
ORF_headers <- bind_rows(read_table("~/Desktop/msc/contextualdata/ORF_headers/malaspina_parsed_headers.txt", col_names = FALSE) %>%
    mutate(Project = "Malaspina") %>%
    rename(label = X1),
  read_tsv("~/Desktop/msc/contextualdata/ORF_headers/TARA_parsed_headers.txt", col_names = FALSE) %>%
    mutate(Project = "TARA") %>%
    rename(label = X1),
  read_tsv("~/Desktop/msc/contextualdata/ORF_headers/OSD_parsed_headers.txt", col_names = FALSE) %>%
    mutate(Project = "OSD") %>%
    rename(label = X1),
  read_tsv("~/Desktop/msc/contextualdata/ORF_headers/GS_parsed_headers.txt", col_names = FALSE) %>%
    mutate(Project = "GOS") %>%
    rename(label = X1))

# load conextual data
contextual_data <- read_tsv("~/Desktop/msc/contextualdata/sample_data/compiled_metadata.txt", col_names = TRUE)

# count number of samples for each project in contextual data
contextual_data %>%
  group_by(Project) %>%
  count()

# count number of samples for each project in ORF data
ORF_headers %>%
  group_by(Project) %>%
  count()

# take the contextual data, filter out the Malaspina section, 
#rename the "sample" column to "label" and 
#merge it with the ORF_headers (malaspina section), 
#then filter out the one that are not matching

contextual_data %>%
  filter(Project == "Malaspina") %>%
  rename(label = Sample) %>%
  full_join(ORF_headers %>% filter(Project == "Malaspina"), by = "label") %>% 
  filter(is.na(Project.y) | is.na(Project.x)) %>%
  View

contextual_data %>%
  filter(Project == "TARA") %>%
  rename(label = Sample) %>%
  full_join(ORF_headers %>% filter(Project == "TARA"), by = "label") %>% 
  filter(is.na(Project.y) | is.na(Project.x)) %>%
  View

contextual_data %>%
  filter(Project == "GOS") %>%
  rename(label = Sample) %>%
  full_join(ORF_headers %>% filter(Project == "GOS"), by = "label") %>% 
  filter(is.na(Project.y) | is.na(Project.x)) %>%
  View

contextual_data %>%
  filter(Project == "OSD") %>%
  rename(label = Sample) %>%
  full_join(ORF_headers %>% filter(Project == "OSD"), by = "label") %>% 
  filter(is.na(Project.y) | is.na(Project.x)) %>%
  View
