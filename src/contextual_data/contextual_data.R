library("tidyverse")

# concatenate all headers from the ORF data
ORF_headers <- bind_rows(read_tsv("~/Desktop/msc/metadata/ICM_parsed.headers", col_names = FALSE) %>%
  mutate(Project = "Malaspina") %>%
  rename(label = X1),
  read_tsv("~/Desktop/msc/metadata/TARA_parsed.headers", col_names = FALSE) %>%
    mutate(Project = "TARA") %>%
    rename(label = X1),
read_tsv("~/Desktop/msc/metadata/OSD_parsed.headers", col_names = FALSE) %>%
  mutate(Project = "OSD") %>%
  rename(label = X1),
read_tsv("~/Desktop/msc/metadata/GS_parsed.headers", col_names = FALSE) %>%
  mutate(Project = "GOS") %>%
  rename(label = X1))

# load conextual data
contextual_data <- read_tsv("~/Desktop/msc/metadata/compiled_metadata.txt", col_names = TRUE)

# count number of samples for each project in contextual data
contextual_data %>%
  group_by(Project) %>%
  count()

# count number of samples for each project in ORF data
ORF_headers %>%
  group_by(Project) %>%
  count()

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
