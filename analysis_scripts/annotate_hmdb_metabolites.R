library(metaboliteIDmapping)
library(AnnotationHub)
ah <- AnnotationHub()
datasets <- query( ah, "metaboliteIDmapping")

datasets[4]



options(timeout = 50000)
annotation_data_1 <- ah[["AH79817"]]
annotation_data_2 <- ah[["AH83115"]]
annotation_data_3 <- ah[["AH91792"]]

annotation_data_3.1 <- annotation_data_3[which(!is.na(annotation_data_3$HMDB)&!is.na(annotation_data_3$KEGG)),]
met2dis <- openxlsx::read.xlsx('L:/Regan Data/HMDB Data/hmdb_disease_data_v2.xlsx')
# dim(unique(met2dis))
length(unique(met2dis$Disease))
met2dis <- openxlsx::read.xlsx('L:/Regan Data/HMDB Data/hmdb_disease_data_v2.xlsx')|>
  dplyr::mutate(KEGGID = annotation_data_3.1$KEGG[match(HMDB_ID,annotation_data_3.1$HMDB)])
met.KEGG <- unique(met2dis$KEGGID)
met.KEGG <- met.KEGG[!is.na(met.KEGG)]
length(met.KEGG)
openxlsx::write.xlsx(met2dis,'L:/Regan Data/HMDB Data/met2dis_Ref.xlsx')
