library("BDgraph") # check overlap of bulk soil core ASVs
library("parallel")
library("phyloseq")


# load data
load("./kathe/phyloseq_objects.RData")

load(file = "./kathe/model_borutized.RData")

# make a new ps object, that only has AVSs classified as important by boruta
network_input_ps<-prune_taxa(model_borutized$coefnames,physeq_goodSamples_rarefied)


# transform phyloseq objects otus into matrixes
matrix_input<-t(as(otu_table(network_input_ps), "matrix"))

# transform phyloseq objects metadata into df
metadata_input<-as(sample_data(network_input_ps), "data.frame")

# force some metadata into  factor
metadata_input$Treatment<-as.factor(metadata_input$Treatment)

#calculate network, removing spurious edges caused by treament metadata
BD_out_750k<-bdgraph.dw(data = matrix_input, 
                   x = metadata_input,
                   algorithm = "bdmcmc",
                   formula = y ~  Treatment,
                   cores = "all",
                   iter_bdw = 750000,
                   iter = 150000,
                   save = TRUE)

# save output externallu
save(BD_out_750k, file = "./kathe/BD_out_750k.RData")
