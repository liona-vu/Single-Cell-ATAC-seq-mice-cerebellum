
##### PART 1. LOAD IN LIBRARIES, FILES, AND EDA ####
#load in libraries
library(ArchR)
library(parallel)

#set parameters to ensure reproducibility
addArchRLocking(locking = TRUE)
set.seed(1)

#load in annotated RDS file, file can be found with the following link:
#https://zenodo.org/records/6097549
anno_proj <- readRDS(file="E10_13_final.rds")

#perform some EDA
class(anno_proj)
#is a Seurat object [1] "Seurat", attr(,"package"), [1] "Seurat"

table(anno_proj$cellType1)
# seems like there are 25 different cell types such as CN.Th, GABA.Pre

sum(table(anno_proj$cellType1))
#19372 cells in total

#get file data
input_files <- c("fragment_files/GSM5393630_E12_snATAC_fragments.tsv.gz",
                 "fragment_files/GSM5393631_E13_snATAC_fragments.tsv.gz",
                 "fragment_files/GSM5393632_E14_snATAC_fragments.tsv.gz" )

#add names to the fragment files
names(input_files) <- c("GSM5393630_E12_snATAC", "GSM5393631_E13_snATAC", "GSM5393632_E14_snATAC")

#check if names are properly added
input_files

#                   GSM5393630_E12_snATAC 
#"GSM5393630_E12_snATAC_fragments.tsv.gz" 
#                   GSM5393631_E13_snATAC 
#"GSM5393631_E13_snATAC_fragments.tsv.gz" 
#                   GSM5393632_E14_snATAC 
#"GSM5393632_E14_snATAC_fragments.tsv.gz" 

#install mice reference genome
# install.packages("/home/lionavu/projects/def-itobias/BINF_6999/R_libs/BSgenome.Mmusculus.UCSC.mm10_1.4.3.tar.gz", 
                #  lib = "/home/lionavu/projects/def-itobias/BINF_6999/R_libs", 
                #  repos = NULL, 
                #  type = "source")

##### PART 2. CREATE ARROW FILES ####
#load mice reference package from computer
library(BSgenome.Mmusculus.UCSC.mm10, 
        lib.loc = "/home/rstudio/R_libs")

#Set the ArchR default genome to mm10 for mice 
addArchRGenome("mm10")

#Set parallel threads to capability to 1
addArchRThreads(threads = 1)

#create arrow files with default parameters 
ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

#Determine doublet scores in arrow file, will be used for doublet removal later
doub_scores <- addDoubletScores(input = ArrowFiles,
                                k = 10,
                                knnMethod = "KNN",
                                LSIMethod = 1)

#create ArchR project
proj_atacseq <- ArchRProject(ArrowFiles = ArrowFiles, 
                             outputDirectory = "scatac_seq_proj",
                             copyArrows = TRUE)
#save ArchR project
proj_atacseq <- saveArchRProject(ArchRProj = proj_atacseq, 
                                 outputDirectory = "scatac_seq_proj",
                                 load = TRUE)

##### PART 3. QUALITY CONTROL #####
#Load in previous ArchRproject
proj_atacseq <- loadArchRProject(path= "./scatac_seq_proj")

#check for metadata and see if doublets were added
head(proj_atacseq@cellColData)

#look at available matrices
getAvailableMatrices(proj_atacseq)

#check samples names and TSS enrichment
head(proj_atacseq$Sample)
head(proj_atacseq$TSSEnrichment)
table(proj_atacseq@cellColData$Sample)

#Generate simplified sample names for plotting 
names <- as.character(proj_atacseq@cellColData$Sample)
names <- gsub("GSM539363[0-2]_", "", names)
names <- gsub("_snATAC", "", names)

#Generate a new column for the simplified sample names 
proj_atacseq@cellColData$ShortSample <- names  
head(proj_atacseq@cellColData) #new column is here
head(row.names(proj_atacseq@cellColData))

#Now we have to check for the quality of the cells, pre-doublet removal
# plot violin plots to check for TSS enrichment scores for each samples pre doublet removal
plot_violin_tss <- plotGroups(ArchRProj = proj_atacseq,
                         groupBy = "ShortSample",
                         colorBy = "cellColData",
                         name = "TSSEnrichment",
                         plotAs = "violin", 
                         alpha = 0.5, #alpha is for shading
                         baseSize = 15,
                         discreteSet = "calm")

#Rename x-axis label
plot_violin_tss <- plot_violin_tss + xlab(label = "Timepoints") +
 theme(axis.text.x = element_text(angle = 0, hjust = 0.25))

#      legend.position = "none") +
#  scale_fill_viridis(discrete = TRUE, alpha = 0.5, option = "C") +
#  scale_color_viridis(discrete = TRUE)

plotPDF(plot_violin_tss, ArchRProj = proj_atacseq, name = "plot_violin_tss_pre",
        height = 5, width = 5, addDOC = FALSE)

#Check for unique nuclear fragments pre doublet removal
plot_violin_nfrag <- plotGroups(ArchRProj = proj_atacseq,
                                groupBy = "ShortSample",
                                colorBy = "cellColData",
                                name="log(nFrags)",
                                plotAs = "violin",
                                discreteSet = "calm",
                                baseSize = 15,
                                alpha = 0.5)

plot_violin_nfrag <- plot_violin_nfrag + xlab(label = "Sample Timepoints")

plotPDF(plot_violin_nfrag, name="plot_violin_unique_fragments", ArchRProj = proj_atacseq,
        width = 5, height = 5)

#filter for doublets 
filt_proj_atacseq <- filterDoublets(ArchRProj = proj_atacseq)

#Check number of cells pre and post doublet removal
proj_atacseq@cellColData@nrows #pre 26709
filt_proj_atacseq@cellColData@nrows #post 24276

#plot TSS scores after doublet removal
plot_violin_tss_post <- plotGroups(ArchRProj = filt_proj_atacseq,
                                   groupBy = "ShortSample",
                                   colorBy = "cellColData",
                                   name = "TSSEnrichment",
                                   plotAs = "violin",
                                   discreteSet = "calm",
                                   alpha = 0.5,
                                   baseSize = 15)

plot_violin_tss_post <- plot_violin_tss_post + xlab(label = "Sample Timepoints")

#save plot
plotPDF(plot_violin_tss_post, name = "plot_violin_tss_post", ArchRProj = proj_atacseq,
        width = 5, height = 5, addDOC = FALSE)

#plot unique fragments after doublet removal
plot_violin_nfrags_post <- plotGroups(ArchRProj = filt_proj_atacseq,
                                   groupBy = "ShortSample",
                                   colorBy = "cellColData",
                                   name = "log(nFrags)",
                                   plotAs = "violin",
                                   discreteSet = "calm",
                                   alpha = 0.5,
                                   baseSize = 15)

plot_violin_nfrags_post <- plot_violin_nfrags_post + xlab(label = "Sample Timepoints")

plotPDF(plot_violin_nfrags_post, name = "plot_violin_nfrags_post", ArchRProj = proj_atacseq,
        width = 5, height = 5, addDOC = FALSE)

#grabs archr calm palette, and rename each colour with each sample names
my_palette <- paletteDiscrete(values = proj_atacseq@cellColData$ShortSample, set = "calm")

#calculate and plot fragment distributions of pre doublet removal
frag_dist_pre <- plotFragmentSizes(ArchRProj = proj_atacseq,
                                   groupBy = "ShortSample",
                                   pal = my_palette)

#Save plot
plotPDF(frag_dist_pre, name = "fragment_sizes_plot_pre", ArchRProj = proj_atacseq, 
        height = 5, width =5, addDOC = FALSE)

#calculate and plot fragment size distributions of post doublet removal
#frag_dist_post <- plotFragmentSizes(ArchRProj = filt_proj_atacseq, 
 #                                   groupBy = "ShortSample",
  #                                  pal = my_palette)

#Save plot
#plotPDF(frag_dist_post, name = "fragment_sizes_plot_post", ArchRProj = proj_atacseq, 
#        height = 5, width =5, addDOC = FALSE)

#plot tss enrichment profiles pre doublet
tss_enrich_pre <- plotTSSEnrichment(ArchRProj = proj_atacseq,
                                    groupBy = "ShortSample",
                                    pal = my_palette)
#Save plot
plotPDF(tss_enrich_pre, name = "tss_enrichment_profile_plot_pre", ArchRProj = proj_atacseq, 
        height = 5, width =8, addDOC = FALSE)

#plot tss enrichment profiles post doublet
#tss_enrich_post <- plotTSSEnrichment(ArchRProj = filt_proj_atacseq,
                        #            groupBy = "ShortSample",
                         #           pal = my_palette)
#Save plot
#plotPDF(tss_enrich_post, name = "tss_enrichment_profile_plot_post", ArchRProj = filt_proj_atacseq, 
#        height = 5, width =5, addDOC = FALSE)


##### PART 4. DIMENSIONALITY REDUCTION WITH LSI ####
#Perform iterative LSI for post doublet removal
filt_proj_atacseq_2 <- addIterativeLSI(ArchRProj = filt_proj_atacseq,
                                  useMatrix = "TileMatrix", 
                                  name = "IterativeLSI", 
                                  iterations = 2, 
                                  clusterParams = list(resolution = c(0.2), 
                                                       sampleCells = 10000, 
                                                       n.start = 10), 
                                                       varFeatures = 25000, 
                                                       dimsToUse = 1:30)

#plot elbow plot to see if 30 dims is enough, grabs the LSI iterations and calculates standard deviations
rd <- getReducedDims(ArchRProj = filt_proj_atacseq_2, reducedDims = "IterativeLSI", dimsToUse = 1:30)
sd <- colSds(rd)

elbow_plot <- ggplot(data.frame(dims = 1:30, stdev = sd[1:30])) + 
  geom_point(mapping = aes_string(x = 'dims', y = 'stdev')) + 
  labs(x = "LSI",y="Standard Deviation") +
  theme_minimal()

#why is the first point so low?

ggsave(filename = "elbow_plot.png", plot = elbow_plot, path = "/home/rstudio/scatac_seq_proj/Plots")










