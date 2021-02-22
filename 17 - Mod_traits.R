
## load the correlation tables

tcga_traits <- read.table("Correlation_matrix_TCGA_traits_eight.txt", sep = "\t", header = TRUE, row.names = 1)
meta_traits <- read.table("Correlation_matrix_METABRIC_traits_eight.txt", sep = "\t", header = TRUE, row.names = 1)

tcga_traits <- as.matrix(tcga_traits)
meta_traits <- as.matrix(meta_traits)

mat_col <- data_frame(Colours = c("black", "blue3","brown","forestgreen","magenta3","red3","tan3","yellow3"))
rownames(mat_col) <- rownames(tcga_traits)


library(colorspace)
##view colour palettes
hsv_palettes(plot = TRUE)
##create a colour palette for later use
q4 <- diverging_hsv(50)
q4 <- diverging_hcl(50, "Blue_Red")

#diverging red-blue colors
swatchplot(
  diverging_hsv(7),
  desaturate(diverging_hsv(7)),
  diverging_hcl(7, c = 100, l = c(50, 90)),
  desaturate(diverging_hcl(7, c = 100, l = c(50, 90))),
  nrow = 2
)


#colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(50)
pdf("Module_Traits_META.pdf",height= 6, width=12, family = "Gill Sans MT")
par(oma = c(1,1,1,1))
pheatmap(
  mat               = meta_traits,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color             = q4,
  border_color      = "grey",
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_colors = mat_col$Colours,
  treeheight_row = 0, 
  treeheight_col = 0,
  drop_levels       = TRUE,
  angle_col = 45,
  fontsize          = 12,
  main              = "Correlation between modules and key clinical traits: METABRIC"
)
dev.off()


mod_traits <- read.table("Module_Trait_Final.txt", sep = "\t", header= TRUE)

levels(mod_traits$Trait)
mod_traits$Trait <- factor(mod_traits$Trait,levels = c("ESR1", "PGR", "AR", "FOXM1", "ERBB2", "ER Metagene", "CIN Metagene"))

p <- ggplot(mod_traits, aes(Module, reorder(Trait, desc(Trait)))) + 
      geom_tile(aes(fill = Correlation), colour = "black") + 
      scale_fill_gradient2() +
      theme_minimal(base_size = 20) +
      facet_wrap(~Dataset)
p

### module vs module
mod_mod <- read.table("Mod_Mod_Final.txt", sep = "\t", header= TRUE)

levels(mod_mod$Module_1)
levels(mod_mod$Module_2)
mod_traits$Trait <- factor(mod_traits$Trait,levels = c("ESR1", "PGR", "AR", "FOXM1", "ERBB2", "ER Metagene", "CIN Metagene"))

p <- ggplot(mod_mod, aes(Module_1, Module_2)) + 
  geom_tile(aes(fill = Correlation), colour = "black") + 
  scale_fill_gradient2() +
  theme_minimal(base_size = 20) +
  facet_wrap(~Dataset)
p
