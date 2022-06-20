################################################################################
################################################################################

# Define functions 

################################################################################
################################################################################


####=========================================================================###
####=========================================================================###

# LOAD PACKAGES

####=========================================================================###
####=========================================================================###


load_packages=function(){
  require(R.utils) # referenced 
  require(binom) # referenced 
  require(data.table) # referenced
  require(stringr) # referenced
  require(formattable) # not referenced yet - used for table plot in genetic hetereogeneity
  require(tidyr) # referenced
  require(grid) # referenced as part of R - used for plotting heterogeneity plot 
  require(gridExtra) # referenced used for heterogeneity plot
  require(meta) # referenced
  require(DT) # referenced
  require(plyr) # referenced (needed for round any)
  require(dplyr) # referenced
  require(ggvis) #referenced
  library(kinship2)
  #require(CoSeg)
  #install.packages("CoSeg", repos="http://R-Forge.R-project.org")
  require(berryFunctions) # referenced needed for path plot 

}

####=========================================================================###
####=========================================================================###

# DATASET MANAGEMENT FUNCTIONS

####=========================================================================###
####=========================================================================###


population_people_dataset_function=function(lit_review_pen){
  # Function to pull out population information individuals from lit_review_pen and split this based on '|'
  population_people_dataset=subset(lit_review_pen[!is.na(lit_review_pen$population_carriers_count) &  lit_review_pen$population_carriers_count>0,],select=c("gene","HGVS","population_carriers_primary_phenotype","population_carriers_continent","population_carriers_nationality","population_carriers_family_history","population_carriers_pmid")) 
  population_people_dataset=separate_rows(population_people_dataset,c("population_carriers_primary_phenotype","population_carriers_continent","population_carriers_nationality","population_carriers_family_history","population_carriers_pmid"),sep="\\|")
  population_people_dataset$population_carriers_family_history<-gsub("F","Familial",population_people_dataset$population_carriers_family_history)
  population_people_dataset$population_carriers_family_history<-gsub("S","Sporadic",population_people_dataset$population_carriers_family_history)
  return(population_people_dataset)
}

render.acmg_pathogenicity_plots=function(input_variant,lit_review_pen,plot_height){
  print(plot_height)
  # create blank plot arbitrarily with limits 100,100
  print(input_variant)
  plot(1,1,     
    lwd=plot_height/100,
    col="white",
    type="l",
    xaxt="n",
    yaxt="n",
    las=1,
    xlab="",
    ylab="",
    main=paste(input_variant,lit_review_pen$acmg_literal[!is.na(lit_review_pen$HGVS) & lit_review_pen$HGVS==input_variant],sep=": "),
    cex.lab=plot_height/600,
    cex.main=plot_height/600,
    col.main=nicegrey,
    col.lab="white",
    frame=FALSE,
    xlim=c(0,100),
    ylim=c(0,100)
    )
  text_start_pos <- 85
  box_width <- 3
  gap <- 1
  pen_left<- 46.5
  ben_left<- 50.5 
  text(pen_left+box_width/2,text_start_pos,"pathogenic",col=nicegrey,adj=0,srt=90,cex=plot_height/600,xpd=NA)
  text(ben_left+box_width/2,text_start_pos,"benign",col=nicegrey,adj=0,srt=90,cex=plot_height/600,xpd=NA)
  category_explanations<- c("very strong / stand alone","strong","moderate","supporting")
  path_colours<-list(nicered,niceorange,darkgreen,lightgreen)
  #ben_colours<-list(reallydarkblue,nicedarkblue,nicelightblue,reallylightblue)
  ben_colours<-list(reallydarkblue,nicedarkblue,"white",nicelightblue)
  j <- text_start_pos-gap # upper y for first boxes
  for (i in 1:4){
    roundedRect(pen_left,j-box_width,pen_left + box_width ,j,border="white",col=path_colours[[i]], bothsame=TRUE, aspcorrect=TRUE,xpd=NA)
    roundedRect(ben_left,j-box_width,ben_left + box_width ,j,border="white",col=ben_colours[[i]], bothsame=TRUE, aspcorrect=TRUE,xpd=NA)
    text(44,j-box_width/2,category_explanations[i],col=nicegrey,adj=1,cex=plot_height/600,xpd=NA)
    j=j-box_width-gap
  }
  pathogenic_text <- c("null variant |","same amino acid change |","de novo |","functional studies |","prevalence v controls |","hot spot |","prevalence v database |","trans/cis testing |","in-frame INDEL |","novel missense at same position |","de novo |","segregation |","variant spectrum |","in silico |","phenotype |","reputable source |")
  pathogenic_category <- c("PVS1","PS1","PS2","PS3","PS4","PM1","PM2","PM3","PM4","PM5","PM6","PP1","PP2","PP3","PP4","PP5")
  benign_text <- c(" | database frequency"," | prevalence v controls"," | observed in healthy adult"," | functional studies"," | segregation"," | variant spectrum"," | trans/cis testing"," | in-frame INDEL"," | in silico"," | alternate locus observations"," | reputable source"," | synonymous variants")
  benign_category <- c("BA1","BS1","BS2","BS3","BS4","BP1","BP2","BP3","BP4","BP5","BP6","BP7")
  pathogenic_colours <- as.list(strsplit(lit_review_pen$path_colours[lit_review_pen$HGVS==input_variant], ",")[[1]])
  benign_colours <- as.list(strsplit(lit_review_pen$ben_colours[lit_review_pen$HGVS==input_variant], ",")[[1]])
  j<-text_start_pos - (6*box_width + 4*gap)
  k<-j-(2*box_width+2*gap)
  for (i in 1:16){
    roundedRect(pen_left,j-box_width,pen_left + box_width,j,border="white",col=pathogenic_colours[[i]], bothsame=TRUE, aspcorrect=TRUE,xpd=NA)
    text(38,j-box_width/2,pathogenic_text[i],col=nicegrey,adj=1,cex=plot_height/600,xpd=NA)
    text(39,j-box_width/2,pathogenic_category[i],col=nicegrey,adj=0,cex=plot_height/600,xpd=NA)
    if (i <= length(benign_category)){
      roundedRect(ben_left,k-box_width,ben_left+box_width,k,border="white",col=benign_colours[[i]], bothsame=TRUE, aspcorrect=TRUE,xpd=NA)
      text(51+box_width+1,k-box_width/2,benign_category[i],col=nicegrey,adj=0,cex=plot_height/600,xpd=NA)
      text(56+box_width+1,k-box_width/2,benign_text[i],col=nicegrey,adj=0,cex=plot_height/600,xpd=NA)
    }
    j=j-box_width-gap
    k=k-box_width-gap
  }
}
# ####=========================================================================###
# ####=========================================================================###

# # SUMMARY PLOT FUNCTIONS

# ####=========================================================================###
# ####=========================================================================###
gene_text=function(input_gene){
  gene_text<- gene_summary_text[[input_gene]]
  return(paste("<div style='color:grey'>",gene_text,sep=""))
}
gene_summary_text <- list(
  "SOD1"="Of 244 reported SOD1 variants, 49 are identified as P or LP causes of ALS. Variants are present in every exon and all are missense, with no B or LB missense variants observed (supplementary figure XXX A). These variants are a major global cause of FALS (11.0%: 95% CI 9.7-12.5%) and a minor cause of SALS (0.9%: 95% CI 0.8-1.2%) (supplementary figure xxx); and while several variants are rare globally, they explain large proportions of cases in a particular region and thus have significant within- or between-continent geographic heterogeneity (supplementary figure xxxx). ALS patients with P or LP variants in SOD1 have an earlier median AOO (48.5: 95% CI 46.5-50) than non-SOD1 variants (55: 95% CI 55-56) (supplementary figure XXX), although this does not appear to be the case for all SOD1 P or LP variants. Carriers of SOD1 VUS also present with moderately early AOO, indicating the presence of further P variants currently with insufficient supporting evidence. While three ALS-FTD SOD1 variant carriers are reported, two of these individuals carry a LB or intronic variant, indicating an alternative genetic cause.",
  "TARDBP"="10 P or LP missense variants are found in the C-terminal glycine rich final exon of TARDBP (supplementary figure XXX ). These variants present with phenotypes spanning the ALS-FTD spectrum; although this may be variant dependent (supplementary figure XXX). P and LP TARDBP variants are major global causes of FALS (4.0%: 95% CI 3.2%-5.0%) and FFTD (2.0% 95% CI 1.0-3.8%) and are also observed in SALS (0.9% 95% CI 0.7-1.1%) and SFTD (0.2% 95% CI 0.03-0.9%). Geographic heterogeneity is observed for TARDBP:c.1144G>A(p.[A382T]) which is present in a large proportion of Sardinian FALS (32% 95% CI 23.0-42.1%) and SALS (20.4% 95% CI 15.8-25.6%) cases but is virtually absent throughout the rest of the world."
  )
scale_bar_plot_linear=function(left,right,text_size,plot_height){
    # Set text size 
  # Plot side by side linear scale bars 
  # Horizontal Lines
  lines(x=c(left,right),y=c(10,10),col=nicegrey,lty=1,lwd=0.6*plot_height/1200)
  # Vertical Lines 
  for (number in seq(0,1,0.2)){
    lines(x=c(left+number*(right-left),left+number*(right-left)),y=c(10,8),col=nicegrey,lty=1,lwd=0.6*plot_height/1200)
  }
  # text for scale bar 
  text(c(left,left+0.2*(right-left),left+0.4*(right-left),left+0.6*(right-left),left+0.8*(right-left),right),5,c("0","20","40","60","80","100%"),col=nicegrey,adj=0.5,cex=text_size*plot_height/1200)

}
log_position=function(numberofbreaks,value){

  # Function to convert a value to it's position on a log scale 
  100*(log10(value)+2)/numberofbreaks
}
scale_bar_plot_log=function(left,right,text_size){

  # Plot side by side log scale bars 
  # Horizontal Lines
  # Need to create axis break 
  # Set text size 
  lines(x=c(left,left+4),y=c(10,10),col=nicegrey,lty=1,lwd=0.6)
  lines(x=c(left+7,right),y=c(10,10),col=nicegrey,lty=1,lwd=0.6)
  # Vertical Lines 
  for (number in seq(0,1,0.25)){
    lines(x=c(left+number*(right-left),left+number*(right-left)),y=c(10,8),col=nicegrey,lty=1,lwd=0.6)
  }
  for (number in c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),seq(1,9,1),seq(10,90,10))){
    lines(x=c(left+log_position(4,number),left+log_position(4,number)),y=c(10,8),col=nicegrey,lty=1,lwd=0.6)
  }
  # text for scale bar 
  text(c(left,left+0.25*(right-left),left+0.5*(right-left),left+0.75*(right-left),right),5,c("0","0.1","1","10","100%") ,col=nicegrey,adj=0.5,cex=text_size)

}


render.summary.proportion_explained.plot=function(lit_review_pen,plot_height){
  par(mfrow=c(1,1), mai = c(0.1,0.1,0.1,0.1))
  # Create Blank Plot 
  blank_plot(375,95)
  # Set Starting Plotting Variables
  text_size <- 1.2
  gap_between   <- 5
  end_pos       <- 100
  top           <- 90
  bottom        <- 85
  als_ftd_gap   <- 20
  left.1        <- 75
  right.1       <- left.1 + 100
  left.2        <- right.1 + als_ftd_gap 
  right.2       <- left.2 + 100
  box_dimension <- 4.5*plot_height/1200
  # Set Starting Variables 
  regions_phenotype_selection="ALS"
  regions_pathogenicity_selection="1" # 1 is P variants, 2 is P or LP variants, 3 is all variants 
  regions_history_selection="Familial"
  # 12 Groups to Plot
  # Variables change as we cycle through each group
  for (number in 1:12){
    start_pos <- left.1
    # Reset Plotting Variables
    if (number %% 2 == 0 ){ #even numbers 
      top<-bottom-gap_between
      bottom<-top-gap_between
      regions_history_selection="Sporadic"
    }
    if (number==3 | number==9){
      top<-bottom-3*gap_between
      bottom<-top-gap_between
      regions_history_selection="Familial"
      regions_pathogenicity_selection="2"
    }
    if (number==5 | number==11){
      top<-bottom-3*gap_between
      bottom<-top-gap_between
      regions_history_selection="Familial"
      regions_pathogenicity_selection="3"
    }
    if (number>=7){
      start_pos <- left.2
      end_pos <- start_pos
    }
    if (number==7){
      top<-90
      bottom<-top-gap_between
      regions_history_selection="Familial"
      regions_pathogenicity_selection="1"
      regions_phenotype_selection="FTD"
  }
  
  # Get list of potential variants 
  # This is clunky as taken from the journALS data browser script so has unneccessary redundancy
  regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
  # If selection is familial or sporadic can calculate these directly 
  if (regions_history_selection=="Familial" | regions_history_selection=="Sporadic"){
    # Filter population individuals: phenotype of interest; history of interest ; variant matches 
    population_people_use    <- population_people_dataset[population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history==regions_history_selection & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
    # Filter population studies: phenotype of interest; history of interest 
    population_studies_use   <- population_studies[population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History==regions_history_selection & !is.na(population_studies$History),]
    # Collapse Files If there is a file to collapse 
    if (nrow(population_studies_use)>0){
      population_studies_use_collapse <- aggregate(population_studies_use$Count,by=list(gene=population_studies_use$Gene),FUN=sum)
      # Create blank column
      population_studies_use_collapse$Carrier_count<-0
      # Sum the carriers from each country 
      for (Gene in population_studies_use_collapse$gene){
        population_studies_use_collapse$Carrier_count <- ifelse(population_studies_use_collapse$gene==Gene,nrow(population_people_use[population_people_use$gene==Gene,]),population_studies_use_collapse$Carrier_count)
      }
      # Get the proportion + CI of carriers for each country 
      population_studies_use_collapse$proportion        <-binom.confint(x=population_studies_use_collapse$Carrier_count,n=population_studies_use_collapse$x,method='wilson')$mean
      population_studies_use_collapse$proportion.lower  <-binom.confint(x=population_studies_use_collapse$Carrier_count,n=population_studies_use_collapse$x,method='wilson')$lower
      population_studies_use_collapse$proportion.upper  <-binom.confint(x=population_studies_use_collapse$Carrier_count,n=population_studies_use_collapse$x,method='wilson')$upper
      # Filter file to pathogenic or LP genes 
      if (regions_pathogenicity_selection == as.numeric(1)) {
        population_studies_use_collapse=population_studies_use_collapse[population_studies_use_collapse$gene %in% path_genes,]
      } 
      else if ((regions_pathogenicity_selection == as.numeric(2))|(regions_pathogenicity_selection == as.numeric(3))) {
        population_studies_use_collapse=population_studies_use_collapse[population_studies_use_collapse$gene %in% p_or_lp_genes,]
      }
      # Reorder
      population_studies_use_collapse <- population_studies_use_collapse[with(population_studies_use_collapse,rev(order(proportion,x,gene))),]
    }
  } 
  # Calculate unexplained proportion
  unexplained<-100*(1-sum(population_studies_use_collapse$proportion))
  # Plot unexplained proportion 
  if (unexplained>0){
    if (number<7){ # i.e if ALS 
      end_pos <- left.1 + unexplained
    } 
    else {
      end_pos <- left.2 + unexplained
    }
    rect(start_pos,bottom,end_pos,top,col="grey70",border=NA)
    text(start_pos+3,(bottom+top)/2,paste(round(unexplained,2),"%",sep=""),col="white",adj=0,cex=text_size*plot_height/1200)
    start_pos<-end_pos
  }
  # Plot remaining proportions 
  if (exists("population_studies_use_collapse")==TRUE){
    number_of_boxes<-nrow(population_studies_use_collapse)
    for (i in 1:number_of_boxes){
      end_pos=start_pos+100*population_studies_use_collapse[i,]$proportion
      gene_use<- population_studies_use_collapse[i,]$gene
      colour<-vangogh_palette[match(gene_use,p_or_lp_genes)]
      rect(start_pos,bottom,end_pos,top,col=colour,border=NA)
      start_pos<-end_pos
    }
  }
}
# Add text annotations to plot 
text(right.2+3,c(87.5,77.5,57.5,47.5,27.5,17.5),c("Familial","Sporadic"),col=nicegrey,adj=0,cex=text_size*plot_height/1200,xpd=NA)
text(left.1-3,c(82.5,52.5,22.5),c("Pathogenic Variants","Pathogenic & Likely Pathogenic Variants","Reported Variants in \nPathogenic & Likely Pathogenic Genes"),col=nicegrey,adj=1,cex=text_size*plot_height/1200,xpd=NA)
text(c((left.1+right.1)/2,(left.2+right.2)/2),95,c("ALS","FTD"),col=nicegrey,adj=0.5,cex=text_size*plot_height/1200,xpd=NA)
# lines to left of plot
lines(x=c(left.1-2,left.1-2),y=c(90,75),col=nicegrey,lty=1,lwd=0.6*plot_height/1200,xpd=NA)
lines(x=c(left.1-2,left.1-2),y=c(60,45),col=nicegrey,lty=1,lwd=0.6*plot_height/1200,xpd=NA)
lines(x=c(left.1-2,left.1-2),y=c(30,15),col=nicegrey,lty=1,lwd=0.6*plot_height/1200,xpd=NA)
# Plot scale bars 
scale_bar_plot_linear(left.1,right.1,text_size,plot_height)
scale_bar_plot_linear(left.2,right.2,text_size,plot_height)
lines(x=c(right.2+30,right.2+30),y=c(87.5,17.5),col=nicegrey,lty=1,lwd=0.6*plot_height/1200,xpd=NA)

# plot grey no variant box
roundedRect(right.2+33,91,right.2+33+box_dimension,91+box_dimension,border="white",col="grey70", bothsame=TRUE, aspcorrect=TRUE,corfactor=6,xpd=NA)
text(right.2+42,91+box_dimension/2,"No Variant",col=nicegrey,adj=0,cex=text_size*plot_height/1200,xpd=NA)
# Plot boxes for remaining genes 
height=86
for (gene in unique(population_studies_use_collapse$gene[population_studies_use_collapse$gene %in% p_or_lp_genes])){
  colour<-vangogh_palette[match(gene,p_or_lp_genes)]
  roundedRect(right.2+33,height,right.2+33+box_dimension,height+box_dimension,border="white",col=colour, bothsame=TRUE, aspcorrect=TRUE,corfactor=6,xpd=NA)   
  text(right.2+42,height+box_dimension/2,gene,col=nicegrey,adj=0,cex=text_size*plot_height/1200,font=3,xpd=NA)
  height=height-5
}

}

# ####=========================================================================###
# ####=========================================================================###

# # AGE FUNCTIONS

# ####=========================================================================###
# ####=========================================================================###

# # Compute p value cutoff for age tests

age_p_value_cutoff=function(all_ages_df){
  # Two things to consider
  #   1:  What is the lowest p value possible from a given number of individuals (e.g. if we take the three individuals with lowest aoo what p value do we get)
  #   2:  Is this p-value lower than the Bonferroni corrected p value for the number of possible tests at that level (e.g. 27 different tests for each variant (phenotype x sex x family history) but there may only be two ages for certain tests so discount these)
  # For each variant carrier with a reported age
  for (number in seq(1:length(sort(all_ages_df$all_carriers_aoo)))){
    print("here")
    number_of_tests<-0
    # group 1
    # Find counts of all carriers for all variants
    df<-as.data.frame(table(all_ages_df$HGVS)) 
    number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
    # group 2
    for (phenotype in c("ALS","FTD")){
      df<-as.data.frame(table(all_ages_df$HGVS[all_ages_df$all_carriers_primary_phenotype==phenotype])) 
      number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
      for (sex in c("M","F")){
        df<-as.data.frame(table(all_ages_df$HGVS[all_ages_df$all_carriers_primary_phenotype==phenotype & all_ages_df$all_carriers_sex==sex])) 
        number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
        for(fam in c("F","S")){
          df<-as.data.frame(table(all_ages_df$HGVS[all_ages_df$all_carriers_primary_phenotype==phenotype & all_ages_df$all_carriers_sex==sex & all_ages_df$all_carriers_family_history==fam])) 
          number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
        }
      }
      for(fam in c("F","S")){
        df<-as.data.frame(table(all_ages_df$HGVS[all_ages_df$all_carriers_primary_phenotype==phenotype & all_ages_df$all_carriers_family_history==fam])) 
        number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
      }
    }
    for (sex in c("M","F")){
      df<-as.data.frame(table(all_ages_df$HGVS[all_ages_df$all_carriers_sex==sex])) 
      number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
      for(fam in c("F","S")){
        df<-as.data.frame(table(all_ages_df$HGVS[all_ages_df$all_carriers_sex==sex & all_ages_df$all_carriers_family_history==fam])) 
        number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
      }
    }
    for(fam in c("F","S")){
      df<-as.data.frame(table(all_ages_df$HGVS[all_ages_df$all_carriers_family_history==fam])) 
      number_of_tests=number_of_tests+nrow(df[df$Freq>=number,])
    }


    # max p-value at this level
    max_p<-0.05/number_of_tests
    # Find the minimum possible p-value from the data at this level (i.e. the n youngest age of onset compared to the rest)
    kruskal<-kruskal_test(sort(all_ages_df$all_carriers_aoo)[1:number],sort(all_ages_df$all_carriers_aoo)[(number+1):length(sort(all_ages_df$all_carriers_aoo))])
    if (as.numeric(kruskal) < as.numeric(max_p)){
      return(max_p)
    }
    # The variants which reach siginifance are:
    # Phe: All | Sex: All | Family History : All
      # "ATXN2:c.(CAG)n"
      # "C9orf72:c.-45+163GGGGCC[>24]"
      # "FUS:c.1486C>T(p.[R496*])"
      # "FUS:c.1512_1513delAG(p.[G505fs])"
      # "FUS:c.1564C>T(p.[R522C])"
      # "FUS:c.1577C>T(p.[P526L])"
      # "SOD1:c.272A>C(p.[D91A])"
      # "VAPB:c.166C>T(p.[P56S])"
    # Phe: ALS | Sex: All | Family History : All
      # "ATXN2:c.(CAG)n"
      # "C9orf72:c.-45+163GGGGCC[>24]"
      # "FUS:c.1486C>T(p.[R496*])"
      # "FUS:c.1564C>T(p.[R522C])"
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: FTD | Sex: All | Family History : All
      # "MAPT:c.1842T>G(p.[N614K])"
    # Phe: All | Sex: M | Family History : All
      # "C9orf72:c.-45+163GGGGCC[>24]"
      # "FUS:c.1486C>T(p.[R496*])"
      # "FUS:c.1564C>T(p.[R522C])"
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: All | Sex: F | Family History : All
      # "ATXN2:c.(CAG)n"
      # "C9orf72:c.-45+163GGGGCC[>24]"
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: All | Sex: All | Family History : S
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: All | Sex: All | Family History : F
      # "C9orf72:c.-45+163GGGGCC[>24]"
      # "FUS:c.1564C>T(p.[R522C])"
      # "GRN:c.813_816delCACT(p.[T272fs])"
      # "VAPB:c.166C>T(p.[P56S])"
    # Phe: ALS | Sex: M | Family History : All
      # "C9orf72:c.-45+163GGGGCC[>24]"
      # "FUS:c.1564C>T(p.[R522C])"
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: ALS | Sex: F | Family History : All
      # "ATXN2:c.(CAG)n"
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: ALS | Sex: S | Family History : All
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: ALS | Sex: All | Family History : F
      # "C9orf72:c.-45+163GGGGCC[>24]"
    # Phe: FTD | Sex: All | Family History : F
      # "MAPT:c.1842T>G(p.[N614K])"
    # Phe: All | Sex: M | Family History : F
      # "C9orf72:c.-45+163GGGGCC[>24]"
    # Phe: All | Sex: F | Family History : S
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: All | Sex: F | Family History : F
      # "C9orf72:c.-45+163GGGGCC[>24]"
    # Phe: ALS | Sex: M | Family History : F
      # "C9orf72:c.-45+163GGGGCC[>24]"
    # Phe: ALS | Sex: F | Family History : S
      # "FUS:c.1577C>T(p.[P526L])"
    # Phe: ALS | Sex: F | Family History : F
      # "C9orf72:c.-45+163GGGGCC[>24]"
    } 
  }

  min_age_p_value=function(variant_variant_browser,all_ages_df,age_p_cutoff){
    variant_ages<-all_ages_df[all_ages_df$HGVS==variant_variant_browser & is.na(all_ages_df$HGVS)==F,]
    all_ages<-all_ages_df[all_ages_df$HGVS!=variant_variant_browser& is.na(all_ages_df$HGVS)==F,]
    if (nrow(variant_ages[!is.na(variant_ages$all_carriers_aoo),])==0 ){
      this_list=list(
        "min_p_value"="NA",
        "min_phenotype"="NA",
        "min_sex"="NA",
        "min_fam"="NA",
        "above"="NA",
        "evidence"="no"
        )
      return(this_list)
    }
    min_p_value<-kruskal_test(variant_ages$all_carriers_aoo,all_ages$all_carriers_aoo)
    min_phenotype<-"All"
    min_sex<-"All"
    min_fam<-"All"
    for (phenotype in c("All","ALS","FTD")){
      if (phenotype != "All"){
        variant_ages_temp<-variant_ages[variant_ages$all_carriers_primary_phenotype==phenotype,]
        all_ages_temp<-all_ages[all_ages$all_carriers_primary_phenotype==phenotype,]
      }
      else {
        variant_ages_temp<-variant_ages
        all_ages_temp<-all_ages
      }
      for (sex in c("All","M","F")){
        if (sex != "All"){
          variant_ages_temp<-variant_ages_temp[variant_ages_temp$all_carriers_sex==sex,]
          all_ages_temp<-all_ages_temp[all_ages_temp$all_carriers_sex==sex,]
        }
        for (fam_hist in c("All","F","S")){
         if (fam_hist != "All"){
           variant_ages_temp<-variant_ages_temp[variant_ages_temp$all_carriers_family_history==fam_hist,]
           all_ages_temp<-all_ages_temp[all_ages_temp$all_carriers_family_history==fam_hist,]
         }
         if(nrow(variant_ages_temp[!is.na(variant_ages_temp$all_carriers_aoo),])>0 & nrow(all_ages_temp[!is.na(all_ages_temp$all_carriers_aoo),])>0){
           test_k_test<-as.numeric(kruskal_test(variant_ages_temp$all_carriers_aoo,all_ages_temp$all_carriers_aoo))
           if (test_k_test <= min_p_value){
             min_p_value<-test_k_test
             min_phenotype<-phenotype
             min_sex<-sex
             min_fam<-fam_hist
           }
         }
       }
     }
   }
   if (min_p_value < age_p_cutoff){
     evidence<-"strong"
     above<-"below"
   }
   else{
     evidence<-"no"
     above<-"above"
   }
   this_list=list(
     "min_p_value"=min_p_value,
     "min_phenotype"=min_phenotype,
     "min_sex"=min_sex,
     "min_fam"=min_fam,
     "above"=above,
     "evidence"=evidence
     )
   return(this_list)
 }

 age_text_line_3=function(variant_variant_browser,all_ages_df,age_p_cutoff){
  if(nrow(all_ages_df[all_ages_df$HGVS==variant_variant_browser & is.na(all_ages_df$HGVS)==F,])==0){
    return(HTML(paste("<b>Note:</b> No ages are recorded for variant ", variant_variant_browser, ". This contributes no evidence of pathogenicity.",sep="")))

  }
  else{
   return(HTML(paste("<b>Note:</b> The most significant p value observed is ",min_age_p_value(variant_variant_browser,all_ages_df,age_p_cutoff)$min_p_value," (Phenotype: ",min_age_p_value(variant_variant_browser,all_ages_df,age_p_cutoff)$min_phenotype," | Sex: ",min_age_p_value(variant_variant_browser,all_ages_df,age_p_cutoff)$min_sex," | Family History : ", min_age_p_value(variant_variant_browser,all_ages_df,age_p_cutoff)$min_fam,"). This is ",min_age_p_value(variant_variant_browser,all_ages_df,age_p_cutoff)$above," the significance threshold of ", formatC(age_p_cutoff, format = "e", digits = 2),". This contributes ",min_age_p_value(variant_variant_browser,all_ages_df,age_p_cutoff)$evidence," evidence of pathogenicity.",sep="")))
 }
}
# AGE FUNCTIONS TO FILTER DATA BASED ON USER SELECTIONS
age_phenotype_selection=function(all_ages,variant_ages_x,gene_ages_x,phenotype_selection){
  this_list=list(
    "all_ages"=all_ages[all_ages$all_carriers_primary_phenotype==phenotype_selection & is.na(all_ages$all_carriers_primary_phenotype) ==FALSE,],
    "variant_ages_x"=variant_ages_x[variant_ages_x$all_carriers_primary_phenotype==phenotype_selection & is.na(variant_ages_x$all_carriers_primary_phenotype) ==FALSE,],
    "gene_ages_x"=gene_ages_x[gene_ages_x$all_carriers_primary_phenotype==phenotype_selection & is.na(gene_ages_x$all_carriers_primary_phenotype) ==FALSE,]
    )
  return(this_list)
}
#
age_phenotype_selection_compare=function(all_ages_1,all_ages_2,variant_ages_x_1,variant_ages_x_2,gene_ages_x_1,gene_ages_x_2,phenotype_selection){
  this_list=list(
    "all_ages_1"=all_ages_1[all_ages_1$all_carriers_primary_phenotype==phenotype_selection & is.na(all_ages_1$all_carriers_primary_phenotype) ==FALSE,],
    "all_ages_2"=all_ages_2[all_ages_2$all_carriers_primary_phenotype==phenotype_selection & is.na(all_ages_2$all_carriers_primary_phenotype) ==FALSE,],
    "variant_ages_x_1"=variant_ages_x_1[variant_ages_x_1$all_carriers_primary_phenotype==phenotype_selection & is.na(variant_ages_x_1$all_carriers_primary_phenotype) ==FALSE,],
    "variant_ages_x_2"=variant_ages_x_2[variant_ages_x_2$all_carriers_primary_phenotype==phenotype_selection & is.na(variant_ages_x_2$all_carriers_primary_phenotype) ==FALSE,],
    "gene_ages_x_1"=gene_ages_x_1[gene_ages_x_1$all_carriers_primary_phenotype==phenotype_selection & is.na(gene_ages_x_1$all_carriers_primary_phenotype) ==FALSE,],
    "gene_ages_x_2"=gene_ages_x_2[gene_ages_x_2$all_carriers_primary_phenotype==phenotype_selection & is.na(gene_ages_x_2$all_carriers_primary_phenotype) ==FALSE,]
    )
  return(this_list)
}
#
age_sex_selection=function(all_ages,variant_ages_x,gene_ages_x,sex_selection){
  this_list=list(
    "all_ages"=all_ages[all_ages$all_carriers_sex==sex_selection & is.na(all_ages$all_carriers_sex) ==FALSE,],
    "variant_ages_x"=variant_ages_x[variant_ages_x$all_carriers_sex==sex_selection & is.na(variant_ages_x$all_carriers_sex) ==FALSE,],
    "gene_ages_x"=gene_ages_x[gene_ages_x$all_carriers_sex==sex_selection & is.na(gene_ages_x$all_carriers_sex) ==FALSE,]
    )
  return(this_list)
}
#
age_sex_selection_compare=function(all_ages_1,all_ages_2,variant_ages_x_1,variant_ages_x_2,gene_ages_x_1,gene_ages_x_2,sex_selection){
  this_list=list(
    "all_ages_1"=all_ages_1[all_ages_1$all_carriers_sex==sex_selection & is.na(all_ages_1$all_carriers_sex) ==FALSE,],
    "all_ages_2"=all_ages_2[all_ages_2$all_carriers_sex==sex_selection & is.na(all_ages_2$all_carriers_sex) ==FALSE,],
    "variant_ages_x_1"=variant_ages_x_1[variant_ages_x_1$all_carriers_sex==sex_selection & is.na(variant_ages_x_1$all_carriers_sex) ==FALSE,],
    "variant_ages_x_2"=variant_ages_x_2[variant_ages_x_2$all_carriers_sex==sex_selection & is.na(variant_ages_x_2$all_carriers_sex) ==FALSE,],
    "gene_ages_x_1"=gene_ages_x_1[gene_ages_x_1$all_carriers_sex==sex_selection & is.na(gene_ages_x_1$all_carriers_sex) ==FALSE,],
    "gene_ages_x_2"=gene_ages_x_2[gene_ages_x_2$all_carriers_sex==sex_selection & is.na(gene_ages_x_2$all_carriers_sex) ==FALSE,]
    )
  return(this_list)
}
#
age_family_selection=function(all_ages,variant_ages_x,gene_ages_x,family_selection){
  this_list=list(
    "all_ages"=all_ages[all_ages$all_carriers_family_history==family_selection & is.na(all_ages$all_carriers_family_history) ==FALSE,],
    "variant_ages_x"=variant_ages_x[variant_ages_x$all_carriers_family_history==family_selection & is.na(variant_ages_x$all_carriers_family_history) ==FALSE,],
    "gene_ages_x"=gene_ages_x[gene_ages_x$all_carriers_family_history==family_selection & is.na(gene_ages_x$all_carriers_family_history) ==FALSE,]
    )
  return(this_list)
}
#
age_family_selection_compare=function(all_ages_1,all_ages_2,variant_ages_x_1,variant_ages_x_2,gene_ages_x_1,gene_ages_x_2,family_selection){
  this_list=list(
    "all_ages_1"=all_ages_1[all_ages_1$all_carriers_family_history==family_selection & is.na(all_ages_1$all_carriers_family_history) ==FALSE,],
    "all_ages_2"=all_ages_2[all_ages_2$all_carriers_family_history==family_selection & is.na(all_ages_2$all_carriers_family_history) ==FALSE,],
    "variant_ages_x_1"=variant_ages_x_1[variant_ages_x_1$all_carriers_family_history==family_selection & is.na(variant_ages_x_1$all_carriers_family_history) ==FALSE,],
    "variant_ages_x_2"=variant_ages_x_2[variant_ages_x_2$all_carriers_family_history==family_selection & is.na(variant_ages_x_2$all_carriers_family_history) ==FALSE,],
    "gene_ages_x_1"=gene_ages_x_1[gene_ages_x_1$all_carriers_family_history==family_selection & is.na(gene_ages_x_1$all_carriers_family_history) ==FALSE,],
    "gene_ages_x_2"=gene_ages_x_2[gene_ages_x_2$all_carriers_family_history==family_selection & is.na(gene_ages_x_2$all_carriers_family_history) ==FALSE,]
    )
  return(this_list)
}
#
# RETURN RESULTS OF FILTER FUNCTIONS
age_data_frame_non_compare_function=function(variant,gene,phenotype_selection,sex_selection,family_history_selection){
  if (sex_selection=="Male"){
    sex_selection<- "M"
  }
  else if(sex_selection=="Female"){
    sex_selection<-"F"
  }
  if (family_history_selection=="Familial"){
    family_history_selection<- "F"
  }
  else if(family_history_selection=="Sporadic"){
    family_history_selection<-"S"
  }
  all_ages<- all_ages_df[all_ages_df$HGVS != variant & is.na(all_ages_df$HGVS)==FALSE,]
  variant_ages_x<- all_ages_df[all_ages_df$HGVS==variant& is.na(all_ages_df$HGVS)==FALSE,]
  gene_ages_x<-all_ages_df[all_ages_df$gene==gene & all_ages_df$HGVS != variant& is.na(all_ages_df$HGVS)==FALSE& is.na(all_ages_df$gene)==FALSE,]
  listicle<- list("all_ages1"=all_ages,"variant_ages_x1"=variant_ages_x,"gene_ages_x1"=gene_ages_x)
  lit_listicle<-c("all_ages","variant_ages_x","gene_ages_x")
  j=0
  for (i in names(listicle)){
    j=j+1
    if (phenotype_selection != "All"){
      listicle[[i]]=age_phenotype_selection(listicle$all_ages1,listicle$variant_ages_x1,listicle$gene_ages_x1,phenotype_selection)[[lit_listicle[j]]]
    }
    if (sex_selection != "All"){ 
      listicle[[i]]=age_sex_selection(listicle$all_ages1,listicle$variant_ages_x1,listicle$gene_ages_x1,sex_selection)[[lit_listicle[j]]]
    }
    if (family_history_selection != "All"){
      listicle[[i]]=age_family_selection(listicle$all_ages1,listicle$variant_ages_x1,listicle$gene_ages_x1,family_history_selection)[[lit_listicle[j]]]    
    }
  }
  this_list <- list("all_ages" = listicle$all_ages1$all_carriers_aoo, "variant_ages_x" = listicle$variant_ages_x1$all_carriers_aoo, "gene_ages_x" = listicle$gene_ages_x1$all_carriers_aoo)
  return(this_list)
}
# MAIN PLOT ECDF
#XXXXXXXXXXXXXXXXXX
render.age_plot <- function(phenotype1,phenotype2,sex1,sex2,fam1,fam2,variant,gene,plot_type,plot_height){
  if (plot_type=="Cumulative Distribution"){
    render.age_plot_ecdf(phenotype1,phenotype2,sex1,sex2,fam1,fam2,variant,gene,plot_height)
  } 
  else {
    render.age_plot_kde(phenotype1,phenotype2,sex1,sex2,fam1,fam2,variant,gene,plot_height)
  }
}
render.age_plot_ecdf <- function(phenotype1,phenotype2,sex1,sex2,fam1,fam2,variant,gene,plot_height){
  if (phenotype1 != "Comparison" & sex1 != "Comparison" & fam1 != "Comparison"){
    non_comparison  <- age_data_frame_non_compare_function(variant,gene,phenotype1,sex1,fam1)
    all_ages    <-non_comparison$all_ages
    variant_ages_x  <-non_comparison$variant_ages_x
    gene_ages_x   <-non_comparison$gene_ages_x
    age_plot_no_compare(
      all_ages,
      gene_ages_x,
      variant_ages_x,
      variant,
      gene,
      phenotype1,
      sex1,
      fam1,
      plot_height
      )
  }
  else if(phenotype1 == "Comparison"){
    phen<- "Comparison"
    sex<-sex2
    fam<-fam2
  }
  else if (sex1 == "Comparison"){
    phen<-phenotype2
    fam<-fam2
    sex<-"Comparison"
  }
  else if (fam1 == "Comparison"){
    sex<-sex2
    phen<-phenotype2
    fam<-"Comparison"
  }
  if(phenotype1 == "Comparison" | sex1 == "Comparison" |fam1 == "Comparison" ){
    comparison      <-age_data_frame_compare_function(variant,gene,phen,sex,fam)
    all_ages_1      <-comparison$all_ages_1
    all_ages_2      <-comparison$all_ages_2
    variant_ages_x_1  <-comparison$variant_ages_x_1
    variant_ages_x_2  <-comparison$variant_ages_x_2
    gene_ages_x_1     <-comparison$gene_ages_x_1
    gene_ages_x_2     <-comparison$gene_ages_x_2
    age_plot_compare(
      all_ages_1,     
      all_ages_2,       
      variant_ages_x_1,   
      variant_ages_x_2, 
      gene_ages_x_1,    
      gene_ages_x_2,        
      variant,
      gene,
      phen,
      sex,
      fam,
      plot_height
      )
  }
}

# MAIN PLOT FOR ECDF NON COMPARE
age_plot_no_compare = function (all_ages_x,gene_ages_x,variant_ages_x,variant,gene,phenotype,sex,family_history,plot_height){
  if (length(variant_ages_x)>0){
    variant_ecdf<-ecdf(variant_ages_x)
    variant_ages_y<- variant_ecdf(variant_ages_x)
  }
  if (length(gene_ages_x)>0){
    gene_ecdf<-ecdf(gene_ages_x)
    gene_ages_y<- gene_ecdf(gene_ages_x)
  }
  if (length(all_ages_x)>0){
    dist_funct<-ecdf(all_ages_x)
  }
  all_ages_y <- dist_funct(all_ages_x)
  # Calculate kruskal test 
  if (length(variant_ages_x >0) & length(all_ages_x >0)){
    kruskal_p <- kruskal_test(variant_ages_x,all_ages_x)
  }
  if(length(all_ages_x)>3){
   smooth_all_ages<-smooth.spline(all_ages_x, all_ages_y, spar=0.35)
   plot_margins_variable
   plot(
    smooth_all_ages,     
    lwd=6,
    col=darkgreen,
    type="l",
    xaxt="n",
    yaxt="n",
    las=1,
    xlab="Age",
    ylab="Cumulative Proportion",
    main="Observed Cumulative Age of Onset",
    cex.lab=1.4*plot_height/600,
    cex.main=1.4*plot_height/600,
    col.main=nicegrey,
    col.lab=nicegrey,
    xlim=c(0,100),
    ylim=c(0,1.3))

    #calculate and plot median 95% CI using bootstraps
  # Plot median and CI for all carriers except P and LP of interest
  if(length(sort(unique(all_ages_x)))>1){
    median_all_ages <- median(all_ages_x)
    all_ages_lower<- median_ci(all_ages_x,"lower")
    all_ages_upper<- median_ci(all_ages_x,"upper")
    lines(x=c(median_all_ages,median_all_ages),y=c(0,1),lty=2,lwd=3*plot_height/600,col=darkgreen)
    rect(all_ages_lower, 0, all_ages_upper, 1,border=NA,col=darkgreen_rgb_fade)
  }
  # Plot median and CI for VUS in gene of interest 
  if(length(sort(unique(gene_ages_x)))>1){
    this_median     <- median(gene_ages_x)
    lines(x=c(this_median,this_median),y=c(0,1),col=nicedarkblue,lty=2,lwd=3*plot_height/600)
    ages_lower<- median_ci(gene_ages_x,"lower")
    ages_upper<- median_ci(gene_ages_x,"upper")
    rect(ages_lower, 0, ages_upper, 1,border=NA,col=nicedarkblue_rgb_fade)
  }
  # Plot median and CI for P and LP variants of interest 
  if(length(sort(unique(variant_ages_x)))>1){
    this_median     <- median(variant_ages_x)
    lines(x=c(this_median,this_median),y=c(0,1),col=nicered,lty=2,lwd=3*plot_height/600)
    ages_lower<- median_ci(variant_ages_x,"lower")
    ages_upper<- median_ci(variant_ages_x,"upper")
    rect(ages_lower, 0, ages_upper, 1,border=NA,col=nicered_rgb_fade)
  }
  # Replot all ages line over CIs
  lines(smooth_all_ages,lwd=6*plot_height/600,col=darkgreen)
  # Plot all age line in legend 
  lines(x=c(-1,0.5),y=c(1.3*0.85,1.3*0.85),lwd=6*plot_height/600,col=darkgreen)
  # Plot all age text in legend 
  text(5,(1.3*0.85),"All Cases Distribution",col=nicegrey,adj=0,cex=1*plot_height/600)
  # Plot kruskal test 
  if (length(variant_ages_x >0)){
    points(46.5,(1.3*0.85),cex=2*plot_height/600,pch=21,col="white",bg=darkgreen)
    points(45,(1.3*0.85),cex=2*plot_height/600,pch=21,col="white",bg=nicered)
    text(50,(1.3*0.85),paste("Test For Difference: p=",kruskal_p,sep=""),col=nicegrey,adj=0,cex=1*plot_height/600)
  }
}
  #create a blank plot if not enough data to calculate spline (likely to only be needed when testing a data subset)
  else if (length(all_ages_x)<=3){  
    plot(50,0.5,
      col="white",       
      pch=21,
      bg="white",      
      type="l",
      xaxt="n",
      yaxt="n",
      las=1,
      xlab="Age",
      ylab="Cumulative Proportion",
      main="Observed Cumulative Age of Onset",
      cex.lab=1.4*plot_height/600,
      cex.main=1.4*plot_height/600,
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,100),
      ylim=c(0,1.3))
    text(40,0.5,"Insufficent Data",col=nicegrey,adj=0,cex=plot_height/600)
  }
  # plot gene and variant points if there are any
  if (length(gene_ages_x) >0){
    points(
      x=gene_ages_x,
      y=gene_ages_y,
      cex=1.2*plot_height/600,
      pch=21,
      col="white",
      bg=nicedarkblue
      )
    points(x=0,y=(1.3*0.9),cex=1.1*plot_height/600,pch=21,col="white",bg=nicedarkblue)
    text(5,(1.3*0.9),gene,col=nicegrey,adj=0,cex=plot_height/600)
  }
  else{
    text(5,(1.3*0.9),paste("No observed ages for",gene," excluding ",variant," carriers",sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  }
  if (length(variant_ages_x)>0){
    points(
      x=variant_ages_x,
      y=variant_ages_y,
      cex=1.6*plot_height/600,
      pch=21,
      col="white",
      bg=nicered
      )
    points(x=0,y=(1.3*0.95),cex=1.6*plot_height/600,pch=21,col="white",bg=nicered)
    text(5,(1.3*0.95),variant,col=nicegrey,adj=0,cex=plot_height/600)
  }
  else {
    text(5,(1.3*0.95),paste("No observed ages for ",variant," carriers",sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  }
  selections_legend(phenotype,sex,family_history,1.3,plot_height)
  axis(side=1,lwd=1,cex.axis=1*plot_height/600,col.axis=nicegrey,col=nicegrey,at = seq(0, 100, by = 20))
  axis(side=2,lwd=1,las=1,cex.axis=1*plot_height/600,col.axis=nicegrey,col=nicegrey,at = seq(0, 1, by = 0.2))
  box(lwd=2*plot_height/600,col=nicegrey)
}
#

# MAIN PLOT FOR ECDF COMPARE 
age_data_frame_compare_function=function(variant,gene,phenotype,sex,family_history){
  if (sex=="Male"){
    sex<- "M"
  }
  else if(sex=="Female"){
    sex<-"F"
  }
  if (family_history=="Familial"){
    family_history<- "F"
  }
  else if(family_history=="Sporadic"){
    family_history<-"S"
  }
  all_ages_1<- all_ages_df[all_ages_df$HGVS != variant & is.na(all_ages_df$HGVS) == FALSE,]
  all_ages_2<- all_ages_df[all_ages_df$HGVS != variant & is.na(all_ages_df$HGVS) == FALSE,]
  variant_ages_x_1<- all_ages_df[all_ages_df$HGVS==variant & is.na(all_ages_df$HGVS) == FALSE,]
  variant_ages_x_2<- all_ages_df[all_ages_df$HGVS==variant & is.na(all_ages_df$HGVS) == FALSE,]
  gene_ages_x_1<-all_ages_df[all_ages_df$gene==gene & all_ages_df$HGVS != variant& is.na(all_ages_df$gene) == FALSE,]
  gene_ages_x_2<-all_ages_df[all_ages_df$gene==gene & all_ages_df$HGVS != variant& is.na(all_ages_df$gene) == FALSE,]
  listicle<- list("all_ages1"=all_ages_1,"all_ages2"=all_ages_2,"variant_ages_x1"=variant_ages_x_1,"variant_ages_x2"=variant_ages_x_2,"gene_ages_x1"=gene_ages_x_1,"gene_ages_x2"=gene_ages_x_2)
  lit_listicle<-c("all_ages_1","all_ages_2","variant_ages_x_1","variant_ages_x_2","gene_ages_x_1","gene_ages_x_2")
  if (phenotype=="Comparison"){
    # if phenotype is comparison set 1= ALS, 2 = FTD, after set sex and fam
    listicle$all_ages1<-        all_ages_1[all_ages_1$all_carriers_primary_phenotype=="ALS"&is.na(all_ages_1$all_carriers_primary_phenotype) ==FALSE,]
    listicle$all_ages2<-        all_ages_2[all_ages_2$all_carriers_primary_phenotype=="FTD"&        is.na(all_ages_2$all_carriers_primary_phenotype) ==FALSE,]
    listicle$variant_ages_x1<-  variant_ages_x_1[variant_ages_x_1$all_carriers_primary_phenotype=="ALS"&  is.na(variant_ages_x_1$all_carriers_primary_phenotype) ==FALSE,]
    listicle$variant_ages_x2<-  variant_ages_x_2[variant_ages_x_2$all_carriers_primary_phenotype=="FTD"&  is.na(variant_ages_x_2$all_carriers_primary_phenotype) ==FALSE,]
    listicle$gene_ages_x1<-     gene_ages_x_1[gene_ages_x_1$all_carriers_primary_phenotype=="ALS"&     is.na(gene_ages_x_1$all_carriers_primary_phenotype) ==FALSE,]
    listicle$gene_ages_x2<-     gene_ages_x_2[gene_ages_x_2$all_carriers_primary_phenotype=="FTD"&     is.na(gene_ages_x_2$all_carriers_primary_phenotype) ==FALSE,] 
    j=0
    for (i in names(listicle)){
      j=j+1
      if (sex != "All"){ 
        listicle[[i]]=age_sex_selection_compare(listicle$all_ages1,listicle$all_ages2,listicle$variant_ages_x1,listicle$variant_ages_x2,listicle$gene_ages_x1,listicle$gene_ages_x2,sex)[[lit_listicle[j]]]
      }
      if (family_history != "All"){
        listicle[[i]]=age_family_selection_compare(listicle$all_ages1,listicle$all_ages2,listicle$variant_ages_x1,listicle$variant_ages_x2,listicle$gene_ages_x1,listicle$gene_ages_x2,family_history)[[lit_listicle[j]]]
      }
    }
  }
  else if (sex=="Comparison"){
    listicle$all_ages1<-all_ages_1[all_ages_1$all_carriers_sex=="M" & is.na(all_ages_1$all_carriers_sex)==FALSE,]
    listicle$all_ages2<-all_ages_2[all_ages_2$all_carriers_sex=="F" & is.na(all_ages_2$all_carriers_sex)==FALSE,]
    listicle$variant_ages_x1<- variant_ages_x_1[variant_ages_x_1$all_carriers_sex=="M" & is.na(variant_ages_x_1$all_carriers_sex) ==FALSE,]
    listicle$variant_ages_x2<- variant_ages_x_2[variant_ages_x_2$all_carriers_sex=="F" & is.na(variant_ages_x_2$all_carriers_sex) ==FALSE,]
    listicle$gene_ages_x1<-gene_ages_x_1[gene_ages_x_1$all_carriers_sex=="M" & is.na(gene_ages_x_1$all_carriers_sex) ==FALSE,]
    listicle$gene_ages_x2<-gene_ages_x_2[gene_ages_x_2$all_carriers_sex=="F" & is.na(gene_ages_x_2$all_carriers_sex) ==FALSE,]
    j=0
    for (i in names(listicle)){
      j=j+1
      if (phenotype!="All"){
        listicle[[i]]=age_phenotype_selection_compare(listicle$all_ages1,listicle$all_ages2,listicle$variant_ages_x1,listicle$variant_ages_x2,listicle$gene_ages_x1,listicle$gene_ages_x2,phenotype)[[lit_listicle[j]]]
      }
      if (family_history != "All"){
        listicle[[i]]=age_family_selection_compare(listicle$all_ages1,listicle$all_ages2,listicle$variant_ages_x1,listicle$variant_ages_x2,listicle$gene_ages_x1,listicle$gene_ages_x2,family_history)[[lit_listicle[j]]]
      }
    }
  }
  else if (family_history=="Comparison"){
    listicle$all_ages1<-all_ages_1[all_ages_1$all_carriers_family_history=="F" & is.na(all_ages_1$all_carriers_family_history)==FALSE,]
    listicle$all_ages2<-all_ages_2[all_ages_2$all_carriers_family_history=="S" & is.na(all_ages_2$all_carriers_family_history)==FALSE,]
    listicle$variant_ages_x1<- variant_ages_x_1[variant_ages_x_1$all_carriers_family_history=="F" & is.na(variant_ages_x_1$all_carriers_family_history) ==FALSE,]
    listicle$variant_ages_x2<- variant_ages_x_2[variant_ages_x_2$all_carriers_family_history=="S" & is.na(variant_ages_x_2$all_carriers_family_history) ==FALSE,]
    listicle$ gene_ages_x1<-gene_ages_x_1[gene_ages_x_1$all_carriers_family_history=="F" & is.na(gene_ages_x_1$all_carriers_family_history) ==FALSE,]
    listicle$ gene_ages_x2<-gene_ages_x_2[gene_ages_x_2$all_carriers_family_history=="S" & is.na(gene_ages_x_2$all_carriers_family_history) ==FALSE,]
    j=0
    for (i in names(listicle)){
      j=j+1
      if (phenotype!="All"){
        listicle[[i]]=age_phenotype_selection_compare(listicle$all_ages1,listicle$all_ages2,listicle$variant_ages_x1,listicle$variant_ages_x2,listicle$gene_ages_x1,listicle$gene_ages_x2,phenotype)[[lit_listicle[j]]]
      }
      if (sex != "All"){ 
        listicle[[i]]=age_sex_selection_compare(listicle$all_ages1,listicle$all_ages2,listicle$variant_ages_x1,listicle$variant_ages_x2,listicle$gene_ages_x1,listicle$gene_ages_x2,sex)[[lit_listicle[j]]]
      }
    }
  }
  this_list<- list("all_ages_1"= listicle$all_ages1$all_carriers_aoo,"all_ages_2"= listicle$all_ages2$all_carriers_aoo,"variant_ages_x_1"= listicle$variant_ages_x1$all_carriers_aoo,"variant_ages_x_2"= listicle$variant_ages_x2$all_carriers_aoo,"gene_ages_x_1"= listicle$gene_ages_x1$all_carriers_aoo,"gene_ages_x_2"= listicle$gene_ages_x2$all_carriers_aoo)
  return(this_list)
}
#
# # CALCULATE ECDF
# age_distribution=function(age,classifier,variant){
#   if (classifier == "variant"){
#     temp_df <- lit_review_pen[lit_review_pen$HGVS==input$variant_variant_browser,]
#   } 
#   else if (classifier == "gene"){
#     temp_df <- lit_review_pen[lit_review_pen$gene==input$variant_variant_browser,]
#   }
#   ages<- as.numeric(na.omit(separate_rows(temp_df,all_carriers_aoo,sep="\\|")$all_carriers_aoo))
#   ages<- sort(ages[complete.cases(ages)])
#   dist<-ecdf(ages)
#   return(round(100*dist(age),0))
# }
#
# MAIN PLOT KDE
render.age_plot_kde = function(phenotype_selection,sex_selection,family_history_selection,phenotype_selection2,sex_selection2,family_history2_selection,input_variant,input_gene,plot_height){
  if (sex_selection=="Male"){
    sex_selection<- "M"
  }
  else if(sex_selection=="Female"){
    sex_selection<-"F"
  }
  if (family_history_selection=="Familial"){
    family_history_selection<- "F"
  }
  else if(family_history_selection=="Sporadic"){
    family_history_selection<-"S"
  }
  all_ages<-all_ages_df
  if (phenotype_selection != "Comparison" & sex_selection != "Comparison" & family_history_selection != "Comparison"){
    variant_ages_x<-all_ages_df[all_ages_df$HGVS == input_variant & is.na(all_ages_df$HGVS)==FALSE,]
    gene_ages_x<-all_ages_df[all_ages_df$gene==input_gene & all_ages_df$HGVS != input_variant & is.na(all_ages_df$gene)==FALSE & is.na(all_ages_df$HGVS)==FALSE,]
    listicle<- list("all_ages1"=all_ages,"variant_ages_x1"=variant_ages_x,"gene_ages_x1"=gene_ages_x)
    lit_listicle<-c("all_ages","variant_ages_x","gene_ages_x")
    j=0
    for (i in names(listicle)){
      j=j+1
      if (phenotype_selection != "All"){
        listicle[[i]]=age_phenotype_selection(listicle$all_ages1,listicle$variant_ages_x1,listicle$gene_ages_x1,phenotype_selection)[[lit_listicle[j]]]
      }
      if (sex_selection != "All"){ 
        listicle[[i]]=age_sex_selection(listicle$all_ages1,listicle$variant_ages_x1,listicle$gene_ages_x1,sex_selection)[[lit_listicle[j]]]
      }
      if (family_history_selection != "All"){
        listicle[[i]]=age_family_selection(listicle$all_ages1,listicle$variant_ages_x1,listicle$gene_ages_x1,family_history_selection)[[lit_listicle[j]]]    
      }
    }
    kde_plot(listicle$variant_ages_x1$all_carriers_aoo,listicle$gene_ages_x1$all_carriers_aoo,listicle$all_ages1$all_carriers_aoo[listicle$all_ages1$HGVS != input_variant & is.na(listicle$all_ages1$HGVS)==F],input_variant,input_gene,phenotype_selection,sex_selection,family_history_selection,plot_height)

  }
  else if(phenotype_selection == "Comparison"){
    phen<- "Comparison"
    sex<-sex_selection2
    fam<-family_history2_selection
  }
  else if (sex_selection == "Comparison"){
    phen<-phenotype_selection2
    fam<-family_history2_selection
    sex<-"Comparison"
  }
  else if (family_history_selection == "Comparison"){
    sex<-sex_selection2
    phen<-phenotype_selection2
    fam<-"Comparison"
  }
  if(phenotype_selection == "Comparison" | sex_selection == "Comparison" |family_history_selection == "Comparison" ){
    comparison      <-age_data_frame_compare_function(input_variant,input_gene,phen,sex,fam)
    all_ages_1      <-comparison$all_ages_1
    all_ages_2      <-comparison$all_ages_2
    variant_ages_x_1  <-comparison$variant_ages_x_1
    variant_ages_x_2  <-comparison$variant_ages_x_2
    gene_ages_x_1     <-comparison$gene_ages_x_1
    gene_ages_x_2     <-comparison$gene_ages_x_2
    kde_plot_compare(all_ages_1,all_ages_2,variant_ages_x_1,variant_ages_x_2,gene_ages_x_1,gene_ages_x_2,phen,sex,fam,input_variant,input_gene,plot_height)
  }
}




# KDE PLOT NO COMPARE
kde_plot=function(variant_ages_x,gene_ages_x,all_ages_x,variant,gene,phenotype,sex,family_history,plot_height){
  ### calculate kde ###
  variant_kde<-dense_function(variant_ages_x)
  gene_kde<-dense_function(gene_ages_x)
  all_ages_kde<-dense_function(all_ages_x)
  ### set plot height, minimum=0.004 
  select_max      <- max(variant_kde$y,gene_kde$y,all_ages_kde$y,0.004)
  y_max<- 1.3*select_max

  ### select medians ###
  if(length(all_ages_x)>0){
    median_all_ages <- median(all_ages_x)
  }
  ### calculate kruskal wallis p value ###
  if (length(variant_ages_x >0) & length(all_ages_x >0)){
    kruskal_p <- kruskal_test(variant_ages_x,all_ages_x)
  }

  #create plot
  if (length(all_ages_x)>0){
    plot_margins_variable
    plot(
      x=c(median_all_ages,median_all_ages),
      y=c(0,select_max),
      lty=2,
      lwd=3,
      col=darkgreen,
      type="l",
      las=1,
      xlab="Age",
      ylab="Density",
      main="Observed Density Age of Onset",
      cex.lab=1.4*plot_height/600,
      cex.main=1.4*plot_height/600,
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,100),
      ylim=c(0,y_max),
      xaxt="n",
      yaxt="n")
    #calculate and plot median 95% CI using bootstraps
    # all ages
    if(length(sort(unique(all_ages_x)))>1){
      all_ages_lower<- median_ci(all_ages_x,"lower")
      all_ages_upper<- median_ci(all_ages_x,"upper")
      rect(all_ages_lower, 0, all_ages_upper, select_max,border=NA,col=darkgreen_rgb_fade)
    }
    # gene
    kde_median(gene_ages_x,select_max,nicedarkblue,nicedarkblue_rgb_fade,plot_height)
    # variant
    kde_median(variant_ages_x,select_max,nicered,nicered_rgb_fade,plot_height)
    # all ages
    lines(all_ages_kde,col=darkgreen,lwd=6*plot_height/600)
    # gene
    if (length(gene_ages_x) > 1){
      lines(gene_kde,col=nicedarkblue,lwd=6*plot_height/600)
    }
    # variant
    if (length(variant_ages_x) > 1){
      lines(variant_kde,col=nicered,lwd=6*plot_height/600)
    }
    selections_legend(phenotype,sex,family_history,y_max,plot_height)
    axis(side=1,lwd=1*plot_height/600,cex.axis=1*plot_height/600,col.axis=nicegrey,col=nicegrey,at = seq(0,100,by=20))
    axis(side=2,lwd=1*plot_height/600,las=1,cex.axis=1*plot_height/600,col.axis=nicegrey,col=nicegrey,at = seq(0, select_max,by=(round(select_max/4,3))))
    box(lwd=2*plot_height/600,col=nicegrey)
    # Create Legend
    # Variant
    var_height<- 0.95
    kde_legend_1(variant_ages_x,y_max,var_height,variant,nicered,nicered_rgb_fade,plot_height)
    gene_height<-0.9
    kde_legend_1(gene_ages_x,y_max,gene_height,gene,nicedarkblue,nicedarkblue_rgb_fade,plot_height)
    # All individuals
    all_height<-0.85
    kde_legend(y_max,all_height,"All cases",darkgreen,darkgreen_rgb_fade,plot_height)
    # kruskal
    if (length(variant_ages_x >0)){
      points(58.5,y_max*all_height,cex=2*plot_height/600,pch=21,col="white",bg=darkgreen)
      points(57,y_max*all_height,cex=2*plot_height/600,pch=21,col="white",bg=nicered)
      text(60,y_max*all_height,paste("Test For Difference: p=",kruskal_p,sep=""),col=nicegrey,adj=0,cex=plot_height/600)
    }
  }
  else if (length(all_ages_x)==0){ 
    plot_margins_variable
    plot(50,0.5,
      col="white",       
      pch=21,
      bg="white",      
      type="l",
      xaxt="n",
      yaxt="n",
      las=1,
      xlab="Age",
      ylab="Density",
      main="Observed Density Age of Onset",
      cex.lab=1.4*plot_height/600,
      cex.main=1.4*plot_height/600,
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,100),
      ylim=c(0,1.3))
    text(40,0.5,"Insufficent Data",col=nicegrey,adj=0,cex=plot_height/600)
  }
}
#
# KDE COMPARE PLOT
kde_plot_compare=function(all_ages_1,all_ages_2,variant_ages_x_1,variant_ages_x_2,gene_ages_x_1,gene_ages_x_2,phenotype,sex,family_history,variant,gene,plot_height){
  if (phenotype == "Comparison"){
    phenotype="ALS|FTD"
  }
  if (sex== "Comparison"){
    sex="Male|Female"
  } 
  if (family_history == "Comparison"){
    family_history="Familial|Sporadic"
  }
    # create density functions
    all_ages_1_kde<-dense_function(all_ages_1)
    all_ages_2_kde<-dense_function(all_ages_2)
    variant_ages_x_1_kde<-dense_function(variant_ages_x_1)
    variant_ages_x_2_kde<-dense_function(variant_ages_x_2)
    gene_ages_x_1_kde<-dense_function(gene_ages_x_1)
    gene_ages_x_2_kde<-dense_function(gene_ages_x_2)
    ### set plot height, minimum=0.004 
    select_max <- max(all_ages_1_kde$y,all_ages_2_kde$y,variant_ages_x_1_kde$y,variant_ages_x_2_kde$y,gene_ages_x_1_kde$y,gene_ages_x_2_kde$y,0.004)
    y_max<- 1.3*select_max
    # create blank plot 
    plot_margins_variable
    plot(
      x=40,
      y=select_max/2,
      col="white",       
      pch=21,
      bg="white",      
      type="l",
      xaxt="n",
      yaxt="n",
      las=1,
      xlab="Age",
      ylab="Density",
      main="Observed Density Age of Onset",
      cex.lab=1.4*plot_height/600,
      cex.main=1.4*plot_height/600,
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,100),
      ylim=c(0,y_max)
      )
    # plot density line (>1) or single line if single point (==1)
    density_line_plot(all_ages_1,all_ages_1_kde,darkgreen,select_max,plot_height)
    density_line_plot(all_ages_2,all_ages_2_kde,lightgreen,select_max,plot_height)
    density_line_plot(variant_ages_x_1,variant_ages_x_1_kde,nicered,select_max,plot_height)
    density_line_plot(variant_ages_x_2,variant_ages_x_2_kde,niceorange,select_max,plot_height)
    density_line_plot(gene_ages_x_1,gene_ages_x_1_kde,nicedarkblue,select_max,plot_height)
    density_line_plot(gene_ages_x_2,gene_ages_x_2_kde,nicelightblue,select_max,plot_height)
    # plot legend
    legend_line_line(y_max,darkgreen,lightgreen,"All Cases Distribution",0.85,plot_height)
    legend_line_line(y_max,nicedarkblue,nicelightblue,gene,0.9,plot_height)
    legend_line_line(y_max,nicered,niceorange,variant,0.95,plot_height)
    # tidy plot 
    selections_legend(phenotype,sex,family_history,y_max,plot_height)
    axis(side=1,lwd=1*plot_height/600,cex.axis=1,col.axis=nicegrey,col=nicegrey,at = seq(0,100,by=20))
    axis(side=2,lwd=1*plot_height/600,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at = seq(0, select_max,by=(round(select_max/4,3))))
    box(lwd=2*plot_height/600,col=nicegrey)
  }
  #
  density_line_plot=function(ages,ages_kde,colour,select_max,plot_height){
    if (length(ages) > 1){
      lines(ages_kde,col=colour,lwd=6*plot_height/600)
    }
    else if (length(ages)==1){
      lines(x=c(ages,ages),y=c(0,select_max),col=colour,lty=2,lwd=3*plot_height/600)
    }
  }
  age_plot_compare = function (     
    all_ages_1,       
    all_ages_2,     
    variant_ages_x_1,
    variant_ages_x_2, 
    gene_ages_x_1,    
    gene_ages_x_2,        
    variant,
    gene,
    phenotype,
    sex,
    family_history,
    plot_height){
    if (phenotype == "Comparison"){
      phenotype="ALS|FTD"
    }
    if (sex== "Comparison"){
      sex="Male|Female"
    } 
    if (family_history == "Comparison"){
      family_history="Familial|Sporadic"
    }
    if (length(variant_ages_x_1)>0){
      variant1_ecdf<-ecdf(variant_ages_x_1)
      variant1_ages_y<- variant1_ecdf(variant_ages_x_1)
    }
    if (length(variant_ages_x_2)>0){
      variant2_ecdf<-ecdf(variant_ages_x_2)
      variant2_ages_y<- variant2_ecdf(variant_ages_x_2)
    }
    if (length(gene_ages_x_1)>0){
      gene1_ecdf<-ecdf(gene_ages_x_1)
      gene1_ages_y<- gene1_ecdf(gene_ages_x_1)
    }
    if (length(gene_ages_x_2)>0){
      gene2_ecdf<-ecdf(gene_ages_x_2)
      gene2_ages_y<- gene2_ecdf(gene_ages_x_2)
    }
    if (length(all_ages_1)>0){
      dist_funct1<-ecdf(all_ages_1)
      all_ages_y_1 <- dist_funct1(all_ages_1)
    }
    if (length(all_ages_2)>0){
      dist_funct2<-ecdf(all_ages_2)
      all_ages_y_2 <- dist_funct2(all_ages_2)
    }
    if(length(all_ages_1)>3){
      smooth_all_ages2<-smooth.spline(all_ages_2, all_ages_y_2, spar=0.35)
      plot(
        smooth_all_ages2,     
        lwd=6*plot_height/600,
        col=lightgreen,
        type="l",
        xaxt="n",
        yaxt="n",
        las=1,
        xlab="Age",
        ylab="Cumulative Proportion",
        main="Observed Cumulative Age of Onset",
        cex.lab=1.4*plot_height/600,
        cex.main=1.4*plot_height/600,
        col.main=nicegrey,
        col.lab=nicegrey,
        xlim=c(0,100),
        ylim=c(0,1.3))
      legend_line_line(1.3,darkgreen,lightgreen,"All Cases Distribution",0.85,plot_height)
    }
  #create a blank plot if not enough data to calculate spline (likely to only be needed when testing a data subset)
  else if (length(all_ages_1)<=3){  
    plot(50,0.5,
      col="white",       
      pch=21,
      bg="white",      
      type="l",
      xaxt="n",
      yaxt="n",
      las=1,
      xlab="Age",
      ylab="Cumulative Proportion",
      main="Observed Cumulative Age of Onset",
      cex.lab=1.4*plot_height/600,
      cex.main=1.4*plot_height/600,
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,100),
      ylim=c(0,1.3))
    text(40,0.5,"Insufficent Data",col=nicegrey,adj=0,cex=plot_height/600)
  }
  if(length(all_ages_2)>3){
    smooth_all_ages1<-smooth.spline(all_ages_1, all_ages_y_1, spar=0.35)
    lines(smooth_all_ages1,lwd=6*plot_height/600,col=darkgreen)
  }
  if (length(gene_ages_x_1) >0){
    points(
      x=gene_ages_x_1,
      y=gene1_ages_y,
      cex=1.2*plot_height/600,
      pch=21,
      col="white",
      bg=nicedarkblue
      )
  }
  if (length(gene_ages_x_2) >0){
    points(
      x=gene_ages_x_2,
      y=gene2_ages_y,
      cex=1.2*plot_height/600,
      pch=21,
      col="white",
      bg=nicelightblue
      )
  }
  if (length(gene_ages_x_1)>0 | length(gene_ages_x_2)>0){
    points(x=0,y=(1.3*0.9),cex=1.2*plot_height/600,pch=21,col="white",bg=nicedarkblue)
    points(x=3,y=(1.3*0.9),cex=1.2*plot_height/600,pch=21,col="white",bg=nicelightblue)
    text(1.3,(1.3*0.9),"|",col=nicegrey,adj=0,cex=plot_height/600)
    text(5,(1.3*0.9),gene,col=nicegrey,adj=0,cex=plot_height/600)
  }
  else {
    text(5,(1.3*0.9),paste("No observed ages for ",gene," carriers",sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  }
  if (length(variant_ages_x_1)>0){
    points(
      x=variant_ages_x_1,
      y=variant1_ages_y,
      cex=1.6*plot_height/600,
      pch=21,
      col="white",
      bg=nicered
      )
  }
  if (length(variant_ages_x_2)>0){
    points(
      x=variant_ages_x_2,
      y=variant2_ages_y,
      cex=1.6*plot_height/600,
      pch=21,
      col="white",
      bg=niceorange
      )
  }
  if (length(variant_ages_x_2)>0 | length(variant_ages_x_1)>0){
    points(x=0,y=(1.3*0.95),cex=1.6*plot_height/600,pch=21,col="white",bg=nicered)
    points(x=3,y=(1.3*0.95),cex=1.6*plot_height/600,pch=21,col="white",bg=niceorange)
    text(1.3,(1.3*0.95),"|",col=nicegrey,adj=0,cex=plot_height/600)
    text(5,(1.3*0.95),variant,col=nicegrey,adj=0,cex=plot_height/600)
  }
  else {
    text(5,(1.3*0.95),paste("No observed ages for ",variant," carriers",sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  }
  selections_legend(phenotype,sex,family_history,1.3,plot_height)
  axis(side=1,lwd=1,cex.axis=1*plot_height/600,col.axis=nicegrey,col=nicegrey,at = seq(0, 100, by = 20))
  axis(side=2,lwd=1,las=1,cex.axis=1*plot_height/600,col.axis=nicegrey,col=nicegrey,at = seq(0, 1, by = 0.2))
  box(lwd=2*plot_height/600,col=nicegrey)
}
#
# Calculate medians for kde
kde_median=function(ages,select_max,colour,col_fade,plot_height){
  if (length(ages) > 0){
    this_median     <- median(ages)
    lines(x=c(this_median,this_median),y=c(0,select_max),col=colour,lty=2,lwd=3*plot_height/600)
  }
  if(length(sort(unique(ages)))>1){
    ages_lower<- median_ci(ages,"lower")
    ages_upper<- median_ci(ages,"upper")
    rect(ages_lower, 0, ages_upper, select_max,border=NA,col=col_fade)
  }
}
#
# Calculate kde
dense_function=function(ages){
  if (length(ages) > 1){
    this     <-density(ages)
  }
  else {
    this     <-list(y=0)
  }
  return(this)
}
#
# LEGEND FUNCTIONS
kde_legend=function(y_max,height,variable,colour,rect_colour,plot_height){
  lines(x=c(0,0),y=c(y_max*(height-0.02),y_max*(height+0.02)),lwd=6*plot_height/600,col=colour)
  rect(1.5,y_max*(height-0.02),6,y_max*(height+0.02),col=rect_colour,border=NA)
  text(0.6,y_max*height,"|",col=nicegrey,adj=0,cex=plot_height/600)
  lines(x=c(1.5,6),y=c(y_max*height,y_max*height),col=colour,lty=2,lwd=3*plot_height/600)    
  text(7,y_max*height,paste(variable," | median and 95% CI",sep=""),col=nicegrey,adj=0,cex=plot_height/600)
}
#
legend_line_line=function(y_max,colour1,colour2,text,height,plot_height){
  lines(x=c(-1,0.5),y=c(y_max*height,y_max*height),lwd=6*plot_height/600,col=colour1)
  lines(x=c(3,4.5),y=c(y_max*height,y_max*height),lwd=6*plot_height/600,col=colour2)
  text(1.3,(y_max*height),"|",col=nicegrey,adj=0,cex=plot_height/600)
  text(5,(y_max*height),text,col=nicegrey,adj=0,cex=plot_height/600)
}
#
selections_legend=function(phenotype,sex,family_history,height,plot_height){
  if (sex=="F"){
    sex <- "Female"
  }
  if (sex=="M"){
    sex <- "Male"
  }
  if (family_history=="F"){
    family_history <- "Familial"
  }
  if (family_history=="S"){
    family_history <- "Sporadic"
  }
  text(0,height,paste("Phenotype: ",phenotype,sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  text(30,height,paste("Sex: ",sex,sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  text(60,height,paste("Family History: ",family_history,sep=""),col=nicegrey,adj=0,cex=plot_height/600)
}
#
kde_legend_1=function(ages_x,y_max,var_height,variable,col,col_fade,plot_height){
  if(length(ages_x)>1){
    kde_legend(y_max,var_height,variable,col,col_fade,plot_height)
  }
  else if (length(ages_x)>1){
    lines(x=c(1.5,6),y=c(y_max*var_height,y_max*var_height),col=col,lty=2,lwd=3*plot_height/600)
    text(7,y_max*var_height,paste("Single observed ",variable," carrier",sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  }
  else{
    text(1,y_max*var_height,paste("No observed ages for ",variable," carriers",sep=""),col=nicegrey,adj=0,cex=plot_height/600)
  }
}
#
# MISC
create_dataframe_of_individuals <- function(lit_review_pen){
  # This creates a dataframe of individuals and their phenotypes. Used for tables on variant browser / phenotype and region browser / individuals
  # Pull out relevant columns and filter to carriers in literature
  #pheno_frame <- subset(lit_review_pen, select=c("HGVS","HGVS2","identifier","all_carriers_nationality","all_carriers_ethnicity","all_carriers_ethnicity_explanation","all_carriers_primary_phenotype","all_carriers_detailed_phenotype","all_carriers_family_history","all_carriers_sex","all_carriers_aoo","all_carriers_disease_duration","all_carriers_onset_site","all_carriers_zygosity","all_carriers_cognitive_impairment","all_carriers_denovo","all_carriers_denovo_confirmed","all_carriers_concurrent_variants","all_carriers_pmid","all_carriers_comment","all_carriers_pedigree","all_carriers_continent"),(lit_review_pen$all_carriers_count>0 & !is.na(lit_review_pen$all_carriers_count)))
  pheno_frame <- subset(lit_review_pen, select=c("HGVS","identifier","all_carriers_nationality","all_carriers_ethnicity","all_carriers_ethnicity_explanation","all_carriers_primary_phenotype","all_carriers_detailed_phenotype","all_carriers_family_history","all_carriers_sex","all_carriers_aoo","all_carriers_disease_duration","all_carriers_onset_site","all_carriers_zygosity","all_carriers_cognitive_impairment","all_carriers_denovo","all_carriers_denovo_confirmed","all_carriers_concurrent_variants","all_carriers_pmid","all_carriers_comment","all_carriers_pedigree","all_carriers_continent"),(lit_review_pen$all_carriers_count>0 & !is.na(lit_review_pen$all_carriers_count)))  
  # Split based on pipe
  pheno_frame <- data.frame(separate_rows(pheno_frame,all_carriers_nationality: all_carriers_continent,sep="\\|"))
  #pheno_frame_colnames<c("HGVS","identifier","Nationality","Ethnicity","Ethnicity_Explanation","Primary_Phenotype","Detailed_Phenotype","Family_History",      "AOO","Duration","Onset","Zygosity","Cognitive_Impairment","De_novo","De_novo_parentage","Concurrent_Variants","PMID","Comments","Pedigree","Continent")
  # Split concurrent variants onto multiple lines for display
  pheno_frame$all_carriers_concurrent_variants<-gsub(";", "\n",pheno_frame$all_carriers_concurrent_variants)
  # Format comments 
  pheno_frame$all_carriers_comment <- gsub("([A-z]\\. )","\\1\n",pheno_frame$all_carriers_comment)
  # Rename
  #names(pheno_frame) <- c("HGVS_use",pheno_frame_colnames)
  names(pheno_frame) <- pheno_frame_colnames

  return(pheno_frame)
}

all_countries<-function(lit_review_pen){
  # Function to find list of countries for drop down menu on region browser 
  # This is a manually curated list, some individuals have these non countries recorded as nationality. We want to exclude these from drop down list
  omit_countries<-c("Asia(unspecified)","Europe(unspecified)", "French_Canada", "Latin_America", "NA", "North_America",  "Pacific_Islands", "The_Balkans", "Korea")
  # Split based on pipe
  countries<-separate_rows(lit_review_pen[lit_review_pen$all_carriers_count>0 & !is.na(lit_review_pen$all_carriers_count),],all_carriers_nationality,sep="\\|")
  # Split based on slash
  countries<-separate_rows(countries,all_carriers_nationality,sep="\\/")$all_carriers_nationality
  # Remove if in omit list
  countries<-unique(countries[!(countries %in%  omit_countries) & !is.na(countries)])
  return(sort(countries))
}
median_ci=function(data,range){
  tempdf<-data.frame(data=data,group=1)
  tempdf_median<- rcompanion::groupwiseMedian(data ~ group,data=tempdf,conf = 0.95,R= 100,percentile = TRUE,bca = FALSE,digits = 30)
  if(range=="lower"){
    return(tempdf_median$Percentile.lower)
  } 
  else if (range == "upper"){

    return(tempdf_median$Percentile.upper)
  }
}
#
kruskal_test=function(data1,data2){
  tempdf<-data.frame(data=data1,group=1)
  tempdf<-rbind(tempdf,data.frame(data=data2,group=2))
  kruskal<-kruskal.test(data ~ group, data = tempdf)$p.value
  if (kruskal >= 0.001){
    return(round(kruskal,4))
  } 
  else {
    return(format(kruskal,digits=3,scientific=T))
  }
}
#

####=========================================================================###
####=========================================================================###
#
####=========================================================================###
####=========================================================================###

# MISCELLANEOUS FUNCTIONS

####=========================================================================###
####=========================================================================###
#nabbed from http://kyrcha.info/2012/07/08/tutorials-fitting-a-sigmoid-function-in-r
sigmoid = function(params, x){
  1 / (1 + exp(-params[2] * (x - params[3])))
}
#
####=========================================================================###
####=========================================================================###

# HOMEPAGE PLOTS

####=========================================================================###
####=========================================================================###

#
####=========================================================================###
####=========================================================================###

# MAIN GENE PLOT

####=========================================================================###
####=========================================================================###
  # PLot for main gene
    render.comparison_plots = function(lit_review_pen, input_highlighted_variant, input_x_axis_selection,    input_y_axis_selection, input_colour_fill,  input_gene, input_freq,input_regression_line,input_acmg_selection,plot_height){
    # Filter dataset to just those where x-value,y-value and colour are not NA 
    # Need to find all variants that are non NA for colour, x-axis,y-axis and pathogenicity selection
    dataset_in_use <- lit_review_pen
    # Filter to gene of interest
    dataset_in_use <- dataset_in_use[dataset_in_use$gene == input_gene & !is.na(dataset_in_use$gene) & !is.na(dataset_in_use$HGVS),]
    # Filter to pathogenicity selection
    dataset_in_use <- dataset_in_use[dataset_in_use$acmg_literal %in% input_acmg_selection,]
    # Filter to X axis selection
    dataset_in_use <- dataset_in_use[!is.na(dataset_in_use[[comparison_plot_selection_options[[input_x_axis_selection]]]]),]
    # Filter to Y axis selection
    dataset_in_use <- dataset_in_use[!is.na(dataset_in_use[[comparison_plot_selection_options[[input_y_axis_selection]]]]),]
    # Filter to Colour Selection
    dataset_in_use <- dataset_in_use[!is.na(dataset_in_use[[comparison_plots_colour_columns[[input_colour_fill]]]]),]


 # Set x and y maxes
  axis_maxes=list(
    input_freq,
    100,
    100,
    100,
    100,
    input_freq,
    input_freq,
    input_freq,
    input_freq,
    input_freq,
    1,
    1,
    1,
    1,
    cadd_max,
    1
    )


  ###################################
  ###       SELECT X VALUES       ###
  ###################################

  use_x_values <- dataset_in_use[[comparison_plot_selection_options[[input_x_axis_selection]]]]
  use_x_max             <- axis_maxes[[input_x_axis_selection]]

  ###################################
  ###       SELECT Y VALUES       ###
  ###################################

  use_y_values          <- dataset_in_use[[comparison_plot_selection_options[[input_y_axis_selection]]]]
  use_y_max             <- axis_maxes[[input_y_axis_selection]]
  x_axis_label <- comparison_plot_axis_labels[[input_x_axis_selection]]
  y_axis_label <- comparison_plot_axis_labels[[input_y_axis_selection]]

  ###################################
  ###         SET SCALE BAR       ###
  ###################################

  # scale bar to appear in top right corner
  rect_r <- use_x_max
  rect_t<- 1.2*use_y_max
  rect_l <- 0.65*use_x_max
  rect_b<- 0.95*1.2*use_y_max

  ###################################
  ###         SET COLOUR          ###
  ###################################
  use_scale_bar_caption=list(c(0,0.5,1,"Penetrance"),c("Sporadic","50%","Familial",""),c(cadd_min,"",cadd_max,"CADD Phred Score"))
  use_colour_points <- dataset_in_use[[comparison_plots_colour_columns[[input_colour_fill]]]]
  highlighter_point_colour <- dataset_in_use[[comparison_plots_colour_columns[[input_colour_fill]]]][dataset_in_use$HGVS == input_highlighted_variant]
  #use_colour_points <- use_colour_point[[input_colour_fill]]
  scale_bar_caption <- use_scale_bar_caption[[input_colour_fill]]
  ###################################
  ###        SET CI POINTS        ###
  ###################################
  x_value_lower <-dataset_in_use[[comparison_plot_lower_bound_columns[[input_x_axis_selection]]]][dataset_in_use$HGVS == input_highlighted_variant]
  x_value_upper <-dataset_in_use[[comparison_plot_upper_bound_columns[[input_x_axis_selection]]]][dataset_in_use$HGVS == input_highlighted_variant]
  x_value_point <-dataset_in_use[[comparison_plot_selection_options[[input_x_axis_selection]]]][dataset_in_use$HGVS == input_highlighted_variant]  
  y_value_lower <-dataset_in_use[[comparison_plot_lower_bound_columns[[input_y_axis_selection]]]][dataset_in_use$HGVS == input_highlighted_variant]
  y_value_upper <-dataset_in_use[[comparison_plot_upper_bound_columns[[input_y_axis_selection]]]][dataset_in_use$HGVS == input_highlighted_variant]
  y_value_point <-dataset_in_use[[comparison_plot_selection_options[[input_y_axis_selection]]]][dataset_in_use$HGVS == input_highlighted_variant]  

  ###################################
  ###     CREATE PLOTS            ###
  ###################################

  plot_margins_variable
  gene_features_plot(
    use_x_values,
    use_y_values,
    use_x_max,
    # 1.2 on ymax adds space for header
    use_y_max*1.2,
    input_gene,
    use_colour_points,
    #use_x_literal_dataset,
    #use_x_data_type,
    #use_y_literal_dataset,
    #use_y_data_type,
    "lit_review_pen_ALS_penetrance",
    "lit_review_pen_ALS_penetrance.lower",
    "lit_review_pen_ALS_penetrance.upper",
    x_value_point, 
    y_value_point,
    #use_colour_point_list[[input_colour_fill]],
    highlighter_point_colour,
    x_value_lower,
    x_value_upper,
    y_value_lower,
    y_value_upper,
    input_highlighted_variant,
    rect_l,
    rect_b,
    rect_r,
    rect_t,
    input_regression_line,
    x_axis_label,
    y_axis_label,
    plot_height)
  scale_bar(rect_l,rect_b,rect_r,rect_t,1.2*use_y_max,100,scale_bar_caption,plot_height)
}
#
gene_features_plot = function (
  xvalues,
  yvalues,
  xmax,
  ymax,
  gene,
  colour,
  #x_literal_dataset,
  #x_data_type,
  #y_literal_dataset,
  #y_data_type,
  x_point_horizontal,
  x_point_vertical_low,
  x_point_vertical_high,
  x_point,
  y_point,
  point_colour,
  x_point_lower,
  x_point_upper,
  y_point_lower,
  y_point_upper,
  HGVS_CLICK,
  rect_l,
  rect_b,
  rect_r,
  rect_t,
  regression_line,
  x_axis_label,
  y_axis_label,
  plot_height){
  # create initial plot with data points 
  plot( 
    xvalues, 
    yvalues,
    cex=3*plot_height/700,
    pch=21,
    xaxt="n",yaxt="n",las=1,
    col="white",
    bg=colour,
    xlim=c(0,xmax),
    ylim=c(0,ymax),
    xlab=x_axis_label,
    ylab=y_axis_label,
    cex.lab=1.4*plot_height/700,
    cex.main=1.4*plot_height/700,
    col.main=nicegrey,
    col.lab=nicegrey,
    main=bquote(bold(italic(.(gene)) ~ " Properties"))
    )
  # plot point and confidence intervals on hover
  segments(
    x0=x_point_lower, 
    y0=y_point,
    x1=x_point_upper,
    y1=y_point,
    col="grey",
    lwd=5*plot_height/700
    )
  # 
  segments(
    x0=x_point,
    y0=y_point_lower,
    x1=x_point,
    y1=y_point_upper,
    col="grey",
    lwd=5*plot_height/700
    )
  # hover point
  points( 
    x_point,
    y_point,
    cex=3*plot_height/700,
    pch=21,
    col="black",
    bg=point_colour
    )
  # plot regression line 
  if(regression_line==TRUE){
    abline(lm(yvalues ~ xvalues),lwd=2*plot_height/700,col=nicegrey)
    # plot accompanying text in top left corner 
    text(0,ymax,"Regression:",col=nicegrey,pos=4,cex=1*plot_height/700)
    text(xmax/8,ymax,paste("R²: ",round(summary(lm(yvalues ~ xvalues))$r.squared,digits=4),sep=""),col=nicegrey,pos=4,cex=1*plot_height/700)
    text(xmax/8,0.95*ymax,paste("p: ",round(lmp(lm(yvalues ~ xvalues)),digits=4),sep=""),col=nicegrey,pos=4,cex=1*plot_height/700)
    text(xmax/8,0.9*ymax,paste("Slope: ",round(coef(lm(yvalues ~ xvalues))[2],digits=4),sep=""),col=nicegrey,pos=4,cex=1*plot_height/700)
  }
  # plot text for hover 
  if (length(HGVS_CLICK > 0)){
    text((rect_l+rect_r)/2,(rect_b-(0.15*ymax)),HGVS_CLICK,col=nicegrey,cex=1*plot_height/700)
  }
  # create axes 
  axis(side=1,lwd=1*plot_height/700,cex.axis=1*plot_height/700,col.axis=nicegrey,col=nicegrey)
  axis(side=2,lwd=1*plot_height/700,las=1,cex.axis=1*plot_height/700,col.axis=nicegrey,col=nicegrey,at = seq(0, ymax/1.2,by=(round((ymax/1.2)/4,3))))
  box(lwd=2*plot_height/700,col=nicegrey)
}
#
####=========================================================================###
####=========================================================================###

# HETEROGENEITY PLOT

####=========================================================================###
####=========================================================================###

render.heterogeneity_plots <- function(gene,variant,lit_review_pen,population_selection,phenotype_selection,family_selection){
  ###
  # Function to assess and plot geographic hetereogeneity for a single variant 
  ###
  # Filter population studies file to gene in question
  population_studies_gene       <- population_studies[population_studies$Gene==gene]
  if (nrow(population_studies_gene)>0){
  # Create fork of all population individuals
  population_people_dataset_use<-population_people_dataset
  # Account for phenotype and family history selection in population studies 
  population_studies_gene<-population_studies_gene[population_studies_gene$Phenotype==phenotype_selection,]
  population_people_dataset_use=population_people_dataset_use[!is.na(population_people_dataset_use$population_carriers_primary_phenotype) & population_people_dataset_use$population_carriers_primary_phenotype==phenotype_selection,]
  population_studies_gene<-population_studies_gene[population_studies_gene$History==family_selection,]
  population_people_dataset_use=population_people_dataset_use[!is.na(population_people_dataset_use$population_carriers_family_history) & population_people_dataset_use$population_carriers_family_history==family_selection,]
  # If population_selection is 1, this is a between continent comparison and needs slightly different treatment than within continent comparisons
  if (population_selection == 1){
    # Collapse population file based on continent 
    forest_file <- aggregate(population_studies_gene$Count, by=list(gene=population_studies_gene$Gene,Eth=population_studies_gene$Continent),FUN=sum)
    # Create blank column
    forest_file$Carrier_count <- 0
    # Sum the carriers from each continent 
    for (continent in unique(population_people_dataset_use$population_carriers_continent)){
      forest_file$Carrier_count[forest_file$Eth==continent] <- nrow(population_people_dataset_use[population_people_dataset_use$HGVS %in% variant & population_people_dataset_use$population_carriers_continent==continent & !is.na(population_people_dataset_use$population_carriers_continent),])
    }
    # Get the proportion of carriers for each continent 
    forest_file$proportion <- (forest_file$Carrier_count / forest_file$x)
    forest_plot<-forest_file 
    forest_plot<-forest_plot[with(forest_plot,rev(order(proportion,x))),] 
    plot_labels<-forest_plot$Eth
    for (number in 1:length(populations)){
      plot_labels<-gsub(populations[number],populations_literal[number],plot_labels)
    }
  # If a within continent study 
} 
else if (population_selection != 1) {
    # Collapse population file based based on country 
    forest_file <- aggregate(population_studies_gene$Count, by=list(Country=population_studies_gene$Country,gene=population_studies_gene$Gene,Eth=population_studies_gene$Continent),FUN=sum)
    # Create blank column
    forest_file$Carrier_count <- 0
    # Sum the carriers from each country 
    for (country in unique(population_people_dataset_use$population_carriers_nationality)){
      forest_file$Carrier_count[forest_file$Country==country] <- nrow(population_people_dataset_use[population_people_dataset_use$HGVS %in% variant & population_people_dataset_use$population_carriers_nationality==country & !is.na(population_people_dataset_use$population_carriers_nationality),])
    }
    # Get the proportion of carriers for each country 
    forest_file$proportion <- (forest_file$Carrier_count / forest_file$x)
    forest_plot<-forest_file 
    population<-populations[[population_selection-1]]
    forest_plot<-forest_file[forest_file$Eth==population,]
    forest_plot<-forest_plot[with(forest_plot,rev(order(proportion,x))),] 
    plot_labels<-forest_plot$Country
  }
  # Sort by highest proportion for better plotting 
  # Scenario 1 : >= 1 variant carrier across all continents or only one country / region studied
  if (nrow(forest_plot)>1 & sum(forest_plot$proportion)>0) {
    forest(metaprop(forest_plot$Carrier_count,forest_plot$x,studlab=plot_labels,control=list(maxiter=1000,stepadj=0.5)),print.tau2=FALSE,print.I2.ci=TRUE,print.pval.Q.LRT=TRUE,print.pval.Q=TRUE,LRT=TRUE,print.pval=FALSE,digits=5,xlim=c(0,round_any(max(forest_plot$proportion), 0.1, f = ceiling)))
    grid.text(paste(variant,"Heterogeneity Forest Plot",sep=" "), .5, 0.9, gp=gpar(cex=1))
  # If 0 variant carriers across regions - this cannot be meta analysed so is plotted as a table 
} 
else if ((nrow(forest_plot)>1 & sum(forest_plot$proportion)==0) | (nrow(forest_plot)==1)) {
  forest_plot<-forest_plot[rev(order(forest_plot$x)),] 
  if (population_selection != 1) {
    rownames<-forest_plot$Country
  } 
  else {
    rownames<-forest_plot$Eth
  }
  forest_plot<-subset(forest_plot,select=c("x","Carrier_count"))
  names(forest_plot)<-c("Combined Studies Size","Carrier Count")
  for (number in 1:length(populations)){
    rownames<-gsub(populations[number],populations_literal[number],rownames)
  }
  plot_margins_variable
  grid.table(forest_plot,rows=rownames,theme=ttheme_minimal())
    # If there are no country / continent studies for a gene selection then separate plot needed 
  } 
  else if (nrow(forest_plot)==0){
    plot_margins_variable
    plot(0.5,0.5,
      col="white",       
      pch=21,
      bg="white",      
      type="l",
      xaxt="n",
      yaxt="n",
      xlab="",
      ylab="",
      main="",
      cex.lab=1.4,
      cex.main=1.4,
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,100),
      ylim=c(0,100),bty='n')
    text(50,50,"No population studies recorded",col=nicegrey,adj=0)

  }
} 
else {
  plot_margins_variable
  plot(0.5,0.5,
    col="white",       
    pch=21,
    bg="white",      
    type="l",
    xaxt="n",
    yaxt="n",
    xlab="",
    ylab="",
    main="",
    cex.lab=1.4,
    cex.main=1.4,
    col.main=nicegrey,
    col.lab=nicegrey,
    xlim=c(0,100),
    ylim=c(0,100),bty='n')
  text(50,50,"No population studies recorded",col=nicegrey,adj=0)

}

}

#
####=========================================================================###
####=========================================================================###

# REGIONS PLOT FUNCTIONS

####=========================================================================###
####=========================================================================###

render.regions_analysis_plot=function(selection_region_browser,region_region_browser,regions_pathogenicity_selection,regions_phenotype_selection,regions_history_selection,fam_prop_slider,lit_review_pen,plot_height){
  # Get list of potential variants 
  regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
  ###########
  # If Selection is familial or sporadic
  ###########
  # If selection is familial or sporadic can calculate these directly 
  if (regions_history_selection=="Familial" | regions_history_selection=="Sporadic"){
    # Filter population individuals: phenotype of interest; history of interest ; variant matches 
    population_people_use    <- population_people_dataset[population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history==regions_history_selection & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
    # Filter population studies: phenotype of interest; history of interest 
    population_studies_use   <- population_studies[population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History==regions_history_selection & !is.na(population_studies$History),]
    #Filter to country or continent
    if (selection_region_browser=="Country"){  
      # People
      population_people_use  <- population_people_use[population_people_use$population_carriers_nationality==region_region_browser & !is.na(population_people_use$population_carriers_nationality),]
      # Studies 
      population_studies_use <- population_studies_use[population_studies_use$Country==region_region_browser & !is.na(population_studies_use$Country),]    
      } else if (selection_region_browser=="Continent"){
      # People
      population_people_use  <- population_people_use[population_people_use$population_carriers_continent==populations[match(region_region_browser,populations_literal)] & !is.na(population_people_use$population_carriers_continent),]
      # Studies 
      population_studies_use <- population_studies_use[population_studies_use$Continent==populations[match(region_region_browser,populations_literal)] & !is.na(population_studies_use$Continent),]     
    } 
    # Collapse Files If there is a file to collapse 
    if (nrow(population_studies_use)>0){
      population_studies_use_collapse <- aggregate(population_studies_use$Count,by=list(gene=population_studies_use$Gene),FUN=sum)
      # Create blank column
      population_studies_use_collapse$Carrier_count<-0
      # Sum the carriers from each country 
      for (Gene in population_studies_use_collapse$gene){
        population_studies_use_collapse$Carrier_count <- ifelse(population_studies_use_collapse$gene==Gene,nrow(population_people_use[population_people_use$gene==Gene,]),population_studies_use_collapse$Carrier_count)
      }
      # Get the proportion + CI of carriers for each country 
      population_studies_use_collapse$proportion <-binom.confint(x=population_studies_use_collapse$Carrier_count,n=population_studies_use_collapse$x,method='wilson')$mean
      population_studies_use_collapse$proportion.lower <-binom.confint(x=population_studies_use_collapse$Carrier_count,n=population_studies_use_collapse$x,method='wilson')$lower
      population_studies_use_collapse$proportion.upper <-binom.confint(x=population_studies_use_collapse$Carrier_count,n=population_studies_use_collapse$x,method='wilson')$upper
      # Filter file to pathogenic or LP genes 
      if (regions_pathogenicity_selection == as.numeric(1)) {
        population_studies_use_collapse=population_studies_use_collapse[population_studies_use_collapse$gene %in% path_genes,]
      } 
      else if ((regions_pathogenicity_selection == as.numeric(2))|(regions_pathogenicity_selection == as.numeric(3))) {
        population_studies_use_collapse=population_studies_use_collapse[population_studies_use_collapse$gene %in% p_or_lp_genes,]

      }
      # Reorder
      population_studies_use_collapse <- population_studies_use_collapse[with(population_studies_use_collapse,rev(order(proportion,x,gene))),]
      # To plot overall proportion of individuals explained need to calculate familial and sporadic separately and then merge the binomial proportions nad confidence intervals 
    }
     # Collapse population studies on gene and familial status 
 plot_margins_variable
 # Create the following plot if user chooses Sporadic or Familial or if they choose overall and there are both sporadic and familial population studies available.
 #if ((regions_history_selection=="Familial" | regions_history_selection=="Sporadic") | (exists('population_studies_use_fam')==T & exists('population_studies_use_spor')==T & nrow(population_studies_use_fam)>0 & nrow(population_studies_use_spor)>0)) {
  # sending up a plot space that goes to 100 on y and 140 on x. 140 on x allows me to put spaces between blocks while maintaining 1-100 scale 
  plot(50,50,
    col="black",       
    pch=21,
    bg="black",      
    type="l",
    xaxt="n",
    yaxt="n",
    xlab="",
    ylab="",
      #main="",
      cex.lab=1.4*plot_height/700,
      cex.main=1.4*plot_height/700,
      main=paste("Percentage of ",regions_history_selection," ",regions_phenotype_selection," cases explained in ",region_region_browser,sep=""),
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,50),
      ylim=c(0,180),bty='n')
    #text(50,50,"No Pedigrees Recorded",col=nicegrey,adj=0)
    gap_between <- 3
    start_pos <- 5
    end_pos   <- 10
    bottom<-180
    top<-180
    unexplained<-100
    #number_of_boxes<-nrow(population_studies_use_collapse[population_studies_use_collapse$proportion>0,])
    if (exists("population_studies_use_collapse")==TRUE){
      number_of_boxes<-nrow(population_studies_use_collapse)
      for (i in number_of_boxes:1){
        bottom <- top-100*population_studies_use_collapse[i,]$proportion
        gene_use<- population_studies_use_collapse[i,]$gene
        colour<-vangogh_palette[match(gene_use,p_or_lp_genes)]
        rect(start_pos,bottom,end_pos,top,col=colour,border="white")
        # plot text
        mid_point<-(bottom+top)/2
        confint<-paste(": ",round(100*population_studies_use_collapse[i,]$proportion,2),"% (",round(100*population_studies_use_collapse[i,]$proportion.lower,2),"-",round(100*population_studies_use_collapse[i,]$proportion.upper,2),")",sep="")
        text(end_pos+1,mid_point,bquote(italic(.(gene_use)) ~ .(confint)),col=nicegrey,adj=0,cex=0.65*plot_height/700)
        top<-bottom-gap_between
      }
      # plot unexplained proportion
      unexplained<-100*(1-sum(population_studies_use_collapse$proportion))
    }
    if (unexplained>0){
      bottom<-top-unexplained
      rect(start_pos,bottom,end_pos,top,col="grey70",border="white")
      mid_point<-(top+bottom)/2
      text(end_pos+1,mid_point,paste("No Variant",round(unexplained,2),"%"),col=nicegrey,adj=0,cex=0.65*plot_height/700)
    }

}
  
  else if (regions_history_selection=="Overall"){ 
    # Filter population individuals: phenotype of interest; history of interest ; variant matches 
    population_people_use_fam    <- population_people_dataset[population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Familial" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
    population_people_use_spor   <- population_people_dataset[population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Sporadic" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
    # Filter population studies: phenotype of interest; history of interest 
    population_studies_use_fam   <- population_studies[population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Familial" & !is.na(population_studies$History),]
    population_studies_use_spor   <- population_studies[population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Sporadic" & !is.na(population_studies$History),]
    #Filter to country or continent
    if (selection_region_browser=="Country"){  
      population_people_use_fam  <- population_people_use_fam[population_people_use_fam$population_carriers_nationality==region_region_browser & !is.na(population_people_use_fam$population_carriers_nationality),]
      population_people_use_spor  <- population_people_use_spor[population_people_use_spor$population_carriers_nationality==region_region_browser & !is.na(population_people_use_spor$population_carriers_nationality),] 
      population_studies_use_fam <- population_studies_use_fam[population_studies_use_fam$Country==region_region_browser & !is.na(population_studies_use_fam$Country),]    
      population_studies_use_spor <- population_studies_use_spor[population_studies_use_spor$Country==region_region_browser & !is.na(population_studies_use_spor$Country),]    
    } 
    else if (selection_region_browser=="Continent"){
      population_people_use_fam  <- population_people_use_fam[population_people_use_fam$population_carriers_continent==populations[match(region_region_browser,populations_literal)] & !is.na(population_people_use_fam$population_carriers_continent),]
      population_people_use_spor  <- population_people_use_spor[population_people_use_spor$population_carriers_continent==populations[match(region_region_browser,populations_literal)] & !is.na(population_people_use_spor$population_carriers_continent),]
      population_studies_use_fam <- population_studies_use_fam[population_studies_use_fam$Continent==populations[match(region_region_browser,populations_literal)] & !is.na(population_studies_use_fam$Continent),]           
      population_studies_use_spor <- population_studies_use_spor[population_studies_use_spor$Continent==populations[match(region_region_browser,populations_literal)] & !is.na(population_studies_use_spor$Continent),]     
    }
    if (exists('population_studies_use_fam')==T & exists('population_studies_use_spor')==T & nrow(population_studies_use_fam)>0 & nrow(population_studies_use_spor)>0){
    # Collapse Files
    if (nrow(population_studies_use_fam)>0){
      population_studies_use_collapse_fam <- aggregate(population_studies_use_fam$Count,by=list(gene=population_studies_use_fam$Gene),FUN=sum)
      colnames(population_studies_use_collapse_fam) <- c("gene","x_fam")
    # Create blank column
    population_studies_use_collapse_fam$Carrier_count_fam<-0 
    for (Gene in population_studies_use_collapse_fam$gene){
      population_studies_use_collapse_fam$Carrier_count_fam <- ifelse(population_studies_use_collapse_fam$gene==Gene,nrow(population_people_use_fam[population_people_use_fam$gene==Gene,]),population_studies_use_collapse_fam$Carrier_count_fam)
    }
    # Calculate binomial proportion, upper and lower bound for familial
    population_studies_use_collapse_fam$proportion_fam <-binom.confint(x=population_studies_use_collapse_fam$Carrier_count_fam,n=population_studies_use_collapse_fam$x_fam,method='wilson')$mean
    population_studies_use_collapse_fam$proportion.lower_fam <-binom.confint(x=population_studies_use_collapse_fam$Carrier_count_fam,n=population_studies_use_collapse_fam$x_fam,method='wilson')$lower
    population_studies_use_collapse_fam$proportion.upper_fam <-binom.confint(x=population_studies_use_collapse_fam$Carrier_count_fam,n=population_studies_use_collapse_fam$x_fam,method='wilson')$upper
  }
  if (nrow(population_studies_use_spor)>0){
    population_studies_use_collapse_spor <- aggregate(population_studies_use_spor$Count,by=list(gene=population_studies_use_spor$Gene),FUN=sum)
    colnames(population_studies_use_collapse_spor) <- c("gene","x_spor")
    # Create blank column
    population_studies_use_collapse_spor$Carrier_count_spor<-0
    for (Gene in population_studies_use_collapse_spor$gene){
      population_studies_use_collapse_spor$Carrier_count_spor <- ifelse(population_studies_use_collapse_spor$gene==Gene,nrow(population_people_use_spor[population_people_use_spor$gene==Gene,]),population_studies_use_collapse_spor$Carrier_count_spor)
    }
    # Calculate binomial proportion, upper and lower bound for sporadic
    population_studies_use_collapse_spor$proportion_spor <-binom.confint(x=population_studies_use_collapse_spor$Carrier_count_spor,n=population_studies_use_collapse_spor$x_spor,method='wilson')$mean
    population_studies_use_collapse_spor$proportion.lower_spor <-binom.confint(x=population_studies_use_collapse_spor$Carrier_count_spor,n=population_studies_use_collapse_spor$x_spor,method='wilson')$lower
    population_studies_use_collapse_spor$proportion.upper_spor <-binom.confint(x=population_studies_use_collapse_spor$Carrier_count_spor,n=population_studies_use_collapse_spor$x_spor,method='wilson')$upper
  }
    # Merge familial and sporadic 
    population_studies_use_collapse <- merge(population_studies_use_collapse_fam,population_studies_use_collapse_spor,all.x=T,all.y=T)
    # Need to make some corrections for NA values after merging: set NA to 0 for proportion explained and lower bound and 100 for upper bound (i.e. completely unknown)
    population_studies_use_collapse$proportion_fam<-ifelse(is.na(population_studies_use_collapse$proportion_fam)==TRUE,0,population_studies_use_collapse$proportion_fam)
    population_studies_use_collapse$proportion.lower_fam<-ifelse(is.na(population_studies_use_collapse$proportion.lower_fam)==TRUE,0,population_studies_use_collapse$proportion.lower_fam)
    population_studies_use_collapse$proportion.lower_fam<-ifelse(population_studies_use_collapse$proportion.lower_fam<0,0,population_studies_use_collapse$proportion.lower_fam)
    population_studies_use_collapse$proportion.upper_fam<-ifelse(is.na(population_studies_use_collapse$proportion.upper_fam)==TRUE,100,population_studies_use_collapse$proportion.upper_fam)    
    population_studies_use_collapse$proportion_spor<-ifelse(is.na(population_studies_use_collapse$x_spor)==TRUE,0,population_studies_use_collapse$proportion_spor)
    population_studies_use_collapse$proportion.lower_spor<-ifelse(is.na(population_studies_use_collapse$proportion.lower_spor)==TRUE,0,population_studies_use_collapse$proportion.lower_spor)
    population_studies_use_collapse$proportion.lower_spor<-ifelse(population_studies_use_collapse$proportion.lower_spor<0,0,population_studies_use_collapse$proportion.lower_spor)
    population_studies_use_collapse$proportion.upper_spor<-ifelse(is.na(population_studies_use_collapse$proportion.upper_spor)==TRUE,100,population_studies_use_collapse$proportion.upper_spor)    
    # Weight the binomial proportions by familial proportions slider 
    population_studies_use_collapse$proportion <-(fam_prop_slider*population_studies_use_collapse$proportion_fam) + ((1-fam_prop_slider)*population_studies_use_collapse$proportion_spor)
    # Weight the lower CIs 
    population_studies_use_collapse$proportion.lower <-(fam_prop_slider*population_studies_use_collapse$proportion.lower_fam) + ((1-fam_prop_slider)*population_studies_use_collapse$proportion.lower_spor)
    # Weight the upper CIs
    population_studies_use_collapse$proportion.upper <-(fam_prop_slider*population_studies_use_collapse$proportion.upper_fam) + ((1-fam_prop_slider)*population_studies_use_collapse$proportion.upper_spor)
    # Restrict to 0-1
    population_studies_use_collapse$proportion.upper[population_studies_use_collapse$proportion.upper>1 & !is.na(population_studies_use_collapse$proportion.upper)] <- 1
    population_studies_use_collapse$proportion.lower[population_studies_use_collapse$proportion.lower<0 & !is.na(population_studies_use_collapse$proportion.lower)] <- 0
    if (regions_pathogenicity_selection == as.numeric(1)) {
      population_studies_use_collapse=population_studies_use_collapse[population_studies_use_collapse$gene %in% path_genes,]
      } else if ((regions_pathogenicity_selection == as.numeric(2))|(regions_pathogenicity_selection == as.numeric(3))) {
        population_studies_use_collapse=population_studies_use_collapse[population_studies_use_collapse$gene %in% p_or_lp_genes,]

      }
      population_studies_use_collapse <- population_studies_use_collapse[with(population_studies_use_collapse,rev(order(proportion,gene))),]

    }
  #}
 # Collapse population studies on gene and familial status 
 plot_margins_variable
 # Create the following plot if user chooses Sporadic or Familial or if they choose overall and there are both sporadic and familial population studies available.
 if (exists('population_studies_use_fam')==T & exists('population_studies_use_spor')==T & nrow(population_studies_use_fam)>0 & nrow(population_studies_use_spor)>0) {
  # sending up a plot space that goes to 100 on y and 140 on x. 140 on x allows me to put spaces between blocks while maintaining 1-100 scale 
  plot(50,50,
    col="black",       
    pch=21,
    bg="black",      
    type="l",
    xaxt="n",
    yaxt="n",
    xlab="",
    ylab="",
      #main="",
      cex.lab=1.4*plot_height/700,
      cex.main=1.4*plot_height/700,
      main=paste("Percentage of ",regions_history_selection," ",regions_phenotype_selection," cases explained in ",region_region_browser,sep=""),
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,50),
      ylim=c(0,180),bty='n')
    #text(50,50,"No Pedigrees Recorded",col=nicegrey,adj=0)
    gap_between <- 3
    start_pos <- 5
    end_pos   <- 10
    bottom<-180
    top<-180
    unexplained<-100
    #number_of_boxes<-nrow(population_studies_use_collapse[population_studies_use_collapse$proportion>0,])
    if (exists("population_studies_use_collapse")==TRUE){
      number_of_boxes<-nrow(population_studies_use_collapse)
      for (i in number_of_boxes:1){
        bottom <- top-100*population_studies_use_collapse[i,]$proportion
        gene_use<- population_studies_use_collapse[i,]$gene
        colour<-vangogh_palette[match(gene_use,p_or_lp_genes)]
        rect(start_pos,bottom,end_pos,top,col=colour,border="white")
        # plot text
        mid_point<-(bottom+top)/2
        confint<-paste(": ",round(100*population_studies_use_collapse[i,]$proportion,2),"% (",round(100*population_studies_use_collapse[i,]$proportion.lower,2),"-",round(100*population_studies_use_collapse[i,]$proportion.upper,2),")",sep="")
        text(end_pos+1,mid_point,bquote(italic(.(gene_use)) ~ .(confint)),col=nicegrey,adj=0,cex=0.65*plot_height/700)
        top<-bottom-gap_between
      }
      # plot unexplained proportion
      unexplained<-100*(1-sum(population_studies_use_collapse$proportion))
    }
    if (unexplained>0){
      bottom<-top-unexplained
      rect(start_pos,bottom,end_pos,top,col="grey70",border="white")
      mid_point<-(top+bottom)/2
      text(end_pos+1,mid_point,paste("No Variant",round(unexplained,2),"%"),col=nicegrey,adj=0,cex=0.65*plot_height/700)
    }

}
else {
  blank_plot(100,100)
  text(0,80,"Familial Studies or Sporadic Studies Unavailable for the Given Region",col=nicegrey,adj=0,cex=1*plot_height/700)
}

  }
}


####=========================================================================###
####=========================================================================###

# PLOTTING FUNCTIONS

####=========================================================================###
####=========================================================================###

# function to get p value from regression line, nabbed from https://gist.github.com/stephenturner/722049
lmp = function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
#
scale_bar = function(rect_l,rect_b,rect_r,rect_t,ymax,n_cols,caption,plot_height){
  for(i in 1:n_cols){
    rect((rect_l+((rect_r-rect_l)*i/n_cols)),rect_b,rect_l+((rect_r-rect_l)*(1+i)/n_cols),rect_t,col=cols(n_cols)[i],border=cols(n_cols)[i])
  }
  text(c(rect_l,(rect_l+rect_r)/2,rect_r,(rect_l+rect_r)/2),c(rep((rect_b-(0.035*ymax)),3),(rect_b-(0.07*ymax))),caption,col=nicegrey,cex=1*plot_height/700)
}
#
# Functions for proportion explained 
blank_plot=function(xlim,ylim){
  # Function to create blank 
  plot(xlim/2,ylim/2,
    col="white",       
    pch=21,
    bg="white",      
    type="l",
    xaxt="n",
    yaxt="n",
    xlab="",
    ylab="",
    #main="",
    cex.lab=1.4,
    cex.main=1.4,
    main="",
    col.main=nicegrey,
    col.lab=nicegrey,
    xlim=c(0,xlim),
    ylim=c(0,ylim),
    bty='n')
}
#
####=========================================================================###
####=========================================================================###

# PEDIGREE PLOT

####=========================================================================###
####=========================================================================###
pedigree_call_variants=function(pedigree,input_variant){
  if (nchar(pedigree)==0){
    return(input_variant)
  }
  else {
    ped<-read.csv(paste("pedigrees/",pedigree,".ped",sep=""),header=T,sep='\t')
    if("V3" %in% colnames(ped)){
      return(paste(sort(unique(ped$V1)),sort(unique(ped$V2)),sort(unique(ped$V3)),sep=';'))
    }
    else if ("V2" %in% colnames(ped)){
      return(paste(sort(unique(ped$V1)),sort(unique(ped$V2)),sep=';'))
    }
    else {
      return(sort(unique(ped$V1)))
    }
  }
}


render.pedigree_plot = function(pedigree,plot_height){
  # Create blank plot if no pedigree
  if (nchar(pedigree)==0){
    plot_margins_variable
    # Create blank plot 
    blank_plot(100,100)
    plot(0.5,0.5,
      col="white",       
      pch=21,
      bg="white",      
      type="l",
      xaxt="n",
      yaxt="n",
      xlab="",
      ylab="",
      main="",
      cex.lab=1.4*plot_height/1200,
      cex.main=1.4*plot_height/1200,
      col.main=nicegrey,
      col.lab=nicegrey,
      xlim=c(0,100),
      ylim=c(0,100),bty='n')
    text(2,50,"No Pedigrees Recorded",col=nicegrey,adj=0,cex=plot_height/1200)
    #box(lwd=2,col=nicegrey)
  }
  # read in and plot pedigree
  else {
    ped<-read.csv(paste("/Users/markdoherty/Documents/GitHub/journALS/pedigrees/",pedigree,".ped",sep=""),header=T,sep='\t')
    ped_df<-data.frame(MND=ped$MND, FTD_Dementia=ped$FTD_Dementia, Other=ped$Other) 
    pedAll<- pedigree(id=ped$id,dadid=ped$dadid,momid=ped$momid,sex=ped$Sex,as.matrix(ped_df))
    legendPlot(pedAll,cex=min(1,2*plot_height/1200),symbolsize=min(1,2*plot_height/1200),col.label=NULL,affected.label=NULL,align=T,id=gsub('NA|NA;NA|NA;NA;NA','',paste(ped$AgeOrAgeOfOnset,ped$Genotype,sep="\n")))
    ###
    # If sex of parents is ambigous need to correct for this
    # Each ped file has an 'actual sex' and an 'assumed sex' column. If these differ for any user, then I need to take action
    ###
    # If there are any that do not match
    if (length(which((ped$Sex==ped$Actual_Sex)==FALSE)) >0){
      # Get plot output
      plot_output<-(plot(pedAll,cex=min(1,2*plot_height/1200),symbolsize=min(1,2*plot_height/1200),col.label=NULL,affected.label=NULL,align=T,id=gsub('NA|NA;NA|NA;NA;NA','',paste(ped$AgeOrAgeOfOnset,ped$Genotype,sep="\n"))))
      legendPlot(pedAll,cex=min(1,2*plot_height/1200),symbolsize=min(1,2*plot_height/1200),col.label=NULL,affected.label=NULL,align=T,id=gsub('NA|NA;NA|NA;NA;NA','',paste(ped$AgeOrAgeOfOnset,ped$Genotype,sep="\n")))
      # Get non matching indices
      non_match_indices<-which((ped$Sex==ped$Actual_Sex)==FALSE)
      for (indice in non_match_indices){
        # centre point 
        x_value<-plot_output$x[[indice]]
        # top point 
        y_value<-plot_output$y[[indice]]
        # plot a white rectangle
        #rect(x_value-(plot_output$boxw/2),y_value+plot_output$boxh,x_value+(plot_output$boxw/2),y_value,col="white",border="white",lwd=min(5,5*plot_height/1200))
                rect(x_value-(plot_output$boxw/2),y_value+plot_output$boxh,x_value+(plot_output$boxw/2),y_value,col="white",border="white",lwd=4)
        # draw diamond
        # draw four diagonals - definitely an easier way to do this 
        lines(x=c(x_value-(plot_output$boxw/2),x_value),c(y_value+(plot_output$boxh/2),y_value),lty=1,col="black",lwd=1)
        lines(x=c(x_value-(plot_output$boxw/2),x_value),c(y_value+(plot_output$boxh/2),y_value+(plot_output$boxh)),lty=1,col="black",lwd=1)
        lines(x=c(x_value,x_value+(plot_output$boxw/2)),c(y_value+(plot_output$boxh),y_value+(plot_output$boxh/2)),lty=1,col="black",lwd=1)
        lines(x=c(x_value,x_value+(plot_output$boxw/2)),c(y_value,y_value+(plot_output$boxh/2)),lty=1,col="black",lwd=1)
        # plot internal lines
        # find mid point height
        centre_x <- x_value
        centre_y <- y_value + plot_output$boxh/2 
        # line1 
        # slope
        m1<- tan(-13*pi/180)
        # y=m1*x-m1*centre_x+centre_y # equation for line 1 
        # line2
        m2<- plot_output$boxh/plot_output$boxw
        y_intersect=(m2*(-m1*centre_y/m2+m1*m2*centre_x/m2-m1*m2*plot_output$boxw/(2*m2)-m1*centre_x+centre_y))/(m2-m1)
        x_intersect=x=(y_intersect-centre_y+m2*centre_x-m2*plot_output$boxw/2)/m2
        x_intersect_right <- centre_x+plot_output$boxw/2-x_intersect+centre_x-plot_output$boxw/2
        # vertical line
        lines(x=c(x_value,x_value),c(y_value,y_value+(plot_output$boxh/2)),lty=1,col="black",lwd=1)
        # left line 
        lines(x=c(x_intersect,centre_x),c(y_intersect,centre_y),lty=1,col="black",lwd=1)
        # right line 
        lines(x=c(x_intersect_right,centre_x),c(y_intersect,centre_y),lty=1,col="black",lwd=1)
        # Fill polygon
        if (ped_df[indice,]$MND==1){
          polygon(x=c(centre_x-plot_output$boxw/2,x_value,centre_x,x_intersect),y=c(centre_y,y_value,centre_y,y_intersect),density=NA,col="black",border=NA,lwd=min(1,plot_height/1200))
        }
        if (ped_df[indice,]$FTD_Dementia==1){
          polygon(x=c(centre_x,x_value,centre_x+plot_output$boxw/2,x_intersect_right),y=c(centre_y,y_value,centre_y,y_intersect),density=NA,col="black",border=NA,lwd=min(1,plot_height/1200))
        }
        if (ped_df[indice,]$Other==1){
          polygon(x=c(x_intersect,centre_x,x_intersect_right,x_value),y=c(y_intersect,centre_y,y_intersect,y_value+plot_output$boxh),density=NA,col="black",border=NA,lwd=min(1,plot_height/1200))
        }
      }
    }

  }
}

pedigree_meioses = function(pedigree,input_variant){
   # Calculate segregation evidence
   # Three methods of calculating segregation
   # 1) Counting meioses:Jarvik GP, Browning BL. Consideration of Cosegregation in the Pathogenicity Classification of Genomic Variants. Am J Hum Genet. 2016;98(6):1077–81
   # 2) Full likelihood bayesian methods: [Thompson D, Easton DF, Goldgar DE. A full-likelihood method for the evaluation of causality of sequence variants from family data. Am J Hum Genet. 2003;73(3):652–5]
   # 3) Cosegregation likelihood ratios: [Mohammadi L, et al. A simple method for co-segregation analysis to evaluate the pathogenicity of unclassified variants; BRCA1 and BRCA2 as an example. BMC Cancer.]
   # FLB and CSLR are superior methods however they require 1) pentrance classes for genes 2) Ages for people in pedigrees. 
   # Both are not available for most of our data. 
   # Counting meioses performs similarly to FLB and CSLR for pathogenic variants but cannot identify benign but may be the only option for our data. [Rañola JMO, Liu Q, Rosenthal EA, Shirts BH. A comparison of cosegregation analysis methods for the clinical setting.]
   # Cutoff settings suggested by Jarvik et. al. 
   #                       Single Family     >1 Families
   # Strong evidence:      >1/16             >1/8 
   # Moderate evidence:    ≤1/16             ≤1/8 
   # Supporting evidence:  ≤1/8              ≤1/4 
   # Given the incomplete and age related penetrance of known ALS related variants, as per Jarvik et al a conservative approach is used which only takes into account affected individuals
   # Assume model is homozygous if all affected genotyped indivudals are homozygous, else assume heterozygous


   ped<-read.csv(paste("pedigrees/",pedigree,".ped",sep=""),header=T,sep='\t')
   # Create dummy columns for software to run, age and y.born don't affect result
   ped$Age<-0
   ped$y.born<-0
   ped$female<-ifelse(ped$Sex=="2",1,0)
   ped$status<-ifelse(ped$MND=="1" | ped$FTD_Dementia =="1", TRUE,FALSE)
   # make sure calculating right variant (pedigrees have up to three segregating variants)
   if (sort(unique(ped$V1))==input_variant){
     # check if homozygous - if row is a case and V1_A1 and V1_A2 are both 1
     if(all(ped$V1_A1[ped$V1_A1=="1" & (ped$MND=="1" | ped$FTD_Dementia=="1")]==ped$V1_A2[ped$V1_A1=="1" & (ped$MND=="1" | ped$FTD_Dementia=="1")])==TRUE){
       #for homozygous genotype can take A2
       inheritance_pattern<-"homozygous"
       ped$genotype<- ped$V1_A2
     } 
     else {
       inheritance_pattern<-"heterozygous"
       ped$genotype<- ped$V1_A1
     }
   } 
   if ("V2" %in% colnames(ped)){
     if (sort(unique(ped$V2))==input_variant){
       if(all(ped$V2_A1[ped$V2_A1=="1" & (ped$MND=="1" | ped$FTD_Dementia=="1")]==ped$V2_A2[ped$V2_A1=="1" & (ped$MND=="1" | ped$FTD_Dementia=="1")])==TRUE){
         inheritance_pattern<-"homozygous"
         ped$genotype<- ped$V2_A2
       } 
       else {
         inheritance_pattern<-"heterozygous"
         ped$genotype<- ped$V2_A1
       }
     }
   } 
   if ("V3" %in% colnames(ped)){
     if (sort(unique(ped$V3))==input_variant){
       if(all(ped$V3_A1[ped$V3_A1=="1" & (ped$MND=="1" | ped$FTD_Dementia=="1")]==ped$V3_A2[ped$V3_A1=="1" & (ped$MND=="1" | ped$FTD_Dementia=="1")])==TRUE){
         inheritance_pattern<-"homozygous"
         ped$genotype<- ped$V3_A2
       }
       else {
         inheritance_pattern<-"heterozygous"
         ped$genotype<- ped$V3_A1
       }
     }
   }
   if (is.null(CalculateLikelihoodRatio(ped=ped, affected.vector=ped$status)$separating.meioses)==TRUE){
     meioses<-"not available"
     denominator<-"not available"
     plural<-"not available"
   }
   else{
     meioses<-CalculateLikelihoodRatio(ped=ped, affected.vector=ped$status)$separating.meioses
     if (inheritance_pattern=="heterozygous"){
       denominator<-2^meioses
     }
     else if (inheritance_pattern=="homozygous"){
       denominator<-4^meioses
     }
     if (meioses==1){
       plural<-"meiosis"
     }
     else {
       plural<- "meioses"
     }
   }
   this_list=list(
     "inheritance_pattern"=inheritance_pattern,
     "meioses"=meioses,
     "denominator"=denominator,
     "plural"=plural
     )
   return(this_list)
 }
 pedigree_line_2_text=function(pedigree,variant_variant_browser){
  #if no pedigree for variant
  if (nchar(pedigree)==0){
    return(HTML(paste("<b>Note:</b> No pedigrees observed for ",variant_variant_browser,sep="")))
  }
  #if pedigrees for variant
  else {
    #if no segregation
    if(pedigree_meioses(pedigree,variant_variant_browser)$meioses=="not available"){
      return(HTML(paste("<b>Note:</b> The observed pattern of inheritance for pedigree ",pedigree," is not consistent with a single origin. Meioses cannot be calculated for this pedigree.",sep="")))
    } 
    else if ((paste(pedigree,variant_variant_browser,sep='_') %in% failing_ped_vars)){
      # if variant does not segregate
      return(HTML(paste("<b>Note:</b> Variant ", variant_variant_browser, " does not segregate with disease in pedigree ",pedigree, ". This provides supporting benign evidence.",sep="")))
    } 
    else {
    #if segregation
    return(HTML(paste("<b>Note:</b> Pedigree ", pedigree, " has ",pedigree_meioses(pedigree,variant_variant_browser)$meioses, " ",pedigree_meioses(pedigree,variant_variant_browser)$plural,". Under a ",pedigree_meioses(pedigree,variant_variant_browser)$inheritance_pattern," model this has a probability of 1/",pedigree_meioses(pedigree,variant_variant_browser)$denominator,sep="")))
  }

}
}

find_failing_ped_vars=function(){
  # list of pedigrees
  files <- list.files(path="pedigrees", pattern="*ped", full.names=FALSE, recursive=FALSE)
  # empty vector for nonsegregating ped/vars
  ped_var_nonseg<-c()
  for (entry in files){
    # read each pedigree
    ped<-read.csv(paste("pedigrees/",entry,sep=""),header=T,sep='\t')
    # filter pedigree to known cases 
    ped1<-ped[(ped$MND==1 | ped$FTD_Dementia) & (!is.na(ped$MND) & !is.na(ped$FTD_Dementia)),]
    # filter pedigree to known cases that are homozygous reference (note: can just filter on V1_A1 as no one would be entered as 0/1)
    ped2<-ped1[ped1$V1_A1==0,]
    # if any rows remain then...
    if (nrow(ped2) > 0){
      # variant
      var<-as.character(unique(sort(ped$V1)))
      ped_var<-gsub(".ped","",paste(entry,var,sep='_'))
      # add to list 
      ped_var_nonseg<-c(ped_var_nonseg,ped_var)
    }
    # if pedigree has a second variant (V2)
    if ("V2" %in% colnames(ped)){
      ped2<-ped1[ped1$V2_A1==0,]
      if (nrow(ped2) > 0){
        var<-as.character(unique(sort(ped$V2)))
        ped_var<-gsub(".ped","",paste(entry,var,sep='_'))
        ped_var_nonseg<-c(ped_var_nonseg,ped_var)
      }
    }
    # if pedigree has a third variant (V3)
    if ("V3" %in% colnames(ped)){
      ped2<-ped1[ped1$V3_A1==0,]
      if (nrow(ped2) > 0){
        var<-as.character(unique(sort(ped$V3)))
        ped_var<-gsub(".ped","",paste(entry,var,sep='_'))
        ped_var_nonseg<-c(ped_var_nonseg,ped_var)
      }
    }
  }
  return(ped_var_nonseg)
}

pedigree_line_3_text=function(pedigree,lit_review_pen,variant_variant_browser){
  if (nchar(pedigree)==0){
    return(HTML(paste("<b>Note:</b> No pedigrees are observed for ",variant_variant_browser,". This contributes no evidence of pathogenicity.",sep="")))
  }
  else {
    return(HTML(paste("<b>Note:</b> A segregation probability of 1/",pedigree_meioses_count(variant_variant_browser,lit_review_pen)$denominator," is observed across ",pedigree_meioses_count(variant_variant_browser,lit_review_pen)$pedigree_count," pedigrees for variant ",variant_variant_browser,". This contributes ",pedigree_meioses_count(variant_variant_browser,lit_review_pen)$pathogenic," evidence of pathogenicity.",sep="")))
  }
}

pedigree_meioses_count = function(input_variant,lit_review_pen){
  # see description on pedigree_meioses() function
  # i= denominator
  i=1
  # j= number of families 
  j=0
  # list all pedigrees assocaited with this variant
  ped_list<-unique(separate_rows(lit_review_pen[lit_review_pen$HGVS==input_variant & is.na(lit_review_pen$HGVS)==FALSE,],all_carriers_pedigree,sep="\\|",convert=T)$all_carriers_pedigree)[!is.na(unique(separate_rows(lit_review_pen[lit_review_pen$HGVS==input_variant & is.na(lit_review_pen$HGVS)==FALSE,],all_carriers_pedigree,sep="\\|",convert=T)$all_carriers_pedigree))]
  # Two clauses here. Clause 1 ensures it isn't a ped/variant combo that is non sergregating. Clause 2 ensures variant is present in ped. This would be relevant if a variant is in a proband but not described in the rest of the pedigree (only relevant to HLTF:c.1199C>T(p.[S400L]) variant in pedigree 32253937_1 )
  for (pedigree in ped_list){
    if ((!(paste(pedigree,input_variant,sep='_') %in% failing_ped_vars)) & (any(grep(input_variant,as.character(as.matrix(read.csv(paste("pedigrees/",pedigree,".ped",sep=""),header=T,sep='\t'))),fixed=T))==TRUE )) {
      denominator<-pedigree_meioses(pedigree,input_variant)$denominator
      if (length(denominator)>0&is.numeric(denominator)){
        i=i*denominator
        j<-j+1
      }
    } 
  }
  #single family
  if (j==1){
   if (i>=8 & i<16){
     pathogenic<-"supporting"
     acmg<-"supporting_pathogenic"
   }
   else if (i>=16 & i<32){
     pathogenic<-"moderate"
     acmg<-"moderate_pathogenic"
   }
   else if (i>32){
     pathogenic<-"strong"
     acmg<-"strong_pathogenic"
   }
   else {
     pathogenic<-"no"
     acmg<-NA
   }
 }
  #multiple families
  else if (j>1){
   if (i>=4 & i<8){
     pathogenic<-"supporting"
     acmg<-"supporting_pathogenic"
   }
   else if (i>=8 & i<16){
     pathogenic<-"moderate"
     acmg<-"moderate_pathogenic"
   }
   else if (i>=16){
     pathogenic<-"strong"
     acmg<-"strong_pathogenic"
   }
   else {
     pathogenic<-"no"
     acmg<-NA
   }
 }
 else {
   pathogenic<-"no"
   acmg<-NA
 }
 if(i>1){
   mei<-"meioses"
 }
 else {
   mei<-"meiosis"
 }
 this_list=list(
   "denominator"=i,
   "pedigree_count"=j,
   "plural"=mei,
   "pathogenic"=pathogenic,
   "acmg"=acmg
   )
 return(this_list)

}



# 
####=========================================================================###
####=========================================================================###

# TABLE PLOTS

####=========================================================================###
####=========================================================================###
table_add_top_and_bottom_lines = function(dataset){
  gtable_add_grob (dataset,
    grobs = grobTree(
      # top line 
      segmentsGrob(
        x0=unit(0,"npc"),
        y0=unit(1,"npc"),
        x1=unit(2,"npc"),
        y1=unit(1,"npc"),
        gp=gpar(lwd=1.0)),
      #bottom line 
      segmentsGrob(
        x0=unit(0,"npc"),
        y0=unit(0,"npc"),
        x1=unit(1,"npc"),
        y1=unit(0,"npc"),
        gp=gpar(lwd=1.0))),
    t=1,b=nrow(dataset),l=1,r=ncol(dataset))
}
#
table_add_lines_two_and_three = function(dataset){
  dataset <- gtable_add_grob (dataset,
    grobs = grobTree(
      segmentsGrob(
        # second line 
        x0=unit(0,"npc"),
        y0=unit(1,"npc"),
        x1=unit(2,"npc"),
        y1=unit(1,"npc"),
        gp=gpar(lwd=1.0)),
      #third line
      segmentsGrob(
        x0=unit(0,"npc"),
        y0=unit(0,"npc"),
        x1=unit(1,"npc"),
        y1=unit(0,"npc"),
        gp=gpar(lwd=1.0))),
    t=2,b=2,l=1,r=ncol(dataset))
}
#
table_plot = function(dataset,title_text,left_adjustment,HGVS_variant){
  table <- tableGrob(dataset,theme=tt,rows=NULL)
  title <- textGrob((paste(HGVS_variant,title_text,sep=" ")),gp=gpar(fontsize=7,fontface="italic"),hjust=left_adjustment)
  padding<-unit(0.5,"line")
  table <- gtable_add_rows(table,heights=grobHeight(title)+padding,pos=0)
  table <- gtable_add_grob(table, list(title),
   t=1, l=1, 
   r=ncol(table))
  table   <- table_add_top_and_bottom_lines(table) 
  table   <- table_add_lines_two_and_three(table) 
}
#
#nabbed from https://stackoverflow.com/questions/34530142/drop-down-checkbox-input-in-shiny
dropdownButton <- function(label = "", status = c("default", "primary", "success", "info", "warning", "danger"), ..., width = NULL) {

  status <- match.arg(status)
  # dropdown button content
  html_ul <- list(
    class = "dropdown-menu",
    style = if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), ";"),
    lapply(X = list(...), FUN = tags$li, style = "margin-left: 10px; margin-right: 10px;")
    )
  # dropdown button apparence
  html_button <- list(
    class = paste0("btn btn-", status," dropdown-toggle"),
    type = "button", 
    `data-toggle` = "dropdown"
    )
  html_button <- c(html_button, list(label))
  html_button <- c(html_button, list(tags$span(class = "caret")))
  # final result
  tags$div(
    class = "dropdown",
    do.call(tags$button, html_button),
    do.call(tags$ul, html_ul),
    tags$script(
      "$('.dropdown-menu').click(function(e) {
        e.stopPropagation();
        });")
    )
}
#
####=========================================================================###
####=========================================================================###

# UI FUNCTIONS

####=========================================================================###
####=========================================================================###

ui.regions_analysis_plot=function(){
  fluidRow(
    br(),
    br(),
    column(9,align="center",
      plotOutput("server.regions_analysis_plot",inline=T),
      ),
    column(3,
      selectInput("regions_pathogenicity_selection", "Filter Variants To:",
        c("Pathogenic Variants" = "1",
          "Pathogenic or Likely Pathogenic Variants"="2",
          "Reported Variants in Pathogenic or Likely Pathogenic Genes"="3")
        ),
      selectInput("regions_phenotype_selection", "Phenotype:",
        c("ALS",
          "FTD"
          )
        ),
      selectInput("regions_history_selection", "Family History:",
        c("Familial",
          "Sporadic",
          "Overall"
          )
        ),
      conditionalPanel(
        "input.regions_history_selection=='Overall'",
        sliderInput("fam_prop_slider", 
          h6("Estimated Proportion of Cases which are Familial:",style="font-weight: bold;"), 
          min = 0,
          max = 1, 
          value = 0.20)),
      downloadButton("download.regions_analysis_plot", label=h6("Save figure"),
        )
      )
    )
}

ui.heterogeneity_plots=function(){
  fluidRow(
    br(),
    br(),
    column(9,align="center",
      plotOutput("server.heterogeneity_plots")
      ),
    column(3,
      selectInput("heterogeneity_selection", "Population:",
        c("Continent Comparison" = "1",
          "Asia"="2",
          "Australasia"="3",
          "Europe"="4",
          "Middle East"="5",
          "North America"="6",
          "South America"="7",
          "Africa"="8")
        ),
      selectInput("heterogeneity_phenotype_selection", "Phenotype:",
        c("ALS"="ALS",
          "FTD"="FTD"),
        ),
      selectInput("heterogeneity_family_selection", "Family History:",
        c("Familial"="Familial",
          "Sporadic"="Sporadic")
        ),
      downloadButton("download.heterogeneity_plots", label=h6("Save figure")
        )
      )
    )
}

ui.age_plot=function(input_gene){
  # Create cumulative or KDE age plot 
  fluidRow(
    br(),
    br(),
    column(9,align="center",
      plotOutput("server.age_plot"),
      ),
    fluidRow(
      column(3,
        conditionalPanel(
          condition = "input.age_sex !='Comparison' &  input.age_fam != 'Comparison'",
          selectInput(
            "age_phenotype", 
          #choices=NULL,
          h6("Phenotype:",style="font-weight: bold;"),
          selected="All",
          c("All","ALS","FTD","Comparison")
          )
          ),
        conditionalPanel(
          condition = "input.age_sex =='Comparison' ||  input.age_fam =='Comparison'",
          selectInput(
            "age_phenotype2", 
            h6("Phenotype:",style="font-weight: bold;"),
            selected="All",
            c("All","ALS","FTD")
            )
          ),
        conditionalPanel(
         condition = "input.age_phenotype !='Comparison' &  input.age_fam != 'Comparison'",
         selectInput(
          "age_sex", 
          #choices=NULL,
          h6("Sex:",style="font-weight: bold;"),
          selected="All",
          c("All","Male","Female","Comparison")
          )
         ),
        conditionalPanel(
          condition = "input.age_phenotype =='Comparison' ||  input.age_fam == 'Comparison'",
          selectInput(
            "age_sex2", 
            h6("Sex:",style="font-weight: bold;"),
            selected="All",
            c("All","Male","Female")
            )
          ),
        conditionalPanel(
         condition = "input.age_phenotype !='Comparison' &  input.age_sex != 'Comparison'",
         selectInput(
          "age_fam", 
        #  choices=NULL,
        h6("Family History:",style="font-weight: bold;"),
        selected="All",
        c("All","Familial","Sporadic","Comparison")
        )
         ),
        conditionalPanel(
          condition = "input.age_phenotype =='Comparison' ||  input.age_sex == 'Comparison'",
          selectInput(
            "age_fam2", 
            h6("Family History:",style="font-weight: bold;"),
            selected="All",
            c("All","Familial","Sporadic")
            )
          ),
        selectInput(
          "age_plot_style",
          h6("Plot Type:",style="font-weight: bold;"),
          selected="Cumulative Distribution",
          c("Cumulative Distribution","Kernal Density Estimate")
          )
        )
      )
    )

}

ui.gene_plot=function(input_gene){
  # Plot the exons and introns of a gene and relevant variants
  fluidRow(
    br(),
    br(),
    column(9,align="center",
      plotOutput("server.gene_plot"),
      ),
    column(3,
     prettyCheckboxGroup(fill=TRUE,shape="square",status="default",icon = icon("check"),
       inputId = "gene_coordinates_plot_dataset_filter", 
       label = "Filter Display To Variants Present In:", 
       choiceNames = c("Literature","ALSdb","ALSVS_FAMILIAL","ProjectMinE"),
       choiceValues = c("Literature","ALSdb","ALSVS_FAMILIAL","ProjectMinE"),
       selected = NULL
       ),
      prettyCheckbox(fill=TRUE,shape="square",status="default",icon = icon("check"),
        inputId = "plot_introns", 
        label = "Display Introns"
        ),
      selectizeInput("transcript_selection", "Select Transcript:",
        choices=NULL,
        options=list(
          maxItems=1,
          placeholder='Select Transcript',
          create=TRUE)
        ),
      downloadButton("downloadcoordPlot", label=h6("Save figure")),
      )
    )
}
render.gene_plot=function(lit_review_pen,input_gene,input_gene_coordinates_plot_dataset_filter,input_plot_introns,input_transcript_selection,plot_height){
  # This code works but could be tidied up a lot 
  # Filter Dataset to Gene For Dot plot 
  dataset2<-lit_review_pen[lit_review_pen$gene==input_gene & !is.na(lit_review_pen$gene),]
  # Filter to variants in the literature for dot plot 
  dataset<-dataset2[dataset2$all_carriers_count>0 & !is.na(dataset2$all_carriers_count),]
  # Get exon coordinates 
  this_gene_coordinates <- gene_coordinates[gene_coordinates$name==input_transcript_selection]
  # 5' coordinate
  x_min<-this_gene_coordinates$txStart
  # 3' coordinate
  x_max<-this_gene_coordinates$txEnd
  # Find exon coordinates
  left<-as.numeric(unlist(strsplit(as.character(this_gene_coordinates$exonStarts), split=",")))
  right<-as.numeric(unlist(strsplit(as.character(this_gene_coordinates$exonEnds), split=",")))
  # Find coding exon coordinates
  coding_left<-left
  coding_right<-right
  remove_exons<-c()
  for (number in 1:length(coding_left)){
    # Find completely non coding exons
    if ((coding_left[number] < this_gene_coordinates$cdsStart & coding_right[number] < this_gene_coordinates$cdsStart) | (coding_left[number] > this_gene_coordinates$cdsEnd & coding_right[number]> this_gene_coordinates$cdsEnd)) {
      remove_exons<-append(remove_exons,number)
    }
    # Correct first exon
    if (coding_left[number]<this_gene_coordinates$cdsStart & coding_right[number]>this_gene_coordinates$cdsStart & coding_right[number]<this_gene_coordinates$cdsEnd){
      coding_left[number]<-this_gene_coordinates$cdsStart
    }
    # Correct last exon
    if (coding_left[number]>this_gene_coordinates$cdsStart & coding_right[number]>this_gene_coordinates$cdsEnd){
      coding_right[number] <- this_gene_coordinates$cdsEnd
    }
    # If only one exon 
    if (coding_left[number]<this_gene_coordinates$cdsStart & coding_right[number]>this_gene_coordinates$cdsEnd){
      coding_right[number] <- this_gene_coordinates$cdsEnd
      coding_left[number]<-this_gene_coordinates$cdsStart
    }
  }
  # Remove non coding exons 
  if (length(remove_exons)>0){
    coding_left<-coding_left[-remove_exons]
    coding_right<-coding_right[-remove_exons]
  }
  # Filter dataset to relevant cohorts 
  if ("Literature" %in% input_gene_coordinates_plot_dataset_filter){
    dataset2<-dataset2[dataset2$all_carriers_count>0 & !is.na(dataset2$all_carriers_count),]
  }
  if ("ALSdb" %in% input_gene_coordinates_plot_dataset_filter){
    dataset2<-dataset2[dataset2$ALSdb_AF>0 & !is.na(dataset2$ALSdb_AF),]
  }   
  if ("ALSVS_FAMILIAL" %in% input_gene_coordinates_plot_dataset_filter){
    dataset2<-dataset2[dataset2$ALSVS_FALS_AF>0 & !is.na(dataset2$ALSVS_FALS_AF),]
  }
  if ("ProjectMinE" %in% input_gene_coordinates_plot_dataset_filter){
    dataset2<-dataset2[dataset2$ProjectMinE_AF_cases>0 & !is.na(dataset2$ProjectMinE_AF_cases),]
  }
  # Find pathogenic/benign variants for bar plot - this is variants in the literature 
  p_variants <- as.numeric(dataset$pos[dataset$acmg_literal=="Pathogenic"])
  lp_variants <- as.numeric(dataset$pos[dataset$acmg_literal=="Likely Pathogenic"])
  vus_variants <-as.numeric(dataset$pos[dataset$acmg_literal=="VUS"])
  lb_variants <- as.numeric(dataset$pos[dataset$acmg_literal=="Likely Benign"])
  b_variants <- as.numeric(dataset$pos[dataset$acmg_literal=="Benign"])
  # Find pathogenic/benign variants for dot plot - this can be filtered by dataset
  p2_variants <- as.numeric(dataset2$pos[dataset2$acmg_literal=="Pathogenic"])
  lp2_variants <- as.numeric(dataset2$pos[dataset2$acmg_literal=="Likely Pathogenic"])
  vus2_variants <-as.numeric(dataset2$pos[dataset2$acmg_literal=="VUS"])
  lb2_variants <- as.numeric(dataset2$pos[dataset2$acmg_literal=="Likely Benign"])
  b2_variants <- as.numeric(dataset2$pos[dataset2$acmg_literal=="Benign"])
  # Adjust all coordinates if not plotting introns
  if (input_plot_introns==FALSE & as.numeric(this_gene_coordinates$exonCount) > 1){
    for (number in 1:(as.numeric(this_gene_coordinates$exonCount)-1)){
      # Only continue if intron is > 40 bp (leaving a 20bp border)
      if (left[number+1]-right[number] > 40){
        # Find intron start, end and length
        intron.start<-right[number]+20
        intron.end <-left[number+1]-20
        intron.length<-intron.end-intron.start
        # Remove variants falling in those introns from barplot and dot plot data 
        # Pathogenic 
        if (length(p_variants)>0){
          p_variants<-p_variants[p_variants<intron.start |  p_variants>intron.end]
          if (length(p_variants)>0){
            for (number in 1:length(p_variants)){
              if (p_variants[number]>intron.start){
                p_variants[number]=p_variants[number]-intron.length
              }
            }
          }
        }
        if (length(p2_variants)>0){
          p2_variants<-p2_variants[p2_variants<intron.start |  p2_variants>intron.end]
          if (length(p2_variants)>0){
            for (number in 1:length(p2_variants)){
              if (p2_variants[number]>intron.start){
                p2_variants[number]=p2_variants[number]-intron.length
              }
            }
          }
        }
        # Likely Pathogenic 
        if (length(lp_variants)>0){
          lp_variants<-lp_variants[lp_variants<intron.start |  lp_variants>intron.end]
          if (length(lp_variants)>0){
            for (number in 1:length(lp_variants)){
              if (lp_variants[number]>intron.start){
                lp_variants[number]=lp_variants[number]-intron.length
              }
            }
          }
        }
        if (length(lp2_variants)>0){
          lp2_variants<-lp2_variants[lp2_variants<intron.start |  lp2_variants>intron.end]
          if (length(lp2_variants)>0){
            for (number in 1:length(lp2_variants)){
              if (lp2_variants[number]>intron.start){
                lp2_variants[number]=lp2_variants[number]-intron.length
              }
            }
          }
        }
        # vus
        if (length(vus_variants)>0){
          vus_variants<-vus_variants[vus_variants<intron.start |  vus_variants>intron.end]
          if (length(vus_variants)>0){
            for (number in 1:length(vus_variants)){
              if (vus_variants[number]>intron.start){
                vus_variants[number]=vus_variants[number]-intron.length
              }
            }
          }
        }
        if (length(vus2_variants)>0){
          vus2_variants<-vus2_variants[vus2_variants<intron.start |  vus2_variants>intron.end]
          if (length(vus2_variants)>0){
            for (number in 1:length(vus2_variants)){
              if (vus2_variants[number]>intron.start){
                vus2_variants[number]=vus2_variants[number]-intron.length
              }
            }
          }
        }
        # Likely Benign 
        if (length(lb_variants)>0){
          lb_variants<-lb_variants[lb_variants<intron.start |  lb_variants>intron.end]
          if (length(lb_variants)>0){
            for (number in 1:length(lb_variants)){
              if (lb_variants[number]>intron.start){
                lb_variants[number]=lb_variants[number]-intron.length
              }
            }
          }
        }
        if (length(lb2_variants)>0){
          lb2_variants<-lb2_variants[lb2_variants<intron.start |  lb2_variants>intron.end]
          if (length(lb2_variants)>0){
            for (number in 1:length(lb2_variants)){
              if (lb2_variants[number]>intron.start){
                lb2_variants[number]=lb2_variants[number]-intron.length
              }
            }
          }
        }
        # Benign 
        if (length(b_variants)>0){
          b_variants<-b_variants[b_variants<intron.start |  b_variants>intron.end]
          if (length(b_variants)>0){
            for (number in 1:length(b_variants)){
              if (b_variants[number]>intron.start){
                b_variants[number]=b_variants[number]-intron.length
              }
            }
          }
        }
        if (length(b2_variants)>0){
          b2_variants<-b2_variants[b2_variants<intron.start |  b2_variants>intron.end]
          if (length(b2_variants)>0){
            for (number in 1:length(b2_variants)){
              if (b2_variants[number]>intron.start){
                b2_variants[number]=b2_variants[number]-intron.length
              }
            }
          }
        }
        # Readjust x_max
        x_max <- x_max-intron.length
        # Readjust left coordinates
        for (number in 1:length(left)){
          if (left[number]>intron.start){
            left[number]=left[number]-intron.length
          }
        }
        for (number in 1:length(coding_left)){
          if (coding_left[number]>intron.start){
            coding_left[number]=coding_left[number]-intron.length
          }
        }
        # Readjust right coordinates
        for (number in 1:length(right)){
          if (right[number]>intron.start){
            right[number]=right[number]-intron.length
          }
        }  
        for (number in 1:length(coding_right)){
          if (coding_right[number]>intron.start){
            coding_right[number]=coding_right[number]-intron.length
          }
        } 
        # Readjust variant positions for barplot 
        dataset$pos<-ifelse(!is.na(dataset$pos) & dataset$pos>intron.start,dataset$pos -intron.length,dataset$pos)
      }
    }
  }
  # Create blank plots with arbitrary y from 0 to 165 
  plot_margins_variable
  plot(1,1,     
    lwd=3*plot_height/700,
    col="white",
    type="l",
    xaxt="n",
    yaxt="n",
    las=1,
    xlab="",
    ylab="",
    cex.lab=1.4,
    cex.main=1.4,
    col.main=nicegrey,
    col.lab="white",
    frame=FALSE,
    xlim=c(x_min,x_max),
    ylim=c(0,165),xpd=NA)
  # Plot central line
  lines(x=c(x_min,x_max),y=c(50,50),col="grey",lty=1,lwd=2*plot_height/700,xpd=NA)
  # Plot All Exons
  rect(xleft=left,ybottom=rep(48,length(left)),xright=right,rep(52,length(left)),border=NA,col="grey",xpd=NA)
  # Plot Coding Exons
  rect(xleft=coding_left,ybottom=rep(45,length(coding_left)),xright=coding_right,rep(55,length(coding_left)),border=NA,col="grey",xpd=NA)
  # max number of variants in an exon
  max_per_exon<-0
  for (number in 1:length(coding_left)){
    #count<-length
    count<-nrow(dataset[!is.na(dataset$pos) & (as.numeric(dataset$pos)>as.numeric(coding_left[number]) & as.numeric(dataset$pos<coding_right[number])),])
    if (count>max_per_exon){
      max_per_exon<-count
    }
  }
  colours_for_bars<-c(nicered,niceorange,lightgreen,nicelightblue,reallydarkblue)
  # plot histogram
  for (number in 1:this_gene_coordinates$exonCount){
    p_count <-length(p_variants[p_variants>as.numeric(coding_left[number]) & p_variants<as.numeric(coding_right[number])])
    lp_count <-length(lp_variants[lp_variants>as.numeric(coding_left[number]) & lp_variants<as.numeric(coding_right[number])])
    vus_count <-length(vus_variants[vus_variants>as.numeric(coding_left[number]) & vus_variants<as.numeric(coding_right[number])])
    lb_count <-length(lb_variants[lb_variants>as.numeric(coding_left[number]) & lb_variants<as.numeric(coding_right[number])])
    b_count <-length(b_variants[b_variants>as.numeric(coding_left[number]) & b_variants<as.numeric(coding_right[number])])
    left_pos<-coding_left[number]
    right_pos<- coding_right[number]
    bottom <- 60
    top <-60
    entry<-0
    for (call in list(p_count,lp_count,vus_count,lb_count,b_count)){
      entry<-entry+1
      height<-50*(call/max_per_exon)
      top<-top+height
      rect(xleft=left_pos,ybottom=bottom,xright=right_pos,ytop=top,border=NA,col=colours_for_bars[[entry]],xpd=NA)
      bottom<-bottom+height
    }
  }
  axis(side=2,lwd=1*plot_height/700,cex.axis=1*plot_height/700,col.axis=nicegrey,col=nicegrey,at = c(60,110),labels=c(0,max_per_exon),xpd=NA)
  title( ylab = "           Variant Count / Exon (in Literature)",cex.lab=1.4*plot_height/700,col.lab=nicegrey,xpd=NA)
  height<-40
  number<-0
  # Dot plot
  for (call in list(p2_variants,lp2_variants,vus2_variants,lb2_variants,b2_variants)){
    height<-height-5
    number=number+1
    points( 
      call,
      rep(height,length(call)),
      cex=1.2*plot_height/700,
      pch=21,
      col="white",
      bg=colours_for_bars[[number]],xpd=NA
      )
  }
  # Plot Legend 
  locations<-(x_max-0.3*(x_max-x_min))
  hop<-(0.02*(x_max-x_min))

  points(
    rep(locations,5),
    c(165,158,151,144,137),
    cex=1.2*plot_height/700,
    pch=21,
    col="white",
    #bg="blue")
    bg=colours_for_bars
    ,xpd=NA)
  height=165
  for (entry in c("Pathogenic","Likely Pathogenic","Variant of Uncertain Significance","Likely Benign","Benign")){
    text(locations+hop,height,entry,col=nicegrey,adj=0,cex=1*plot_height/700,xpd=NA)
    height=height-7
  }
  # Plot 5' & 3'
  if (gene_coordinates$strand[gene_coordinates$name==input_transcript_selection]=="+"){
    text(x=x_min-(x_max-x_min)/80,y=50,"5'",col=nicegrey,cex=1*plot_height/700,xpd=NA)
    text(x=x_max+(x_max-x_min)/80,y=50,"3'",col=nicegrey,cex=1*plot_height/700,xpd=NA)
  }
  if (gene_coordinates$strand[gene_coordinates$name==input_transcript_selection]=="-"){
    text(x=x_min-(x_max-x_min)/80,y=50,"3'",col=nicegrey,cex=1*plot_height/700,xpd=NA)
    text(x=x_max+(x_max-x_min)/80,y=50,"5'",col=nicegrey,cex=1*plot_height/700,xpd=NA)
  }
  
}

ui.comparison_plots=function(){
  fluidRow(
    br(),
    br(),
    column(9,align="center",
                plotOutput("server.comparison_plots",
      ),
                br(),
                br(),                br(),
                br(),                br(),
                br(),
                br(),
                br(),
                h6(mainPanel(htmlOutput("comparison_plot_caveat1"),width=9,align="left")),
                h6(mainPanel(htmlOutput("comparison_plot_caveat2"),width=9,align="left"))
                ),
    column(3,
      selectInput("gnomAD_selection", h6("X Axis:",style="font-weight: bold;"),
        c("Overall Allele Frequency (Literature)" = "1",
          "Familial ALS Case Frequency (%) (Literature)" = "2",
          "Sporadic ALS Case Frequency (%) (Literature)" = "3",
          "Familial FTD Case Frequency (%) (Literature)" = "4",
          "Sporadic FTD Case Frequency (%) (Literature)" = "5",
          "gnomAD Allele Frequency" = "6",
          "Project MinE Case Allele Frequency" = "7",
          "Project MinE Control Allele Frequency" = "8",
          "ALSdb Allele Frequency" = "9",
          "ALSVS Famlial Allele Frequency" = "10",
          "ALS Lifetime Penetrance (Population)" = "11",
          "FTD Lifetime Penetrance (Population)" = "12",
          "ALS & FTD Combined Lifetime Penetrance (Population)" = "13",
          "ALS Lifetime Penetrance (Familial)" = "14",
          "Adjusted CADD Score" = "15",
          "Proportion of Cases with Positive Family History"="16")
        ),
      selectInput("dataset_selection", h6("Y Axis:",style="font-weight: bold;"),
        c("Overall Allele Frequency (Literature)" = "1",
          "Familial ALS Case Frequency (%) (Literature)" = "2",
          "Sporadic ALS Case Frequency (%) (Literature)" = "3",
          "Familial FTD Case Frequency (%) (Literature)" = "4",
          "Sporadic FTD Case Frequency (%) (Literature)" = "5",
          "gnomAD Allele Frequency" = "6",
          "Project MinE Case Allele Frequency" = "7",
          "Project MinE Control Allele Frequency" = "8",
          "ALSdb Allele Frequency" = "9",
          "ALSVS Famlial Allele Frequency" = "10",
          "ALS Lifetime Penetrance (Population)" = "11",
          "FTD Lifetime Penetrance (Population)" = "12",
          "ALS & FTD Combined Lifetime Penetrance (Population)" = "13",
          "ALS Lifetime Penetrance (Familial)" = "14",
          "Adjusted CADD Score" = "15",
          "Proportion of Cases with Positive Family History"="16")
        ),

      selectInput("colour_fill", h6("Colour by:",style="font-weight: bold;"),
        c("Penetrance" = "1",
          "Familial Proportion"="2",
          "Adjusted CADD Score"="3")
        ),
      selectizeInput("variant_highlighter", h6("Highlight Variant:",style="font-weight: bold;"),
        choices=NULL,
        options=list(          
          maxItems=1,
          placeholder='Select Variant in Plot',
          create=TRUE)

        ),
      br(),
      prettyCheckboxGroup(fill=TRUE,shape="square",status="default",icon = icon("check"),
        inputId = "gene_plot_pathogenicity_filter", 
        label = "Variants to Include:", 
        choiceNames = list("Pathogenic","Likely Pathogenic","Variant of Uncertain Significance","Likely Benign","Benign"),
        choiceValues = list("Pathogenic","Likely Pathogenic","VUS","Likely Benign","Benign"),
        selected = c("Pathogenic","Likely Pathogenic","VUS")
        ),
      prettyCheckbox(fill=TRUE,shape="square",status="default",icon = icon("check"),
        inputId = "regression_line", 
        label = "Display Regression Line"
        ),
      conditionalPanel(
        #condition = "input.gnomAD_selection =='1'|| input.gnomAD_selection =='3' || input.dataset_selection =='1'",
        condition = "input.gnomAD_selection =='2' || input.dataset_selection =='1'",
        sliderInput("freq", 
          h6("Allele Frequency Range:",style="font-weight: bold;"), 
          min = 0,
          max = 1, 
          value = 0.01)
        ),
      downloadButton("download.comparison_plots", label=h6("Save figure")),
      )

)
}

gene_table_input=function(){
  fluidRow(
    column(3,
      dropdownButton(
        label = "Select Columns To Display", 
        status = "default", 
        width = 300,
        checkboxGroupInput(
          inputId = "visible_columns", 
          label = "Choose", 
          choices = columns_of_interest,
          selected = preselected_columns_gene
          )
        )
      ),
    column(3,
      dropdownButton(
        label = "Select Variant Type", 
        status = "default", 
        width = 300,
        checkboxGroupInput(
          inputId = "variant_type", 
          label = "Choose", 
          choices = impacts,
          selected = impacts
          )
        )
      ),
    column(3,
      dropdownButton(
        label = "Must Be Present In:", status = "default", width = 300,
        checkboxGroupInput(
          inputId = "dataset_type", 
          label = "Choose", 
          choices = c("Literature","ALSdb","ALSVS_FAMILIAL","ProjectMinE"),
          selected = "Literature"
          )
        )
      )
    )
}
#

region_region_choices=function(selection_region_browser){
  if (selection_region_browser=="Country"){
    return(dropdown_countries)
  } 
  else if (selection_region_browser=="Continent"){
    return(populations_literal)
  }
  else if (selection_region_browser=="Global"){
    return("Global")
  }
}
###XXXXXX
variant_highlighter_choices=function(lit_review_pen,input_x_axis_selection,input_gene,input_y_axis_selection,input_colour_fill,input_acmg_selection){
  # Need to find all variants that are non NA for colour, x-axis,y-axis and pathogenicity selection
  dataset_in_use <- lit_review_pen
  # Filter to gene of interest
  dataset_in_use <- dataset_in_use[dataset_in_use$gene == input_gene & !is.na(dataset_in_use$gene),]
  # Filter to pathogenicity selection
  dataset_in_use <- dataset_in_use[dataset_in_use$acmg_literal %in% input_acmg_selection,]
  # Filter to X axis selection
  dataset_in_use <- dataset_in_use[!is.na(dataset_in_use[[comparison_plot_selection_options[[input_x_axis_selection]]]]),]
  # Filter to Y axis selection
  dataset_in_use <- dataset_in_use[!is.na(dataset_in_use[[comparison_plot_selection_options[[input_y_axis_selection]]]]),]
  # Filter to Colour Selection
  dataset_in_use <- dataset_in_use[!is.na(dataset_in_use[[colour_selection_options[[input_colour_fill]]]]),]

  return(dataset_in_use$HGVS)
}
select_region_browser_input=function(){
  fluidRow(
    column(3,
      selectInput("selection_region_browser","Country or Continent:",
        c("Country","Continent","Global"))
      ),
    column(3,
      selectizeInput("region_region_browser", "Select Region of Interest:",
        choices=NULL,
        options=list(
          maxItems=1,
          placeholder='Select Region',
          create=TRUE)
        )
      )
    )
}
#
select_variant_browser_input=function(){
  fluidRow(
    column(3,
      selectizeInput("gene_variant_browser", "Select Gene of Interest:",
        choices=NULL,
        options=list(maxItems=1,placeholder='Select Gene',create=TRUE)
        )
      ),
    column(3,
      selectizeInput("variant_variant_browser", "Select Variant of Interest:",
        choices=NULL,
        options=list(
          maxItems=1,
          placeholder='Select Variant',
          create=TRUE)
        )
      ),
    column(3,
      br(),
      dropdownButton(
        label = "Filter Search To Variants Present In:", status = "default", width = 300,
        checkboxGroupInput(
          inputId = "filter_variant_browser", 
          label = "Choose", 
          choices = c("Literature","ALSdb","ALSVS_FAMILIAL","ProjectMinE"),
          selected = "Literature"
          )
        )
      ),
    column(3,
      br(),
      dropdownButton(
        label = "Variants to Include:", status = "default", width = 300,
        checkboxGroupInput(
          inputId = "pathogenicity_filter_variant_browser", 
          label = "Choose", 
          choices = c("Pathogenic","Likely Pathogenic","VUS","Likely Benign","Benign"),
          selected = c("Pathogenic","Likely Pathogenic","VUS")
          )
        )
      )

    )

}
#
variant_browser_general_information_input=function(){
  fluidRow(
    column(3,
      dropdownButton(
        label = "Select Columns To Display", status = "default", width = 300,
        checkboxGroupInput(inputId = "visible_columns_single_variant", label = "Choose", 
          choices = columns_of_interest,
          selected = preselected_columns
          )
        ),
      style="display:-webkit-inline-box;")
    )
}
#
variant_pedigree_tab_input=function(){
  fluidRow(
    column(3,
      selectizeInput("pedigree_selection", "Select Pedigree:",
        choices=NULL,
        options=list(
          maxItems=1,
          placeholder='Select Pedigree',
          create=TRUE)
        )
      )
    )
}
#
variant_browser_geography_input=function(){
  fluidRow(
    column(3,
      selectInput("geography_selection", "Phenotype:",
        c("ALS","FTD","ALS_FTD")
        )
      )
    )
}


variant_browser_phenotype_information_input=function(){
  fluidRow(
    column(3,
      dropdownButton(
        label = "Select Columns To Display", status = "default", width = 300,
        checkboxGroupInput(inputId = "visible_columns_single_variant_phenotype", label = "Choose", 
          choices = pheno_frame_colnames,
          selected = pheno_frame_colnames
          )
        ),
      style="display:-webkit-inline-box;")
    )
}


# Summary Tab 
summary_gene_tab_input=function(){
  fluidRow(
    column(3,
      selectizeInput("summary_gene_selection", "Select Gene:",
        choices=NULL,
        options=list(
          maxItems=1,
          placeholder='Select Gene',
          create=TRUE)
        )
      )
    )
}

#
####=========================================================================###
####=========================================================================###

# CREATE SOME USEFUL VARIABLES AND VECTORS

####=========================================================================###
####=========================================================================###
# Variable for plotting 
plot_margins_variable<<-"par(mar=c(6,6,2,1),las=0,mgp = c(4, 1, 0))"
###
# Some nice colours 
nicegrey        <<- "grey30"
niceblack       <<- "#3b444b"
nicedarkblue    <<- "#146D97"
nicered         <<- "#C02906"
niceorange      <<- "#FF7400"
nicelightblue   <<- "#1E96B9"
reallydarkblue  <<- "#01184E"
reallylightblue <<- "#CFF6F5"
darkgreen       <<- "#CFAA28"
lightgreen      <<- "#EDDC58"
# Colour fades
darkgreen_rgb_fade    <<- rgb(0.8117647,0.6666667,0.1568627,alpha=0.2)
nicegrey_rgb_fade     <<- rgb(0.3019608,0.3019608,0.3019608,alpha=0.2)
nicered_rgb_fade      <<- rgb(0.7529412, 0.1607843, 0.02352941,alpha=0.2)
nicedarkblue_rgb_fade <<- rgb(0.07843137,0.427451,0.5921569,alpha=0.2)
# Extra colours for use in regions analysis plot 
pink<<- "#f4ddf1"
green<<-"#4897a5"
lime<<-"#c1ef37"
medium_green<<-"#5bcfa5"
brown<<-"#7f4827"
#bright_yellow<<-"#f4ff00"
bright_yellow<<-"#f5db37"
lavendar<<-"#8da3d3"
algae<<-"#526708"
yellow1<<-"#cad874"
sunset<<-"#c24802"
purple<<-"#523386"
peach<<-"#dca995"
light_purple<<-"#624f94"
pink_purple<<-"#9783dc"
hull<<-"#CAA575"
peach2<<-"#d45150"
# This pallete is used for regions analysis plot 
#vangogh_palette<<-c(nicered,hull,peach,sunset,niceorange,lightgreen,darkgreen,reallylightblue,nicelightblue,nicedarkblue,reallydarkblue,pink,green,lime,medium_green,brown,bright_yellow,lavendar,algae,yellow1,purple,light_purple,pink_purple)
#vangogh_palette<<-c(nicered,hull,peach,sunset,niceorange,lightgreen,darkgreen,reallylightblue,nicelightblue,nicedarkblue,reallydarkblue,pink,green,lime,brown,bright_yellow,lavendar,algae,yellow1,purple,light_purple,pink_purple)
vangogh_palette<<-c(nicelightblue,nicedarkblue,reallydarkblue,pink,green,lime,brown,bright_yellow,lavendar,algae,yellow1,purple,light_purple,pink_purple,nicered,hull,peach,sunset,niceorange,lightgreen,darkgreen,reallylightblue)

# Baseline Risk
# baseline_risk <- 2.5e-3
# baseline_risk_familial <- 2.5e-4 #(1/10? need to think about this)
cols<<-colorRampPalette(c(nicedarkblue,"grey70",nicered))
# # Penetrance Columns
#pen_lower_cols  <-    c("lit_review_pen_colour_lower","ALSdb_pen_colour_lower",   "ALSVS_FALS_pen_colour_lower",  "ProjectMinE_pen_colour_lower")
#pen_point_cols  <-    c("lit_review_pen_colour",      "ALSdb_pen_colour",         "ALSVS_FALS_pen_colour",        "ProjectMinE_pen_colour",   "ProjectMinE_pen_colour") #check this extra PM?
#pen_upper_cols  <-    c("lit_review_pen_colour_upper","ALSdb_pen_colour_upper",   "ALSVS_FALS_pen_colour_upper",  "ProjectMinE_pen_colour_upper")
#pen_points      <-    c("lit_review_pen_ALS_penetrance",       "ALSdb_pen_point",          "ALSVS_FALS_pen_point",         "ProjectMinE_pen_point")
#pen_lowers      <-    c("lit_review_pen_ALS_penetrance.lower",           "ALSdb_lower",              "ALSVS_FALS_lower",             "ProjectMinE_lower")
#pen_uppers      <-    c("lit_review_pen_ALS_penetrance.upper",           "ALSdb_upper",              "ALSVS_FALS_upper",             "ProjectMinE_upper")
# Comparison Plot labels and columns
#literal_datasets<-    c("Literature",                 "ALSdb",                    "ALSVS_Familial",               "Project MinE")
#AC_columns      <-    c("population_carriers_count",  "ALSdb_Allele_Count",       "ALSVS_FALS_AC",                "ProjectMinE_AC_cases")
#AF_columns      <-    c("lit_review_AF",              "ALSVS_FALS_AF","ALSdb_AF", "ProjectMinE_AF_cases",         "ProjectMinE_AF_controls")
#af_faf_ylabel   <-    c("Literature Allele Frequency","ALSdb Allele Frequency" ,  "ALSVS FALS Allele Frequency",  "ProjectMinE Case Allele Frequency","ProjectMinE Control Allele Frequency") #check this extra PM?
#####
# Universal Variables for comparison plot 
#####
# Axis Labels 
x_data_types    <-    c("Allele Count",             "Allele Frequency",             "Penetrance",         "CADD Phred Score")
y_data_types=c("Allele Frequency","Allele Count","Penetrance","Familial Proportion")
# Dataset Labels 
x_literal_datasets=list(
    #"gnomAD",
    "gnomAD",
    "gnomAD",
    "Literature",
    ""
    )
comparison_plots_colour_columns=list(
  "lit_review_pen_colour",
  "familial_proportion_colour",
  "cadd_colour"
  )


comparison_plot_lower_bound_columns<-list(
  "lit_review_AF",
  "FALS_Case_Frequency",
  "SALS_Case_Frequency",
  "FFTD_Case_Frequency",
  "SFTD_Case_Frequency",
  "gnomAD_non_neuro_use_AF",
  "ProjectMinE_AF_cases",
  "ProjectMinE_AF_controls",
  "ALSdb_AF",
  "ALSVS_FALS_AF",
  "lit_review_pen_ALS_penetrance.lower",
  "lit_review_pen_FTD_penetrance.lower",
  "lit_review_pen_ALSFTD_penetrance.lower",
  "penetrance_famhist.lower",
  "cadd_scaled",
  "familial_proportion_lower"
  )

comparison_plot_upper_bound_columns<-list(
  "lit_review_AF",
  "FALS_Case_Frequency",
  "SALS_Case_Frequency",
  "FFTD_Case_Frequency",
  "SFTD_Case_Frequency",
  "gnomAD_non_neuro_use_AF",
  "ProjectMinE_AF_cases",
  "ProjectMinE_AF_controls",
  "ALSdb_AF",
  "ALSVS_FALS_AF",
  "lit_review_pen_ALS_penetrance.upper",
  "lit_review_pen_FTD_penetrance.upper",
  "lit_review_pen_ALSFTD_penetrance.upper",
  "penetrance_famhist.upper",
  "cadd_scaled",
  "familial_proportion_upper"
  )

comparison_plot_axis_labels<-list(
  "Overall Allele Frequency (Literature)",
  "Familial ALS Case Frequency (%) (Literature)",
  "Sporadic ALS Case Frequency (%) (Literature)",
  "Familial FTD Case Frequency (%) (Literature)",
  "Sporadic FTD Case Frequency (%) (Literature)",
  "gnomAD Allele Frequency",
  "Project MinE Case Allele Frequency",
  "Project MinE Control Allele Frequency",
  "ALSdb Allele Frequency",
  "ALSVS Famlial Allele Frequency",
  "ALS Lifetime Penetrance (Population)",
  "FTD Lifetime Penetrance (Population)",
  "ALS & FTD Combined Lifetime Penetrance (Population)",
  "ALS Lifetime Penetrance (Familial)",
  "Adjusted CADD Score",
  "Proportion of Cases with Positive Family History")


comparison_plot_selection_options<-list(
  "lit_review_AF",
  "FALS_Case_Frequency",
  "SALS_Case_Frequency",
  "FFTD_Case_Frequency",
  "SFTD_Case_Frequency",
  "gnomAD_non_neuro_use_AF",
  "ProjectMinE_AF_cases",
  "ProjectMinE_AF_controls",
  "ALSdb_AF",
  "ALSVS_FALS_AF",
  "lit_review_pen_ALS_penetrance",
  "lit_review_pen_FTD_penetrance",
  "lit_review_pen_ALSFTD_penetrance",
  "penetrance_famhist_point",
  "cadd_scaled",
  "familial_proportion"
  )

# Lists for HGVS dropdown functions
# gnomAD_selection_options<-list(
#   "gnomAD_non_neuro_use_AC",
#   "gnomAD_non_neuro_use_AF",
#   "lit_review_pen_ALS_penetrance",
#   "cadd_scaled"
#   )
# dataset_selection_options<-list(
#   "lit_review_AF",
#   "population_carriers_count",
#   "lit_review_pen_ALS_penetrance",
#   "familial_proportion")
colour_selection_options<-list(
  "lit_review_pen_colour",
  "familial_proportion_colour",
  "cadd_colour"
  )
###
# Universal Variables for heterogeneity
###
# Populations for heterogeneity plot 
# Names of continents in the files
populations     <-    c("ASIA","AUS","EUR","ME","NAM","SAM","AFR")
# Literal names to use for plotting 
populations_literal <- c("Asia","Australasia","Europe","Middle East","North America","South America","Africa")

pheno_frame_colnames<<-c("HGVS","identifier","Nationality","Ethnicity","Ethnicity: Explanation","Primary Phenotype","Detailed Phenotype","Family History","Sex","Age of Onset","Disease Duration","Onset Site","Zygosity","Cognitive Impairment (Y/N)","de novo","de novo: confirmed parentage","Concurrent Variants","PMID","Comments","Pedigree ID","Continent")
preselected_columns <- c(
  "HGVS",

  #"HGVS2",
  "identifier",
  "impact",
  "gene",
  #"gene2",
  "ProjectMinE_AF_cases",
  "gnomAD_non_neuro_use_AF",
  #"gnomAD_controls_use_AF",
  "ALSVS_FALS_AF",
  "ALSdb_AF",
  "lit_review_AF",
  "acmg_literal")

preselected_columns_gene <- c(
  "HGVS",
  #"HGVS2",
  "identifier",
  "impact",
  "gene",
  #"gene2",
  "ProjectMinE_AF_cases",
  "gnomAD_non_neuro_use_AF",
  #"gnomAD_controls_use_AF",
  "ALSVS_FALS_AF",
  "ALSdb_AF",
  "lit_review_AF",
  "acmg_literal",
  "all_carriers_pmid",
  "all_carriers_aoo",
  "all_carriers_primary_phenotype",
  "all_carriers_nationality")
  #"cadd_scaled",
  #"all_carriers_count",
  #"lit_review_pen_point",
  #"ALSdb_pen_point",
  #"ALSVS_FALS_pen_point",
  #"ProjectMinE_pen_point")
# Scale bar labels
#use_title_point=list("Point Estimate","Point Estimate","CADD Phred Score")
#use_title_upper=list("Upper Confidence Bound","Upper Confidence Bound","")
#use_title_lower=list("Lower Confidence Bound","Lower Confidence Bound","")



columns_of_interest = c(
  "HGVS",
  #"HGVS2",
  "identifier",
  "impact",
  "familial_count",
  "sporadic_count",
  "ProjectMinE_AC_cases",
  "ProjectMinE_Alleles_cases",
  "ProjectMinE_AF_cases",
  "gnomAD_non_neuro_use_AC",
  "gnomAD_non_neuro_use_AN",
  "gnomAD_non_neuro_use_AF",
  #"gnomAD_controls_use_AC",
  #"gnomAD_controls_use_AN",
  #"gnomAD_controls_use_AF",
  "ALSVS_FALS_AC",
  "ALSVS_FALS_AN",
  "ALSVS_FALS_AF",
  "ALSdb_Allele_Count","ALSdb_AN","ALSdb_AF",
  "population_carriers_count",
  "lit_review_AN",
  "lit_review_AF",
  "gene",
  #"gene2",
  "cadd_scaled",
  "all_carriers_count",
  "all_carriers_pmid",
  "all_carriers_aoo",
  "all_carriers_disease_duration",
  "all_carriers_primary_phenotype",
  "all_carriers_detailed_phenotype",
  "all_carriers_denovo",
  "all_carriers_denovo_confirmed",
  "all_carriers_nationality",
  "all_carriers_ethnicity",
  "all_carriers_ethnicity_explanation",
  "all_carriers_sex",
  "all_carriers_onset_site",
  "all_carriers_onset_site_explanation",
  "all_carriers_family_history",
  "all_carriers_cognitive_impairment",
  "all_carriers_zygosity",
  "all_carriers_concurrent_variants",
  "all_carriers_pedigree",
  "all_carriers_comment",
  "population_carriers_nationality",
  "population_carriers_continent",
  "population_carriers_family_history",
  "population_carriers_zygosity",
  "ProjectMinE_AF_all",
  "ProjectMinE_AF_controls",
  "ProjectMinE_Alleles_all",
  "ProjectMinE_Alleles_controls",
  "ProjectMinE_allHOM_A1_PARTICIPANTS",
  "ProjectMinE_allHET_PARTICIPANTS",
  "ProjectMinE_allHOM_A2_PARTICIPANTS",
  "ProjectMinE_allMISSING_PARTICIPANTS",
  "ProjectMinE_allTOTAL_PARTICIPANTS",
  "ProjectMinE_casesHOM_A1_PARTICIPANTS",
  "ProjectMinE_casesHET_PARTICIPANTS",
  "ProjectMinE_casesHOM_A2_PARTICIPANTS",
  "ProjectMinE_casesMISSING_PARTICIPANTS",
  "ProjectMinE_casesTOTAL_PARTICIPANTS",
  "ProjectMinE_controlsHOM_A1_PARTICIPANTS",
  "ProjectMinE_controlsHET_PARTICIPANTS",
  "ProjectMinE_controlsHOM_A2_PARTICIPANTS",
  "ProjectMinE_controlsMISSING_PARTICIPANTS",
  "ProjectMinE_controlsTOTAL_PARTICIPANTS",
  "ALSVS_FALS_het",
  "ALSVS_FALS_homAlt",
  "ALSVS_FALS_homRef",
  "ALSdb_Sample_Count",
  "ALSdb_Hom_Count",
  "ALSdb_Heteroz_Count",
  "ALSdb_Ref_Count",
  "gnomAD_genomes_controls_AC",
  "gnomAD_genomes_controls_AN",
  "gnomAD_exomes_controls_AC",
  "gnomAD_exomes_controls_AN",
  "gnomAD_exomes_controls_use_AC",
  "gnomAD_genomes_controls_use_AC",
  "gnomAD_exomes_controls_use_AN",
  "gnomAD_genomes_controls_use_AN",
  "population_carriers_HOM_cases",
  "population_carriers_HET_cases",
  "all_carriers_HOM_cases",
  "all_carriers_HET_cases",
  "ProjectMinE_AC_controls",
  "lit_review_pen_ALS_penetrance",
  "lit_review_pen_ALS_penetrance.lower",
  "lit_review_pen_ALS_penetrance.upper",
  "ALSdb_pen_point",
  "ALSdb_lower",
  "ALSdb_upper",
  "ALSVS_FALS_pen_point",
  "ALSVS_FALS_lower",
  "ALSVS_FALS_upper",
  "ProjectMinE_pen_point",
  "ProjectMinE_lower",
  "ProjectMinE_upper",
  "familial_proportion",
  "familial_proportion_lower",
  "familial_proportion_upper",
  "all_carriers.HOM.cases",
  "all_carriers.HET.cases",
  "acmg_literal")


impacts<<-c("3_prime_UTR_variant","5_prime_UTR_variant","bidirectional_gene_fusion","conservative_inframe_deletion","conservative_inframe_insertion","disruptive_inframe_deletion","disruptive_inframe_insertion","downstream_gene_variant","exon_loss_variant","frameshift_variant","gene_fusion","initiator_codon_variant","intron_variant","missense_variant","non_coding_transcript_exon_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","start_lost","stop_gained","stop_lost","stop_retained_variant","structural_interaction_variant","synonymous_variant","upstream_gene_variant","Other")
