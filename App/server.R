# run as runApp("App/")
setwd("/journALS_Databrowser")


################################################################################
################################################################################

# Clear Memory

################################################################################
################################################################################

rm(list = ls())
source("App/functions.R")
################################################################################
################################################################################

#Load Packages

################################################################################
################################################################################

# load packages
load_packages()
# read in gnomAD rare variant counts per gene 
gnomAD_rare_variant_count <<- fread("Datafiles/gnomAD_ped_gene_rare_variant_count.gz",header=TRUE,sep='\t',na.strings=c("","NA")) 
# population file loads details of population studies
population_studies      <<- fread("Datafiles/rev14_population_studies_long_format.tsv",header=TRUE,sep='\t',na.strings=c("","NA"))
# read in coordinates of genes and exons
gene_coordinates      <<- fread("Datafiles/rev14_gene_coordinates.tsv",header=TRUE,sep='\t',na.strings=c("","NA"))
lit_review_pen<- fread("Datafiles/journALS_post_analysis.tsv.gz",header=TRUE,sep='\t',na.strings=c("","NA"))


# create dataset
all_ages_df <<- data.frame(fread("Datafiles/rev14_journALS_all_ages_df.tsv",header=TRUE,sep='\t',na.strings=c("","NA")))
failing_ped_vars <<- find_failing_ped_vars()

age_p_cutoff <<-  data.frame(fread("Datafiles/rev14_p_cutoff_df.tsv",header=TRUE,sep='\t',na.strings=c("","NA")))
age_p_cutoff <<- age_p_cutoff$p_cutoff

all_ages_df <<- data.frame(fread("Datafiles/rev14_journALS_all_ages_df.tsv",header=TRUE,sep='\t',na.strings=c("","NA")))
population_people_dataset <<- population_people_dataset_function(lit_review_pen)

dropdown_countries <<- all_countries(lit_review_pen)

# Dataframe of phenotypes of individuals found in the literature 
pheno_frame <<- create_dataframe_of_individuals(lit_review_pen)

# These are needed for regions analysis plot
# Find list of pathogenic variants
pathogenic_variants<<-lit_review_pen$HGVS[lit_review_pen$acmg_literal=="Pathogenic" & !is.na(lit_review_pen$acmg_literal)]
#
# Find list of pathogenic or likely variants
p_or_lp_variants<<-lit_review_pen$HGVS[(lit_review_pen$acmg_literal=="Pathogenic" | lit_review_pen$acmg_literal=="Likely Pathogenic") & !is.na(lit_review_pen$acmg_literal)]
# Reported variants in P ir LP genes 
path_genes<<- unique(lit_review_pen$gene[(lit_review_pen$acmg_literal=="Pathogenic") & !is.na(lit_review_pen$acmg_literal)])
p_or_lp_genes<<-unique(lit_review_pen$gene[(lit_review_pen$acmg_literal=="Pathogenic" | lit_review_pen$acmg_literal=="Likely Pathogenic") & !is.na(lit_review_pen$acmg_literal)])
variants_in_p_or_lp_genes<<-lit_review_pen$HGVS[lit_review_pen$population_carriers_count>0 & !is.na(lit_review_pen$population_carriers_count) & lit_review_pen$gene %in% p_or_lp_genes & !is.na(lit_review_pen$gene)]
list_of_pathogenic_variants<<-list(pathogenic_variants,p_or_lp_variants,variants_in_p_or_lp_genes)
# Maximum CADD values (necessary for comparison plot)
# max values for comparison plot 
cadd_max <<- max(as.numeric(as.character(lit_review_pen$cadd_scaled)),na.rm=TRUE)
cadd_min <<- min(as.numeric(as.character(lit_review_pen$cadd_scaled)),na.rm=TRUE)

####=========================================================================###
####=========================================================================###
server <- function(session,input, output) {


	####=========================================================================###
	####=========================================================================###

	# ACMG PATHOGENICITY PLOT 

	####=========================================================================###
	####=========================================================================###
	
	screen_width <- reactive({ifelse(!is.null(input$innerWidth),input$innerWidth,0)})

	output$server.acmg_pathogenicity_plots <- renderPlot({render.acmg_pathogenicity_plots(selected_variant$selected_variant,lit_review_pen,min(600,screen_width()))}, height=reactive(min(600,screen_width())),width=reactive(min(600,screen_width())))

	####=========================================================================###
	####=========================================================================###

	# HETEROGENEITY PLOT

	####=========================================================================###
	####=========================================================================###
	output$server.heterogeneity_plots <- renderPlot({render.heterogeneity_plots(input$gene_variant_browser,selected_variant$selected_variant,lit_review_pen,as.numeric(input$heterogeneity_selection),input$heterogeneity_phenotype_selection,input$heterogeneity_family_selection)}, width=800,height=400)
	####=========================================================================###
	####=========================================================================###

	# Regions Analysis Plot

	####=========================================================================###
	####=========================================================================###

	output$server.regions_analysis_plot <- renderPlot({render.regions_analysis_plot(input$selection_region_browser,input$region_region_browser,input$regions_pathogenicity_selection,input$regions_phenotype_selection,input$regions_history_selection,input$fam_prop_slider,lit_review_pen,min(700,screen_width()))}, width=reactive(min(700,screen_width())),height=reactive(min(700,screen_width())*9/7))

	####=========================================================================###
	####=========================================================================###

	# GENE COMPARISON PLOTS

	####=========================================================================###
	####=========================================================================###
	
	output$server.comparison_plots <- renderPlot({render.comparison_plots(lit_review_pen,input$variant_highlighter,as.numeric(input$gnomAD_selection),as.numeric(input$dataset_selection),as.numeric(input$colour_fill),input$gene,input$freq,as.numeric(input$regression_line),input$gene_plot_pathogenicity_filter,min(700,screen_width()))}, width=reactive(min(700,screen_width())),height=reactive(min(700,screen_width())*5/7))
	
	# Gene Layout Plot 

	output$server.gene_plot <- renderPlot({render.gene_plot(lit_review_pen,input$gene,input$gene_coordinates_plot_dataset_filter,input$plot_introns,input$transcript_selection,min(700,screen_width()))}, width=reactive(min(700,screen_width())),height=reactive(min(700,screen_width())*5/7))

	####=========================================================================###
	####=========================================================================###

	# AGE PLOT

	output$server.age_plot<-renderPlot({render.age_plot(input$age_phenotype,input$age_phenotype2,input$age_sex,input$age_sex2,input$age_fam,input$age_fam2,selected_variant$selected_variant,input$gene_variant_browser,input$age_plot_style,min(600,screen_width()))},height=reactive(min(600,screen_width())*2/3),width=reactive(min(600,screen_width())))
	
	####=========================================================================###
	####=========================================================================###

	# PEDIGREE PLOTS

	####=========================================================================###
	####=========================================================================###
  
	output$server.pedigree_plot <- renderPlot({render.pedigree_plot(input$pedigree_selection,min(1200,screen_width())*0.95)}, height=reactive((min(1200,screen_width())*6/12)*0.95),width=reactive(min(1200,screen_width())*0.95))

	####=========================================================================###
	####=========================================================================###

	# Summary Plots

	####=========================================================================###
	####=========================================================================###
  
  	output$server.summary.proportion_explained.plot <- renderPlot({render.summary.proportion_explained.plot(lit_review_pen,min(1200,screen_width()))}, width=reactive(min(1200,screen_width())),height=reactive(min(1200,screen_width())*25/120))	#output$server.summary.proportion_explained.plot <- renderPlot({render.summary.proportion_explained.plot(lit_review_pen)}, width=1200,height=250)


	####=========================================================================###
	####=========================================================================###

	# OUTPUT TEXT

	####=========================================================================###
	####=========================================================================###

	output$phenotype_table_text_line1<-renderUI({HTML("<b>Note:</b> This table may include the same individual more than once if reported across multiple studies.")})
	output$age_text_line1<-renderUI({HTML(paste("<b>Note:</b> The population ALS penetrance for ", selected_variant$selected_variant, " is ", lit_review_pen$lit_review_pen_ALS_penetrance_literal[lit_review_pen$HGVS==selected_variant$selected_variant],". These images aim to be a useful comparison between variants but may not accurately represent age related penetrance",sep=""))})
	output$age_text_line2<-renderUI({HTML(paste("<b>Note:</b> ",input$gene_variant_browser," and all carrier distributions exclude ", selected_variant$selected_variant,sep=""))})
	output$age_text_line3<-renderUI({age_text_line_3(selected_variant$selected_variant,all_ages_df,age_p_cutoff)})
	output$regions_analysis_text_line1<-renderUI({HTML("<b>Note:</b> if a gene is absent from the plot this indicates that there is either not an available population study for the given region or there are no Pathogenic or Likely Pathogenic variants identfified in that gene.")})
	output$regions_analysis_text_line2<-renderUI({HTML("<b>Note:</b> The 'Familial Proportion' slider is present as the proportion of ALS cases with a positive family history differs across studies. Our analysis is performed separately on familial and sporadic cases.")})
	output$home_header_text<-renderUI({HTML("<b>Welcome to the journALS Data Browser</b>")})
	output$upload_header_text<-renderUI({HTML("<b>Annotate your Data</b>")})
	output$upload_body_text<-renderUI({HTML("<div style='color:grey'>Here you can annotate your own data. The upload file should be a .txt file in chr:pos:ref:alt format (see example download below). <br><br>Data should be aligned to the hg19 version of the human reference genome. Please note, we perform no realignment or <a href='https://genome.sph.umich.edu/wiki/Variant_Normalization'>normalization</a> of your data. If your identifier does not match ours exactly, your annotation may be missed. Our data are split, left-align and trimmed using <a href='https://genome.sph.umich.edu/wiki/Vt#Normalization'>vt normalize</a>. A file explaining the columns in the output data can be downloaded below.")})
	output$home_body_text<-renderUI({HTML("<div style='color:grey'>This data browser represents an amalgamation and uniform analysis of 30 years of genetic screening in amyotrophic lateral sclerosis (ALS) and frontotemporal dementia (FTD). Data is included from 1,028 studies conducted over the past 30 years. 3,112 variants identified in 356 genes are analysed in a uniform manner assessing for penetrance, prevalence and pathogenicity. Please visit the gene and variant browsers for further information.")})
	output$pedigree_text_line0<-renderUI({HTML(paste("<b>Variant: </b>",pedigree_call_variants(input$pedigree_selection,selected_variant$selected_variant),sep=""))})
	output$pedigree_text_line1<-renderUI({HTML("<b>Note:</b> Age indicates age of onset for cases or last known age for unaffected family members ")})
	output$pedigree_text_line2<-renderUI({pedigree_line_2_text(input$pedigree_selection,selected_variant$selected_variant)})
	output$pedigree_text_line3<-renderUI({pedigree_line_3_text(input$pedigree_selection,lit_review_pen,selected_variant$selected_variant)})
	output$heterogeneity_caveat1<-renderUI({HTML("<b>Note:</b> Significant heterogeneity could result from real geographic heterogeneity or may indicate a difference in reporting e.g. for common variants.")})
	output$heterogeneity_caveat2<-renderUI({HTML("<b>Note:</b> Absence of a variant may indicate real absence or may indicate lack of reporting e.g. for synonymous, common or intronic variants.")})
	output$comparison_plot_caveat1<-renderUI({HTML("<b>Note:</b> For full methods on calculating allele frequencies and penetrance estimates please refer to XXX.")})	
	output$comparison_plot_caveat2<-renderUI({HTML("<b>Note:</b> 95% confidence intervals are displayed for penetrance estimates of highlighted variants.")})	
	output$about.FAQ.header<-renderUI({HTML("<b>FAQ</b>")})
	output$about.FAQ.citation.header<-renderUI({HTML("<div style='color:nicedarkblue'>Is there more information available?")})
	output$about.FAQ.citation.body<-renderUI({HTML("<div style='color:grey'>More information can be found in our <a href='https://google.com'>publication</a>. <br><br> Please cite as: XXX<br><br>")})
	output$about.FAQ.gnomAD_AFs.header<-renderUI({HTML("<div style='color:nicedarkblue'>Why do allele frequencies differ from the gnomAD database?")})
	output$about.FAQ.gnomAD_AFs.body<-renderUI({HTML("<div style='color:grey'>The gnomAD database contains ALS samples. We use the gnomAD non-neuro subset, a collection of 104,068 exomes and 10,636 genomes.<br><br>")})
	output$about.FAQ.resources.header<-renderUI({HTML("<div style='color:nicedarkblue'>What resources are used in this database?")})
	output$about.FAQ.resources.body<-renderUI({HTML("<div style='color:grey'>This work wouldn't have been possible without relying on several excellently collected, generated and curated resources; including: <br><a href='https://gnomad.broadinstitute.org/'>gnomAD</a><br><a href='http://databrowser.projectmine.com/'>Project MinE</a><br><a href='http://als.umassmed.edu/#:~:text=Welcome%20to%20the%20ALS%20Variant,genes%20and%20reduce%20false%2Dpositives./'>ALS Variant Server</a><br><a href='http://alsdb.org/'>ALSdb</a><br><a href='https://sites.google.com/site/jpopgen/dbNSFP'>dbNSFP</a><br><br>")})
	output$summary.header.text<-renderUI({HTML("<b>Results Overview</b>")})
	output$summary.overview.text1<-renderUI({HTML("<div style='color:grey'>Our analysis identifies XXX pathogenic or likely pathogenic variants in XXX genes. Please see below for a a summary of the individual properties of these genes. While non benign variants in these genes are present in 67% of familial ALS and 51% of familial FTD, this is reduced to 39% and 38% respectively when considering strictly pathogenic variants.")})
	output$summary.gene.text1<-renderUI({HTML(gene_text(input$summary_gene_selection))})


	###################################
	###     VARIANT TABLES          ###
	###################################
	output$variant_table <- DT::renderDataTable({
		columns<-input$visible_columns
		dataset<-lit_review_pen[lit_review_pen$gene==input$gene,]
		if ("Literature" %in% input$dataset_type){
			dataset<-dataset[dataset$all_carriers_count>0 & !is.na(dataset$all_carriers_count),]
		}
		if ("ALSdb" %in% input$dataset_type){
			dataset<-dataset[dataset$ALSdb_AF>0 & !is.na(dataset$ALSdb_AF),]
		}		
		if ("ALSVS_FAMILIAL" %in% input$dataset_type){
			dataset<-dataset[dataset$ALSVS_FALS_AF>0 & !is.na(dataset$ALSVS_FALS_AF),]
		}
		if ("ProjectMinE" %in% input$dataset_type){
			dataset<-dataset[dataset$ProjectMinE_AF_cases>0 & !is.na(dataset$ProjectMinE_AF_cases),]
		}
		dataset=data.frame(dataset[dataset$impact %in% input$variant_type,])[input$visible_columns]
		DT::datatable(
			#colnames=c("HGVS"="HGVS2"),
			data=dataset,
			extensions = 'FixedColumns',
			options = list(
				scrollX = TRUE,
				paging=FALSE,
				scrollY = 300,
				pageLength = 50,
				deferRender = TRUE,
				autoWidth = TRUE,
				columnDefs = list(list(width = '20px', targets = columns))
				),
			rownames= FALSE,escape=FALSE
			)
		})
	#
	output$variant_table_single_variant <- DT::renderDataTable({
		columns<-input$visible_columns_single_variant
		DT::datatable(
			#colnames=c("HGVS"="HGVS2","gene"="gene2"),
			data=data.frame(lit_review_pen[lit_review_pen$HGVS==selected_variant$selected_variant,])[input$visible_columns_single_variant],
			extensions = 'FixedColumns',
			options = list(
				scrollX = TRUE,
				paging=FALSE,
				scrollY = 90,
				pageLength = 50,
				deferRender = TRUE,
				autoWidth = TRUE,
				columnDefs = list(list(width = '20px', targets = columns))
				),
			rownames= FALSE,escape=FALSE
			)
		})
	#
	output$variant_table_single_variant_phenotype <- DT::renderDataTable({
		columns <- pheno_frame_colnames
		pheno_frame_use <- pheno_frame[pheno_frame$HGVS==selected_variant$selected_variant,][input$visible_columns_single_variant_phenotype]
		DT::datatable(
			data=pheno_frame_use,
			extensions = 'FixedColumns',
			options = list(
				scrollX = TRUE,
				paging=FALSE,
				scrollY = 300,
				pageLength = 50,
				deferRender = TRUE,
				autoWidth = TRUE,
				columnDefs = list(list(width = '20px', targets = columns))
				),
			rownames= FALSE,escape=FALSE
			)
		})
	#

	output$region_browser_individuals_table <- DT::renderDataTable({
		if (input$selection_region_browser=="Country"){
			if (input$region_region_browser=="Japan"){
				pheno_frame_use<-pheno_frame[grep("Japan|Kii_Peninsula", pheno_frame$Nationality), ]
			} 
			else if (input$region_region_browser=="UK"){
				pheno_frame_use<-pheno_frame[(grep("UK|England|Scotland", pheno_frame$Nationality)),]
			} 
			else if (input$region_region_browser=="Italy"){
				pheno_frame_use<-pheno_frame[grep("Italy|Sicily|Sardinia", pheno_frame$Nationality), ]
			} 
			else {
				pheno_frame_use<-pheno_frame[grep(input$region_region_browser, pheno_frame$Nationality), ]		
			}
		} 
		else if (input$selection_region_browser=="Continent"){
			pheno_frame_use<-pheno_frame[pheno_frame$Continent==populations[match(input$region_region_browser,populations_literal)] & !is.na(pheno_frame$Continent), ]
		}
		else if (input$selection_region_browser=="Global"){
			pheno_frame_use<-pheno_frame
		}
		columns <- pheno_frame_colnames
		DT::datatable(
			data=pheno_frame_use,
			extensions = 'FixedColumns',
			options = list(
				scrollX = TRUE,
				paging=FALSE,
				scrollY = 300,
				pageLength = 50,
				deferRender = TRUE,
				autoWidth = TRUE,
				columnDefs = list(list(width = '20px', targets = columns))
				),
			rownames= FALSE,escape=FALSE
			)
		})
	#

	output$region_browser_pop_studies_table <- DT::renderDataTable({
		if (input$selection_region_browser=="Global"){
			dataset=population_studies
		}
		else if (input$selection_region_browser=="Continent"){
			dataset=population_studies[population_studies$continent==input$region_region_browser,]
		}
		else if (input$selection_region_browser=="Country"){
			dataset=population_studies[population_studies$Country==input$region_region_browser,]
		}
		# Remove unneccesary columns
		columns<- colnames(dataset)
		DT::datatable(
			data=dataset,
			extensions = 'FixedColumns',
			options = list(
				scrollX = TRUE,
				paging=FALSE,
				scrollY = 300,
				pageLength = 50,
				deferRender = TRUE,
				autoWidth = TRUE,
				columnDefs = list(list(width = '20px', targets = columns),list(className = 'dt-center', targets = 0:6))
				),
			rownames= FALSE,escape=FALSE
			)
		})



	####=========================================================================###
	####=========================================================================###

	# OBSERVATIONS

	####=========================================================================###
	####=========================================================================###
	# These are dropdown options etc that get updated based on other selections
	observe({
		updateSelectizeInput(
			session,
			"gene",
			choices=unique(lit_review_pen$gene[!is.na(lit_review_pen$gene)]),
			server=TRUE,
			selected="TARDBP"
			)
		})
	observe({
		updateSelectizeInput(
			session,
			"summary_gene_selection",
			#choices=gsub("_post_analysis.tsv","",list.files('Datafiles/Genes')),
			choices=sort(p_or_lp_genes),
			server=TRUE,
			selected="TARDBP"
			)
		})
	observe({
		updateSelectizeInput(
			session,
			"region_region_browser",
			choices=region_region_choices(input$selection_region_browser),
			server=TRUE#,
			#selected="Ireland"
			)
		})
	observe({
		updateSelectizeInput(
			session,
			"variant_highlighter",
			choices=variant_highlighter_choices(lit_review_pen,as.numeric(input$gnomAD_selection),input$gene,as.numeric(input$dataset_selection),as.numeric(input$colour_fill),input$gene_plot_pathogenicity_filter),
			server=TRUE,
			selected="Ireland"
			)
		})

	observe({
		updateSelectizeInput(
			session,
			"transcript_selection",
			choices=unique(gene_coordinates$name[gene_coordinates$name2==input$gene]),
			server=TRUE
			)
		})
	observe({
		updateSelectizeInput(
			session,
			"gene_variant_browser",
			choices=unique(lit_review_pen$gene[!is.na(lit_review_pen$gene)]),
			server=TRUE,
			selected="TARDBP"
			)
		})



	observe({
		variant_drop_down=lit_review_pen[lit_review_pen$gene==input$gene_variant_browser,][order(pos),]
		if ("Literature" %in% input$filter_variant_browser){
			variant_drop_down<-variant_drop_down[variant_drop_down$all_carriers_count>0 & !is.na(variant_drop_down$all_carriers_count),]
		}
		if ("ALSdb" %in% input$filter_variant_browser){
			variant_drop_down<-variant_drop_down[variant_drop_down$ALSdb_AF>0 & !is.na(variant_drop_down$ALSdb_AF),]
		}		
		if ("ALSVS_FAMILIAL" %in% input$filter_variant_browser){
			variant_drop_down<-variant_drop_down[variant_drop_down$ALSVS_FALS_AF>0 & !is.na(variant_drop_down$ALSVS_FALS_AF),]
		}
		if ("ProjectMinE" %in% input$filter_variant_browser){
			variant_drop_down<-variant_drop_down[variant_drop_down$ProjectMinE_AF_cases>0 & !is.na(variant_drop_down$ProjectMinE_AF_cases),]
		}
		for (entry in c("Pathogenic","Likely Pathogenic","VUS","Likely Benign","Benign")){
			if (!(entry %in% input$pathogenicity_filter_variant_browser)){
				variant_drop_down<-variant_drop_down[variant_drop_down$acmg_literal != entry,]
			}
		}
		variant_drop_down$all_carriers_count[is.na(variant_drop_down$all_carriers_count)] <- 0
		variant_drop_down=variant_drop_down[variant_drop_down$gene==input$gene_variant_browser,]
		variant_drop_down_selection=paste0(variant_drop_down$HGVS, " (", variant_drop_down$all_carriers_count, ")")
		# Preselect the variant with the most carriers in the literature
		pre_selection=paste0(head(variant_drop_down$HGVS[variant_drop_down$all_carriers_count==max(variant_drop_down$all_carriers_coun)],1), " (", head(variant_drop_down$all_carriers_count[variant_drop_down$all_carriers_count==max(variant_drop_down$all_carriers_count)],1), ")")
		#variant_drop_down=variant_drop_down$HGVS[variant_drop_down$gene==input$gene_variant_browser]
		updateSelectizeInput(
			session,
			"variant_variant_browser",
			choices=variant_drop_down_selection,
			# What follows is some javascript to left align variant and right align count
			options = list(render = I(
				'{
					option: function(item, escape) {
						const x = item.value.split(" ");
						return `<p style=\"text-align:left;\">
						${x[0]}
						<span style=\"float:right;\">
						${x[1]}
						</span>
						</p>`
					}
					}')),
			server=TRUE,
			selected=pre_selection
			)
		})

	#######
	#######
	# Tidy the variant selection i.e. remove: (number)
	# Both of these are needed for this to work
	selected_variant <- reactiveValues()
	observe({
		print(gsub(" .*","",input$variant_variant_browser))
		selected_variant$selected_variant<-gsub(" .*","",input$variant_variant_browser)
		})
	#######
	#######


	observe({
		updateSelectizeInput(
			session,
			"pedigree_selection",
			choices=unique(separate_rows(lit_review_pen[lit_review_pen$HGVS==selected_variant$selected_variant & is.na(lit_review_pen$HGVS)==FALSE,],all_carriers_pedigree,sep="\\|",convert=T)$all_carriers_pedigree)[is.na(unique(separate_rows(lit_review_pen[lit_review_pen$HGVS==selected_variant$selected_variant & is.na(lit_review_pen$HGVS)==FALSE,],all_carriers_pedigree,sep="\\|",convert=T)$all_carriers_pedigree))==FALSE],server=TRUE)
		})

	####=========================================================================###
	####=========================================================================###

	# DOWNLOADS

	####=========================================================================###
	####=========================================================================###

	# TAB = Gene browser / comparison plots 

	output$download.comparison_plots <- downloadHandler(
		filename = function(){
			paste(input$gene,"_features.png",sep='')
			},
			contentType = 'image/png',
			content=function(file){
				png(file,width = 2800, height = 2000, units = 'px',res=300)
				par(mfrow=c(1,1),
					mar=c(6,6,2,2),
					las=0,
					mgp = c(4, 1, 0))
				render.comparison_plots(lit_review_pen,input$variant_highlighter,as.numeric(input$gnomAD_selection),as.numeric(input$dataset_selection),as.numeric(input$colour_fill),input$gene,input$freq,as.numeric(input$regression_line),input$gene_plot_pathogenicity_filter)				
				
				dev.off()
				})
	# TAB = Regions Plot

	output$download.regions_analysis_plot <- downloadHandler(
		filename = function(){
			"cases.explained.png"
			},
			contentType = 'image/png',
			content=function(file){
				png(file,width = 1960, height = 2520, units = 'px',res=300)
				par(mfrow=c(1,1),
					mar=c(6,6,2,2),
					las=0,
					mgp = c(4, 1, 0))
				render.regions_analysis_plot(input$selection_region_browser,input$region_region_browser,input$regions_pathogenicity_selection,input$regions_phenotype_selection,input$regions_history_selection,input$fam_prop_slider,lit_review_pen)
				dev.off()
				})

	# TAB = Variant Browser / Geographic Heterogeneity

	output$download.heterogeneity_plots <- downloadHandler(
		filename = function(){
			paste(selected_variant$selected_variant,input$heterogeneity_phenotype_selection,input$heterogeneity_family_selection,"heterogeneity.png",sep='_')
			},
			contentType = 'image/png',
			content=function(file){
				png(file,width = 2240, height = 1600, units = 'px',res=300)
				par(mfrow=c(1,1),
					mar=c(6,6,2,2),
					las=0,
					mgp = c(4, 1, 0))
				render.heterogeneity_plots(input$gene_variant_browser,selected_variant$selected_variant,lit_review_pen,as.numeric(input$heterogeneity_selection),input$heterogeneity_phenotype_selection,input$heterogeneity_family_selection)
				dev.off()
				})



	# TAB = Variant Browser / Pedigrees

	output$download.pedigree_plot <- downloadHandler(
		filename = { 
			paste(input$pedigree_selection,".ped.png",sep="")
			},
			contentType='image/png',
			content=function(file){
				png(file,
					width=24,
					height=8,
					units='in',
					res=300)
				render.pedigree_plot(input$pedigree_selection)
				dev.off()
			}
			)

	output$downloadPedFile <- downloadHandler(
		filename = function(){ 
			paste(input$pedigree_selection,".ped.tsv",sep="")
			},
			contentType='text/tsv',
			content=function(file){
				write.table(read.csv(paste("pedigrees/",input$pedigree_selection,".ped",sep=""),header=T,sep='\t'),file,row.names=F,quote=F)
			}
			)

	# TAB = Upload 

	output$downloadUploadResults <- downloadHandler(
		filename = function(){ 
			"journALS_results.tsv"
			},
			contentType='text/tsv',
			content=function(file){
				write.table(lit_review_pen[!is.na(lit_review_pen$identifier) & lit_review_pen$identifier %in% read.table(input$upload_file$datapath,header=FALSE)$V1,],file,sep='\t',col.names=T,row.names=F,quote=F)
			}
			)	

	output$downloadExampleTxt <- downloadHandler(
		filename="example.txt",
		contentType='text/tsv',
		content=function(file){
			write.table(read.csv("Datafiles/example.txt",header=F,sep='\t'),file,row.names=F,quote=F,sep='\t')
		}
		) 

	output$downloadHeaderExplanation <- downloadHandler(
		filename="HeaderExplanation.txt",
		contentType='text/tsv',
		content=function(file){
			write.table(read.csv("Datafiles/rev14_header_explanation.tsv",header=T,sep='\t'),file,row.names=F,quote=F,sep='\t')
		}
		)   

}
