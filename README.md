# journALS
### To Start:

**Step1** :    		`git clone https://github.com/dohertymark/journALS_Databrowser.git` <br />
**Step2** : 		Navigate to your local journALS directory       	<br />
**Step3** :		`bash Data_processing/correct_directories.sh pwd`  	<br /><br /> <br />
### TO RUN journALS  <br />
**Step1** : 		`R` <br />
**Step2** : 		install.packages("shinyjs") <br />
**Step3** : 		install.packages("shinyWidgets") <br />
**Step2** : 		`library(shiny)`  <br />
**Step3** : 		`runApp("App/")`  <br />
			(may need to install packages) <br />
			databrowser should open in your browser <br /> <br /><br />
### TO DO DATA PREPROCESSING YOURSELF  <br />
**Step1** : 		from pg4 download : /st1/hdd/ctg/Mark/variant_analysis/data_files/ALS_whole_data_rev14.gz to journALS/Datafiles folder 	<br />
**Step2** :		`Rscript Data_processing/journALS_data_analysis_and_filtering.R` 						<br />
			(may need to install packages) <br />
			(this script will generate the file 'journALS_post_analysis.tsv.gz' in the Datafiles folder)
