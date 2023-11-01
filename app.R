## Omixon CLI Explore Genotyping Portal - 
# v0.2
# By: Livia Tran 
# March 30, 2021

library(shiny)
library(shinyjs)
library(shinybusy)
library(disambiguateR)
library(jsonlite)
library(httr)
library(R.utils)
library(data.table)
library(doMC)
library(foreach)

options(shiny.maxRequestSize=100000*1024^2)

#current AWS instance max core = 4
doMC::registerDoMC(3)

#create log file if not already in existence 
if(!file.exists("/var/omixon_output/geno_logfile.log")){
    file.create("/var/omixon_output/geno_logfile.log")
}

#batchDF - split approved subejcts into batches - v 0.1 LT 12/16/20
batchSplit<-function(df, batchNum){
    batchDF<-NULL
    
    #based number of loops on rounded up number of rows in dataframe divided by desired number of batch splits
    for(i in 1:ceiling(nrow(df)/batchNum)){
        
        #take 1:batchNum subjects if # rows in dataframe is greater than or equal to the batch number
        if(nrow(df) >= batchNum){  
            batchDF[[i]]<-df[c(1:batchNum),]
        }
        
        #if remaining rows are less than the batch number, take remaining rows instead 
        else{
            batchDF[[i]]<-df[1:nrow(df),]
        }
        
        #remove rows previously subset
        df<-df[-c(1:batchNum),]
    }
    return(batchDF)
}

#writeLog
#singularizes parameters for log message appending 
writeLog<-function(mess){
    write.table(paste(Sys.time(), sep=" ", mess), file="/var/omixon_output/geno_logfile.log", append = T, row.names=F, col.names = F)
}

userValidity<-function(useremail, project, originid){
    
    h <- c("hibagger@hlacovid19.org", "712ZTRhZBKPzZmmGxD56", "application/json")
    
    names(h) <- c("X-User-Email", "X-User-Token", "Content-type")
    
    content <-fromJSON(content(POST("https://database-hlacovid19.org/query/hibag_preflight",
                                    body = paste('{"email":"',useremail,'","project_name":"',project,'","origin_identifiers":', toJSON(originid),'}', sep = ""),
                                    add_headers(.headers = h), encode='json'), "text"))
    
    if(content$user_approved == TRUE & content$project_found == TRUE){
        projValidity<-TRUE
        
        writeLog(paste0("The e-mail ", useremail, " and ", project, " project are valid."))
        
    }
    
    if(content$user_approved == TRUE & content$project_found == FALSE){
        projValidity<-paste0("The ",  project, " project has not been registered in the COVID-19|HLA Database. A registered project name is required for imputation.")
        
        writeLog(paste0("The e-mail ", useremail, " is valid, but the ", project, " project is not valid."))
    }
    
    
    if(content$user_approved == FALSE){
        projValidity<-paste0("The ",  useremail, " email address has not been registered in the COVID-19|HLA Database. A registered e-mail is required for imputation.")
        
        writeLog(paste0("The e-mail ", useremail, " is not valid."))
    }
    
    return(projValidity)
}


#set filename to most recently generated .csv file with FASTQ file name
CLIEx2GLString <- function(FASTQfile){
    
    if(length(Filter(function(x) grepl(gsub(".fastq.gz", "", FASTQfile), x),  list.files(path="/var/omixon_output",pattern=".csv")))==0){
        return("FALSE")
    }
    
    else{
        CLIgeno<-read.table(file=row.names(file.info(paste("/var/omixon_output/", Filter(function(x) grepl(gsub(".fastq.gz", "", FASTQfile), x),  list.files(path="/var/omixon_output", pattern=".csv")), sep=""))[file.info(paste("/var/omixon_output/", Filter(function(x) grepl(gsub(".fastq.gz", "", FASTQfile), x),  list.files(path="/var/omixon_output",pattern=".csv")), sep=""))$mtime == max(file.info(paste("/var/omixon_output/", Filter(function(x) grepl(gsub(".fastq.gz", "", FASTQfile), x),  list.files(path="/var/omixon_output",pattern=".csv")), sep=""))$mtime),]),header = TRUE,sep = ";",quote = "",as.is = TRUE,colClasses = "character")
        
        ## Remove NA columns and the trailing blank column
        toRemove <- colnames(CLIgeno)[colSums(is.na(CLIgeno)) > 0]
        CLIgeno <- CLIgeno[ , !(names(CLIgeno) %in% append("X",toRemove,length(toRemove)))]
        
        glString <- vector()
        ## Check to make sure everything matches 
        for(i in 2:ncol(CLIgeno)) {
            if(lengths(regmatches(CLIgeno[[i]][1], gregexpr("/", CLIgeno[[i]][1]))) != lengths(regmatches(CLIgeno[[i]][2], gregexpr("/", CLIgeno[[i]][2])))) {
                cat("Error: mismatch between allele string lengths in", colnames(CLIgeno)[i],".\n")
            }
            ## An error does not currently halt execution 
            # else {
            glString <- append(glString, paste(paste(unlist(strsplit(CLIgeno[[i]][1],"/",fixed=TRUE)),unlist(strsplit(CLIgeno[[i]][2],"/",fixed = TRUE)),sep="+"),collapse="|"),length(glString))
            # }
        }
        return(paste(glString,collapse="^"))
    }
}


### Extract Read Counts Function - v 1.1 SJM 12/01/2020
##
##  Identifies the read counts for each HLA locus in an Omixon CLI Explore .qcresult file (from the unzipped "HLA Typing Result" (.htr) archive)
##  Returns a 1-row data-frame of read-counts for each locus

extract_read_counts <- function(qc_file){
    
    
    qc_text <- readLines(qc_file,warn=FALSE)
    
    qc_blocks_start <- unlist(gregexpr(pattern="READ_COUNT",qc_text))[1:length(unlist(gregexpr(pattern="READ_COUNT",qc_text)))]
    qc_blocks_end <- unlist(gregexpr(pattern='(intergenic ambiguity)',qc_text))[1:length(unlist(gregexpr(pattern='(intergenic ambiguity)',qc_text)))]
    
    qc_table <- data.frame(locus=character(0),read_count=integer(0))
    
    for(i in 1:length(qc_blocks_start)){
        qc_block <- substr(qc_text,qc_blocks_start[i],qc_blocks_end[i])
        val_start <- unlist(gregexpr(pattern = "value",qc_block))
        val_end <- unlist(gregexpr(pattern = "CROSSMAPPING",qc_block))
        locus_start <- unlist(gregexpr(pattern = "locus",qc_block))[2]
        qc_table[i,1] <- unlist(substr(qc_block,locus_start+8,val_start-4))
        qc_table[i,2] <- as.numeric(substr(qc_block,val_start+7,val_end-4))
    }
    qc_final <- as.data.frame(t(qc_table[,2]))
    colnames(qc_final) <- t(qc_table[,1])
    rownames(qc_final) <- "read count"
    qc_final
}


ui = fluidPage(
    useShinyjs(),
    add_busy_spinner(spin = "fading-circle"),
    
    titlePanel("Omixon CLI Explore Genotyping Portal"),
    
    sidebarPanel(
        fileInput("userinfo", "User Information .csv",
                  accept = ".csv", multiple=F),
        
        fileInput("fastqfiles", "FASTQ.GZ files",
                  accept = ".gz", multiple=T),
        
        checkboxInput("filterCheck", label="These FASTQ.GZ files have been filtered so that they include only HLA-gene reads."),
        
        
        actionButton("run", "Genotype!"),
        HTML("<br>"),
        HTML("<br>"),
        actionButton("reset", "Clear")),
    mainPanel(
        tabsetPanel(type="tabs",
                    
                    tabPanel("Main", 
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             textOutput("text"), 
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),
                             HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"),
                             
                             fluidRow(column(12,img(src=base64enc::dataURI(file="logo.png", mime="image/png"), align="right", width="300"), height = 200)),
                             
                    ),
                    
                    tabPanel("About", uiOutput("aboutText"),     
                             HTML("<br>"), HTML("<br>"),HTML("<br>"),HTML("<br>"),HTML("<br>"), HTML("<br>"), HTML("<br>"), HTML("<br>"),
                             fluidRow(column(12,img(src=base64enc::dataURI(file="logo.png", mime="image/png"), align="right", width="300"), height = 200)),)),
    )
    
)

server <- function(input, output, session) {
    
    failCount<-0
    
    userinfoRun<-reactiveValues(userinfo=NULL)
    fastqfilesRun<-reactiveValues(fastqfiles=NULL)
    
    disable("run")
    
    isolate(observeEvent(input$userinfo, {
        userinfoRun$userinfo<-input$userinfo
        
    }))
    
    isolate(observeEvent(input$fastqfiles, {
        fastqfilesRun$fastqfiles<-input$fastqfiles
        
    }))
    
    #enable run button once FASTQfiles AND user info.csv are uploaded
    isolate(observe({
        if(length(userinfoRun$userinfo!=0) & length(fastqfilesRun$fastqfiles)!=0 & input$filterCheck ==T){
            enable("run")
        }
        else{
            disable("run")
        }
        
    }))
    
    
    isolate(observeEvent(input$run, {  
        
        withCallingHandlers({
            shinyjs::html("text", "")
            tempPath <- dirname(fastqfilesRun$fastqfiles$datapath)[1]
            
            file.rename(from=fastqfilesRun$fastqfiles$datapath,to=paste(tempPath,fastqfilesRun$fastqfiles$name,sep="/"))
            
            writeLog(paste(input$userinfo$name, "uploaded for user info .csv"))
            
            userInfo<-read.csv(userinfoRun$userinfo$datapath, header=F, sep=",", stringsAsFactors = F, na.strings=c("", NA))
            
            useremail<-userInfo[1,1]
            project<-userInfo[1,2]
            
            originid<-userInfo[2:nrow(userInfo),1]
            
            #put subject:fastq files into a different variable 
            subjectRuns<-userInfo[2:nrow(userInfo),]
            
            projectValidity<-userValidity(useremail, project, originid)
            
            if(projectValidity != "TRUE"){
                message(projectValidity)
            }
            
            else{
                
                message(paste0("The e-mail ", useremail, " and ", project, " project are valid."))
                
                #static configuration files and base command to run CLI Explore
                config <- "/home/ubuntu/cli_ex/user2.conf"
                conf <- paste("-conf", config, sep=" ")
                baseCommand <- "java -jar /home/ubuntu/cli-explore-app-0.1.0-SNAPSHOT-pkg.jar"
                
                failCountTotal<-qcpath<-readMetricsFASTQ<-FASTQ1Uploaded<-subject<-NULL
                
                writeLog(paste(nrow(fastqfilesRun$fastqfiles), "FASTQ files uploaded."))
                
                failSubj<-apprSubj<-GLstringFASTQ1<-NULL
                
                for(i in 1:nrow(subjectRuns)){
                    if((!subjectRuns[i,2] %in% fastqfilesRun$fastqfiles$name)==TRUE | (!is.na(subjectRuns[i,3]) & (!subjectRuns[i,3]  %in% fastqfilesRun$fastqfiles$name)==TRUE)){
                        writeLog(paste("The FASTQ file", fastqfilesRun$fastqfiles$name[[i]], "uploaded does not match any of the FASTQ files provided in the .csv file."))
                        message(paste("The FASTQ file", fastqfilesRun$fastqfiles$name[[i]], "uploaded does not match any of the FASTQ files provided in the User Information .csv file. Please ensure the FASTQ filenames match the ones provided in the User Information .csv file. Moving on to next uploaded FASTQ file."))
                        
                        failSubj[[i]]<-subjectRuns[i,]
                        
                    }
                    
                    else{
                        apprSubj[[i]]<- subjectRuns[i,]
                        
                    }
                } 
                
                failSubjdf<-rbindlist(failSubj)
                
                #bind all approved subjects into one data frame
                apprSubjdf<-data.frame(rbindlist(apprSubj), stringsAsFactors = F)
                
                
                apprSubjSplit<-batchSplit(apprSubjdf, 2)
                
                out<-paste("-out /var/omixon_output")
                
                CLIresults <- sapply(apprSubjSplit, function(x) NULL)
                
                for(i in 1:length(apprSubjSplit)){
                    
                    CLIresults[[i]]<-sapply(1:nrow(apprSubjSplit[[i]]), function(x) NULL)
                    
                    panelMessages<-foreach(a = 1:nrow(apprSubjSplit[[i]])) %dopar% {
                        
                        
                        #if NA is present in third column, singular FASTQ file 
                        if(is.na(apprSubjSplit[[i]][a,3])){
                            
                            
                            writeLog(paste("Genotyping HLA alleles for subject ", apprSubjSplit[[i]][a,1], "....", sep=""))
                            
                            panelMessages<-paste("Genotyping HLA alleles for subject ", apprSubjSplit[[i]][a,1], "....", sep="")
                            
                            CLIresults[[i]][[a]]<-system(paste(baseCommand, paste("-in1", paste(tempPath, apprSubjSplit[[i]][a,2], sep="/")) ,conf, out, sep=" "), intern=T)
                            
                            if(length(CLIresults[[i]][[a]]) < 600){
                                writeLog(paste("Genotyping for subject ", apprSubjSplit[[i]][a,1], " failed due to corrupt FASTQ file(s)."))
                                list(1, append(panelMessages, paste("The specified FASTQ file(s) for ", apprSubjSplit[[i]][a,1], " were not found or are corrupt. Genotyping halted.")))
                            }
                            
                            else{
                                writeLog(paste("Genotyping for subject ", apprSubjSplit[[i]][a,1], " completed."))
                                list(0, append(panelMessages, paste("Genotyping for subject ", apprSubjSplit[[i]][a,1], " completed.")))
                            }
                        }
                        
                        #for subjects with paired FASTQ files
                        else{
                            
                            writeLog(paste("Genotyping HLA alleles for subject ", apprSubjSplit[[i]][a,1], "....", sep=""))
                            panelMessages<-paste("Genotyping HLA alleles for subject ", apprSubjSplit[[i]][a,1], "....", sep="")
                            
                            CLIresults[[i]][[a]]<- system(paste(baseCommand,paste("-in1", paste(tempPath, apprSubjSplit[[i]][a,2], sep="/")), paste("-in2", paste(tempPath, apprSubjSplit[[i]][a,3], sep="/")), conf, out, sep=" "), intern=T)
                            
                            if(length(CLIresults[[i]][[a]]) < 600){
                                writeLog(paste("Genotyping for subject ", apprSubjSplit[[i]][a,1], " failed due to corrupt FASTQ file(s)."))
                                
                                list(1, append(panelMessages, paste("The specified FASTQ file(s) for ", apprSubjSplit[[i]][a,1], " were not found or are corrupt. Genotyping halted.")))
                            }
                            
                            else{
                                writeLog(paste("Genotyping for subject ", apprSubjSplit[[i]][a,1], " completed."))
                                list(0, append(panelMessages, paste("Genotyping for subject ", apprSubjSplit[[i]][a,1], " completed.")))
                            }
                        }
                    }
                    message(paste(sapply(panelMessages, "[[", 2), collapse = "\n"))
                    
                    failCountTotal[[i]]<-as.numeric(sapply(panelMessages, "[[", 1))
                }
                
                for(j in 1:nrow(apprSubjdf)){
                    
                    GLstringFASTQ1[[j]]<-CLIEx2GLString(apprSubjdf[j,2])
                    
                    if(GLstringFASTQ1[[j]]!="FALSE"){
                        
                        write.table(canonicalize(GLstringFASTQ1[[j]])$glstring, file = paste("/var/omixon_output/", paste(paste(paste(paste("__", gsub(".", "_dot_", gsub("@", "_at_", useremail), fixed=T), sep=""), "__", sep=""), paste(project, "__", sep=""), sep=""), paste(originid[[j]], "__", Sys.Date(), ".txt", sep=""), sep=""), sep=""), row.names = F, col.names = F)
                        
                        writeLog(paste("GL string formatted HLA genotypes have been loaded to the database for subject ", apprSubjdf[j,1], ".", sep=""))
                        
                        message(paste("GL string formatted HLA genotypes have been loaded to the database for subject ", apprSubjdf[j,1], ".", sep=""))
                        
                        #create temp directory to unzip individual .htr files
                        system("mkdir /var/omixon_output/temp")
                        
                        unzip(row.names(file.info(paste("/var/omixon_output/", Filter(function(x) grepl(gsub(".fastq.gz", "", apprSubjdf[j,2]), x),  list.files(path = "/var/omixon_output", pattern=".htr")), sep=""))[file.info(paste("/var/omixon_output/", Filter(function(x) grepl(gsub(".fastq.gz", "", apprSubjdf[j,2]), x),  list.files(path = "/var/omixon_output", pattern=".htr")), sep=""))$mtime == max(file.info(paste("/var/omixon_output/", Filter(function(x) grepl(gsub(".fastq.gz", "", apprSubjdf[j,2]), x),  list.files(path= "/var/omixon_output",pattern=".htr")), sep=""))$mtime),]),exdir = "/var/omixon_output/temp")
                        
                        qcpath[[j]]<-row.names(file.info(paste("/var/omixon_output/temp/", Filter(function(x) grepl(gsub(".fastq.gz", "", apprSubjdf[j,2]), x),  list.files(path = "/var/omixon_output/temp/", pattern=".qcresult")), sep=""))[file.info(paste("/var/omixon_output/temp/", Filter(function(x) grepl(gsub(".fastq.gz", "", apprSubjdf[j,2]), x),  list.files(path = "/var/omixon_output/temp/", pattern=".qcresult")), sep=""))$mtime == max(file.info(paste("/var/omixon_output/temp/", Filter(function(x) grepl(gsub(".fastq.gz", "", apprSubjdf[j,2]), x),  list.files(path= "/var/omixon_output/temp/",pattern=".qcresult")), sep=""))$mtime),])
                        
                        readMetricsFASTQ[[j]]<-extract_read_counts(qcpath[[j]])
                        
                        write.csv(readMetricsFASTQ[[j]], file = paste("/var/omixon_output/" ,paste(paste(paste(paste("__", gsub(".", "_dot_", gsub("@", "_at_", useremail), fixed=T), sep=""), "__", sep=""), paste(project, "__", sep=""), sep=""), paste(originid[[j]], "__read-counts__", Sys.Date(), ".csv", sep=""), sep=""), sep=""))
                        
                        system("rm -r /var/omixon_output/temp")
                        
                        writeLog(paste("Read counts for each HLA locus extracted and loaded to the database for subject ", apprSubjdf[j,1], ".", sep=""))
                        
                        #remove all generated files
                        file.remove(paste("/var/omixon_output/" ,list.files(path = "/var/omixon_output", pattern= gsub(".fastq.gz", "", apprSubjdf[j,2])), sep=""))
                        
                        file.remove(paste(tempPath, list.files(path = tempPath, pattern= gsub(".fastq.gz", "", apprSubjdf[j,2])), sep="/"))
                        
                        message("Read counts for each HLA locus extracted and loaded to the database for subject ", apprSubjdf[j,1], ".", sep="")
                        
                    }}
                message(paste("Successful HLA Genotyping completed for", nrow(subjectRuns) - sum(unlist(failCountTotal)) - nrow(failSubjdf), " subjects."))
            }
        },
        message = function(m) {
            shinyjs::html(id = "text", html = paste("<br>", m$message), add = TRUE)
        }) 
        
    })
    )
    
    isolate(observeEvent(input$reset, {
        shinyjs::reset("userinfo")
        shinyjs::reset("fastqfiles")
        shinyjs::reset("filterCheck")
        output$text<-NULL
        fastqfilesRun$fastqfiles<-NULL
        userinfoRun$userinfo<-NULL
        disable("run")
        
    }))
    
    output$aboutText<-renderUI({
        HTML("<b>Version 0.2</b>
     <br>
     <br>
     <font size='4'><b><u>Authors:</b></font></u><br>
   <br>
   Livia Tran - <i>livia.tran@ucsf.edu</i> 
   <br>
   Steven Mack - <i>steven.mack@ucsf.edu </i>
   <br>
   <br>
   This portal uses <b> Omixon CLI Explore </b> software to perform HLA genotyping using whole genome or whole exome sequencing data for the HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRB1, HLA-DRB3, HLA-DRB4, HLA-DRB5 genes.
   The user manual for Omixon CLI Explore can be found <a href='http://omixon-ucsf.s3.amazonaws.com/Explore_CLI/Omixon_CLI_Explore_User_Manual.pdf' target='_blank'> here.</a>       
   <br>
   <br>
   After 'Genotype!' is clicked, a busy spinner on the top right corner will appear to indicate that the application is running.
   <br>
   <br>
  <i>Note: Uploaded fastq.gz files are deleted when genotyping is complete, and will not be stored in the database. </i>
   <br>
   <br>
   <b>REQUIRED FILES</b>
   <br>
   <br>
  1) <i>A User Information .csv file containing the following:</i>
  <br>
  <br>
   -<u>Submitter's e-mail address</u> : Please provide the email address associated with the submitter's database account. 
   <br>
   -<u>Project Name:</u> Please provide the database Project Name associated with the origin identifiers (defined in the 'origin_identifier' field of the data dictionary) for the individuals being genotyped.
   <br>
   -<u>Origin Identifiers: </u> Please provide the origin identifiers associated with the the individuals being genotyped.
   <br>
   -<u>FASTQ file names: </u> Please provide the name of the FASTQ.GZ file or names of the paired FASTQ.GZ files for each individual being genotyped. It is possible to upload more files than are present in the user .csv file, but that will slow down the file upload time, and delay the intiation of genotyping. Please ensure the number of files uploaded corresponds with the number of subjects in the .csv file. 
   <br>
   <br>
   Below is an example of correct formatting for the User Information .csv file, where subject 1 has one associated FASTQ file, and subject 2 has a set of paired FASTQ files: 
          <br>
          <br>
          <i><b>user e-mail, project name, 
          <br>
           subject 1, FASTQFILE1,
           <br>
           subject 2, FASTQFILE2, FASTQFILE3 </i></b>
           <br>
           <br>
           Note there are two commas in each row to separate 3 total values. The third value will be populated if a subject has a paired FASTQ file; otherwise, the third value should be blank. The third column in the header can contain any data, or can be left blank. This field of the header is not read.
           <br>
           <br>
   2) <i>A <b>.gz (Gnu zipped archive file)</b> of FASTQ files for each subject in the User Information .csv file. Please check that the uploaded FASTQ file names match those that have been entered in the User Information .csv file. All fastq.gz files should be stored in the same directory. </i>     
          <br>
          <br>
          To minimize upload and processing time, all submitted fastq.gz files must be restricted to reads that map to HLA genes. The wgsHLAfiltR (Whole Genome Sequence HLA Filter) R package (available <a href='https://github.com/COVID-HLA/wgsHLAfiltR' target='_blank'> here</a>) generates fastq.gz files of reads that map to the classical HLA loci (HLA-A, -C, -B, -DRB1, -DRB3/4/5, -DQA1, -DQB1, -DPA1 and -DPB1) from paired or individual whole-genome or whole-exome sequencing (WGS/WES) fastg.gz files. We recommend using wgsHLAfiltR to generate HLA gene-only fastq.gz files.
          <br>
          <br>
          ")})
}

# Run the application 
shinyApp(ui = ui, server = server)

