#
# Pathway activity based discriminative analysis
#


#
# Convert gene expression matrix to summarized pathway activity matrix using different dimension reduction or averaging approaches
#
summarized_pathway_activity = function(exprs, gsets=NULL, type="median", minsize = 10)
{
  
  if(is.null(gsets) && is.null(database))
  {
    stop("Either the gsets or the database parameter must be specified!")
  }
  
  exprs <- as.matrix(exprs)

	genenames = rownames(exprs)
	
	# list of pathway expression matrices
	pathmat = matrix(0, nrow=length(gsets), ncol=ncol(exprs))
	rownames(pathmat) = rep("", nrow(pathmat))

	count = 0	
	cat('\n',length(gsets),' pathways read.\n')
	for(j in 1:length(gsets))
	{
	  if(j %% 100 == 0){
	  	cat(paste("Current iteration:",j,"\n"))
	  }
	  
	  gset = gsets[[j]]
	  
	  mapid = match(gset, genenames)
	  
	  notna = which(!is.na(mapid))
	  
	  if(length(notna) <  minsize)
	  	next
	
	  curpathmat = exprs[mapid[notna],]
	  
	  meanpathvec = NULL
	  if(type == "mean") {
	  	meanpathvec = apply(curpathmat, 2, mean)
	  } else if(type=="min"){
	    meanpathvec = apply(curpathmat, 2, min)
	  } else if(type=="max"){
	    meanpathvec = apply(curpathmat, 2, max)	
	  } else if(type=="sd"){
	    meanpathvec = apply(curpathmat, 2, sd)
	  } else if(type=="pca"){
	    rem = which(apply(curpathmat, 1, var)==0)
	    curpathmatfilt = curpathmat
	    if(length(rem))
	    	curpathmatfilt = curpathmat[-rem,]
	    if(length(curpathmatfilt))
	    {
	    	pca    <- prcomp(t(curpathmatfilt), retx=T, scale=T) # scaled pca 
				scores <- pca$x[,1]
				meanpathvec = scores
			} else {
				meanpathvec = rep(0, ncol(exprs))
			}
		} else if(type=="mds") {
		  meanpathvec <- as.vector(cmdscale(dist(t(curpathmat)), k = 1))
	  } else {
	  	meanpathvec = apply(curpathmat, 2, median)
	  } 	  
	  	  
	  count = count + 1
	  pathmat[count,] = meanpathvec
	  rownames(pathmat)[count] = names(gsets)[j]
	}
	
	pathmat = pathmat[1:count,]

  return(pathmat)
}


#
# Load pathway data from the MSIGDB database
#
load_pathway_database = function(database = c("CANONIC","GO","KEGG","REACTOME","PID","BIOCARTA","CHRPOSITION"), vnum="7.4")
{

	msigdb_pathways = NULL
	if(database=="CANONIC")
		msigdb_pathways = sapply(readLines(paste("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",vnum,"/c2.cp.v",vnum,".symbols.gmt",sep="")), function(x) strsplit(x, "\t")[[1]])
	if(database=="GO")
		msigdb_pathways = sapply(readLines(paste("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",vnum,"/c5.all.v",vnum,".symbols.gmt",sep="")), function(x) strsplit(x, "\t")[[1]])
	if(database=="KEGG")
		msigdb_pathways = sapply(readLines(paste("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",vnum,"/c2.cp.kegg.v",vnum,".symbols.gmt",sep="")), function(x) strsplit(x, "\t")[[1]])
	if(database=="REACTOME")
		msigdb_pathways =sapply(readLines(paste("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",vnum,"/c2.cp.reactome.v",vnum,".symbols.gmt",sep="")), function(x) strsplit(x, "\t")[[1]])
	if(database=="PID")
		msigdb_pathways = sapply(readLines(paste("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",vnum,"/c2.cp.pid.v",vnum,".symbols.gmt",sep="")), function(x) strsplit(x, "\t")[[1]])
	if(database=="BIOCARTA")
		msigdb_pathways = sapply(readLines(paste("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",vnum,"/c2.cp.biocarta.v",vnum,".symbols.gmt",sep="")), function(x) strsplit(x, "\t")[[1]])
	if(database=="CHRPOSITION")
		msigdb_pathways = sapply(readLines(paste("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",vnum,"/c1.all.v",vnum,".symbols.gmt",sep="")), function(x) strsplit(x, "\t")[[1]])

	return(msigdb_pathways)
}

#
# Run an example analysis on Parkinson's disease data from the GEO database (requires Affymetrix annotation data in file: HG-U133A.na36.annot.csv)
#
run_example = function()
{

		#
		# Load pathway data
		#
		
		pathdat = load_pathway_database("KEGG")

		# convert to required list format
		pathlst = sapply(pathdat, function(x) x[3:length(x)])
		
		names(pathlst) = sapply(pathdat, function(x) x[1])

		
		#
	 	# Load example Parkinson's disease case/control gene expression dataset from GEO
	 	#
		#    Dataset GSE8397: L. B. Moran et al., Neurogenetics, 2006, SN + frontal gyrus, post mortem,	PD (29), healthy (18)
		#    Array platform: Affymetrix HG-U133A
		#

		# Download the data into the current working directory

		# for Windows - manually via the web-browser using this url:
		#
		# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/matrix/GSE8397-GPL96_series_matrix.txt.gz
		#

		# for Mac/Linux - automatically via R command line
		system('wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/matrix/GSE8397-GPL96_series_matrix.txt.gz')

		# Read the data into R
		morandatgeo = read.table(gzfile("GSE8397-GPL96_series_matrix.txt.gz"), header=T, comment.char="!", sep="\t")

		# Use the labels in the first column as row names
		morandat = morandatgeo[,2:ncol(morandatgeo)]
		rownames(morandat) = morandatgeo[,1]
		
		# Filter out tissue samples which are not from the midbrain / substantia nigra region
		moran_tissues = as.matrix(read.table(gzfile("GSE8397-GPL96_series_matrix.txt.gz"), header=F, nrows=1, skip=36, sep="\t"))
		moran_tissues = moran_tissues[2:length(moran_tissues)]
		
		nigra_ind = grep("substantia nigra",moran_tissues)
		
		moran_outcome = as.matrix(read.table(gzfile("GSE8397-GPL96_series_matrix.txt.gz"), header=F, nrows=1, skip=28, sep="\t"))
		moran_outcome = moran_outcome[2:length(moran_outcome)]
		moran_outcome[grep("control",moran_outcome)] = rep("control",length(grep("control",moran_outcome)))
		moran_outcome[grep("Parkinson",moran_outcome)] = rep("parkinson",length(grep("Parkinson",moran_outcome)))
		
		moranfilt = morandat[,nigra_ind]
		#dim(moranfilt)
		#[1] 22283    39
		
		moran_outcomefilt = moran_outcome[nigra_ind]
		#table(moran_outcomefilt)
		#moran_outcomefilt
		#  control Parkinson 
		#       15        24

		# Unzip annotation file (in Mac/Linux, needs to be done manually in Windows)
		system('unzip HG-U133A.na36.annot.csv.zip')
		
		# read annotations file (ignoring comments)
		annot = read.csv("HG-U133A.na36.annot.csv", comment.char="#")
		head(annot)
		
		# map probes to microarray rownames
		mapids = match(rownames(moranfilt), annot$Probe.Set.ID)
		
		# check if all IDs were mapped successfully
		any(is.na(mapids))
		#[1] FALSE
		# ok, no missing IDs
		
		# extract gene symbols corresponding to microarray Probe IDs (take always the first symbol mapped)
		mapped_symbols = sapply( as.character(annot$Gene.Symbol[mapids]) , function(x) strsplit(x, " /// ")[[1]][1])
		

			
		#
		# Convert expression matrix with Affymetrix IDs to Gene Symbol matrix (if multiple probes match to a gene, take the max. average value probe as representative for the gene)
		#
		
		# Function to convert probe-based expression matrix to gene-based expression matrix
		# Parameters:
		#   matdat = input matrix with probe rownames,
		#   mat_conv = vector with gene symbols corresponding to probe rownames (NA for missing conversions)
		probe2genemat <- function(matdat, mat_conv)
		{
		
			if(nrow(matdat) != length(mat_conv))
			{
			  stop("Matrix does not have the same number of rows as the gene name vector")
			}
		
			# take max expression vector (max = maximum of mean exp across samples), if gene occurs twice among probes
			unq <- unique(mat_conv)
			if(any(is.na(unq))){
				unq <- unq[-which(is.na(unq))]
			}
			mat <- matrix(0.0, nrow=length(unq), ncol=ncol(matdat))
			for(j in 1:nrow(mat))
			{
			  ind <- which(unq[j]==mat_conv)
		
			  # show conversion progress, every 1000 probes
			  if(j %% 1000 == 0){
			    print(j)
			  }
		
			  # 1-to-1 probe to gene symbol matching
			  if(length(ind) == 1)
			  {
			    mat[j,] = as.numeric(as.matrix(matdat[ind,]))
			  } else if(length(ind) > 1){
		
			    # multiple probes match to one gene symbol
			    curmat = matdat[ind,]
		
			    # compute average expression per row -> select row with max. avg. expression
			    avg = apply(curmat, 1, mean)
			    mat[j,] = as.numeric(as.matrix(matdat[ind[which.max(avg)],]))
			  }
			}
			rownames(mat) = unq
		
		  return(mat)
		}
		
		# Run the conversion from probe matrix to gene matrix (Moran data)
		moran_symb = probe2genemat(moranfilt, mapped_symbols)
		colnames(moran_symb) = colnames(moranfilt)		
		# dim(moran_symb)

	  
	  	# Extract pathway activities
		moran_path = summarized_pathway_activity(moran_symb, gsets=pathlst, type="median", minsize = 10)			
		print(moran_path[1:5,1:5])


		# Random Forest classification analysis		
		
		# install R-packages for classification
		if(!require('randomForest'))
		{
			install.packages('randomForest')
			require('randomForest')
		}
		
		
		# set seed number for reproducibility
		set.seed(1234) 
		
		# Build Random Forest sample classification model for Zhang et al. data using 250 decision trees		
		rfmod_moran= randomForest(t(moran_path), factor(moran_outcomefilt), ntree=250, keep.forest=TRUE)
		
		# show model evluation based on out-of-bag samples
		print(rfmod_moran)
		# OOB estimate of  error rate: 10.26%
		
		# which pathways were most informative for the prediction (multivariate feature selction):
		print(head(rfmod_moran$importance[order(rfmod_moran$importance, decreasing=T),]))		
		
}

