run_gemma<-function(
                    plinkbase,
                    plinkfolder='.',
                    out=".",
                    type='bslmm',
                    maf=0.0,
                    background=TRUE,
                    dryrun=FALSE
                    ){
background=ifelse(background==F, " ", " &")
if(type=='bslmm'){
  command= paste0('nice ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -bslmm ')
}else if(type=='lm'){
  command= paste0('nice  ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -lm 4  ')
}else if(type=='lmm'){
  command0= paste0('nice  ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -gk 1  -o ' ,out)
  command= paste0('nice  ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -k ', paste0(out,'.cXX.txt') , '  -lmm 4 ')
}
if( !type %in% c('bslmm','lm','lmm')) stop('type of GWA not recognized')
# Add maf
command<-paste(command, "-maf", maf)
# Add out
command<-paste(command, ' -o ' ,out)
# Add background
command<-paste(command, background)
# Run gemma
if(dryrun==TRUE){
  message("Dry run, only printing command")
  if(type=='lmm') print(command0)
  print(command)
}else{
message(paste('running GEMMA command: ',command))
  if(type=='lmm') system(command0)
  system(command)
}
return(TRUE)
}
pred_gemma<-function(
  predict=1,
  bfile="../gemma/rFitness_mhi/515g",
  epm="output/rFitness_mhi.param.txt",
  emu="output/rFitness_mhi.log.txt",
  o="rFitness_mhi-predict",
  dryrun=F
){
  command<-paste("nice  ~/bin/gemma -bfile",bfile,
                 "-predict",predict,
                 "-epm",epm,
                 "-emu",emu,
                 "-o",o
  )
  # gemma -bfile ../gemma/rFitness_mhi/515g -epm output/rFitness_mhi.param.txt -emu output/rFitness_mhi.log.txt  -o rFitness_mhi -predict 1
  if(!dryrun) system(command)
  print(command)
  return(TRUE)
}
.read_gemma<-function(folder="output", name, what='heritability'){
    .read_hyperparameters<-function(outfile, hyperparameter=1, MCMCchain=1000,quantiles= c(0.5,0.025,0.975) ){
    hyp<-read.table(outfile,header=TRUE, stringsAsFactors = FALSE) %>% tail(n=MCMCchain)
    # hist(hyp[,hyperparameter], main = outfile)
    quantile(hyp[,hyperparameter],p=quantiles)
  }

  .read_assoc<-function(outfile="output/fecundity_mhi.assoc.txt"){
    d <- data.table::fread(outfile,header=TRUE,stringsAsFactors = FALSE)
    # d <- read.table(outfile,header=TRUE,stringsAsFactors = FALSE)
    d<-dplyr::rename(d,pos=ps)
    d<-mutate(d, SNP= paste(chr,pos,sep="_"),
              env= .cleanname(outfile),
              cumpos=1:nrow(d))

    if (grepl('.param.', outfile)){
      d<-mutate(d,effect= alpha + (gamma*beta) )
      d<-mutate(d,sparseeffect= (gamma*beta) )
      d$BAY<- d$gamma !=0

    }else if(grepl('.assoc.', outfile)){
      d<-mutate(d,effect= beta)
      d$FDR<- d$p_lrt < fdr_level(d$p_lrt)
      d$BONF<- d$p_lrt < 1/nrow(d)
    }

    return(d)
  }
  ##Run
  if(what=='heritability'){
    res=.read_hyperparameters(file.path(folder,paste0(name,".hyp.txt")))
  }else if(what=='lm'){
    res=.read_assoc(file.path(folder,paste0(name,".assoc.txt")))
  }else if(what=="bslmm"){
    res=.read_assoc(file.path(folder,paste0(name,".param.txt")))
  }else if(what=="bv"){
    res=read.table(file.path(folder,paste0(name,".bv.txt")))
  }
    return(res)
}

.cleanname<-function(x){
  sapply(x, function(i){
    strsplit(i, split =  '/', fixed=TRUE) %>% unlist %>% tail(1) %>%
    strsplit(split =  '.', fixed=TRUE) %>% unlist %>%
    head(1)
  })
}

accuracies<-function(y,x){
  lmo<-lm(y~x)
  a<-fn(coefficients(lmo)[1])
  b<-fn(coefficients(lmo)[2])
  lmos<-summary(lmo)
  r2<-fn(lmos$adj.r.squared)
  rho<-fn(cor(y,x,method="spearman"))
  r<-fn(cor(y,x,method="pearson"))
  return(unlist(list("a"=a,"b"=b,"R2"=r2,"rho"=rho,"r"=r)))
}
read_and_top_gwa<-function(parfile,mapfile,nloci=500){
    gammas<-.read_gemma(folder=dirname(parfile),
                        name=gsub(pattern = ".param.txt", "",basename(parfile)),
                        what = "bslmm")
    if(file.exists(mapfile)){
      map<-data.table::fread(mapfile,verbose = F)
    }else{
      map<-data.table::fread(gsub(x=mapfile, pattern = ".map",replacement = ".bim"))
    }

    try(mapping<-match(moiR::fc(gammas[,2]), moiR::fc(map[,2]) ),silent = T)
    if(!exists("mapping")) stop(".map/.bim file do not share SNP names with .param.txt")

    bslmm<-rep(0,nrow(map))
    bslmm[mapping]<-gammas$effect

    m<-which(rank(max(abs(bslmm))-abs(bslmm)) <= nloci )

    # Propose starting point of selection coefficients based on BSLMM
    s<-rep(0,length(m))
    s<-bslmm[m]
    return(list(s=s,m=m,bslmm=bslmm))
}

ptr<-function(p) sapply(p,function(i) 1+(-log((1-i)/(i))))
iptr<-function(h) sapply(h,function(i) (1/(1+exp(1-i))))
ptr05<-function(p) sapply(p,function(i) 1+(-log((0.5-i)/(i))))
iptr05<-function(h) sapply(h,function(i) (0.5/(1+exp(1-i))))
ptr2<-function(p) sapply(p,function(i) 1+(-log((2-i)/(i))))
iptr2<-function(h) sapply(h,function(i) (2/(1+exp(1-i))))
# seltr<-function(s) sapply(s,function(i) 1+ log(1+i))
# iseltr<-function(x) sapply(x,function(i) exp(i-1)-1)
# Mar 29 2020 what if we constrain the selection coefficients to -1 to 1 because it creates problems of extrapolation? 
seltr<-function(s) sapply(s,function(y) -log( (2/(y+1))-1) )
iseltr<-function(x) sapply(x,function(i) (2/(1+exp(-i)))-1 )




nap<-function(bedfile,famfile,mapfile,
                  myrows=NULL,mycols=NULL,
                  s=NULL,mod=1,epi=1,iter=20,k=1,
                  log=NULL,
                  verbose=1){
  # cat("Input files\n")
  # cat(mapfile,"\n")
  # cat(famfile,"\n")
  # cat(bedfile,"\n")
  cat(" Reading .map file\n")
  if(file.exists(mapfile)){
      map<-data.table::fread(mapfile,verbose = F)
    }else{
      map<-data.table::fread(gsub(x=mapfile, pattern = ".map",replacement = ".bim"))
    }
  cat(" Reading .fam file\n")
  fam<-data.table::fread(famfile)

  y<-moiR::fn(fam[,6])
  N<-nrow(fam)
  p=nrow(map)

  if(is.null(mycols)) mycols<-1:p
  if(is.null(myrows)) myrows<-1:N
  nsnps<-length(mycols)

  if(any(is.na(y)) | any(y == -9)){
    cat("All NAs to zero \n")
    y[is.na(y)]<-0
    y[ y == -9 ]<-0
  }

  
  # starting parameters
  # if(is.null(s)) 
  s<-rnorm(nsnps,0,0.01)
  parstart<-par<-list(
                      "s"=seltr(s),
                      "b"=ptr(0.1),
                      "a"=ptr(0.1),
                      "p"=ptr05(0.1),
                      "mu"=ptr2(1)
                      )
  parstart<-unlist(parstart)
  
  # Optimization
  if(!is.null(log)) sink(log)
  start_time <- Sys.time()
    ### Cross-Validation design ###
    if(k==2){
      cat("K=2, using a 90% train 10% test crossvalidation \n" )
      cv<-sample(0:1,size = length(myrows),prob=c(0.1,0.9),replace = T)
      cat("=>Runing on 90% \n" )
      rcv<-napspgC(
            bedfile,
            N, p,myrows[which(cv==1)],mycols,
            y[which(cv==1)],
            parstart,
            epi,
            mod,
            iter,
            verbose=F)
      r_<-list()
      r_$w<-rcv$w
        r_$w[cv==1]<-NA # to just test the 10 %
      r_$y<-y
        r_$y[cv==1]<-NA # to just test the 10 %
    }else if(k>2){
      cat("Dividing data for",k, "fold crossvalidation \n" )
      cv<-sample(1:k,size = length(myrows),replace = T)
      rcv<-lapply(1:k,function(i){
                                  cat("=>Runing group",i," \n" )
                                  napspgC(
                                  bedfile,
                                  N, p,myrows[which(cv!=i)],mycols,
                                  y[which(cv!=i)],
                                  parstart,
                                  epi,
                                  mod,
                                  iter,
                                  verbose=F)
                                    })
      r_<-list()
      r_$w<-lapply(1:k,function(i){rcv[[i]]$w[which(cv==i)]}) %>% unlist
      r_$y<-sapply(1:k, function(i) y[which(cv==i)]) %>% unlist()
    }
    ### Cross-Validation design ###
    r<-napspgC(
              bedfile,
              N, p,myrows,mycols,
              y[myrows],
              parstart,
              epi,
              mod,
              iter,
              verbose)
    r$AIC <- 2*length(parstart) -2 *(-r$f)
    r$y<-y[myrows]
    if(k>1){
      r$y_in<-r$y
      r$w_in<-r$w
      r$y<-r_$y
      r$w<-r_$w
    }else{
      r$y_in<-r$y
      r$w_in<-r$w
    }
    r$myrows<-myrows
    r$mycols<-mycols

  diftime<-Sys.time() - start_time
  cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
  if(!is.null(log)) sink(log)
  ###=======================================================================###
  # Accuracy
  r$w[r$w<0] <-0
  r$rspearman<-cor.test((r$w),(r$y),method='p',na.rm=T)$estimate
  r$rnonzero<-cor.test(r$w[r$y!=0],r$y[r$y!=0],method='p',na.rm=T)$estimate
  r$rspearman_in<-cor.test((r$w_in),(r$y_in),method='p',na.rm=T)$estimate
  r$rnonzero_in<-cor.test(r$w_in[r$y_in!=0],r$y_in[r$y_in!=0],method='p',na.rm=T)$estimate
  
  ###=======================================================================###
  return(r)
}
# napcall<-function(bedfiles,famfiles,mapfiles,parfiles,
#                    es=c(0.9,1,1.1),mods=c(1,2),loci=100,
#                    iter=100,
#                    k=10,
#                    dryrun=T,
#                    override=F,
#                    background=T,
#                    maxbackground=10,
#                    totbackground=50
#                   ){
#   counter<-1
#   for(i in 1:length(bedfiles)){
#     for(e in es){
#       for(mod in mods){
#         for(l in loci){
#         # background<-ifelse(system("ps | grep R | wc -l")<totbackground,
#         #                    " &"," ")
#         # if((counter %% maxbackground)==0) background<-" "
#         command<-paste("srun -c 8 Rscript napcall.R",
#                        "--param", parfiles[i],
#                        "--fam", famfiles[i],
#                        "--map", mapfiles[i],
#                        "--bed", bedfiles[i],
#                        "--e", e,
#                        "--m", mod,
#                        "--l", l,
#                        "--o", override,
#                        "--i", iter ,
#                        "--k", k ,
#                        " ", background
#                        )
#         message(command)
#         if(!dryrun){system(command)}
#         counter=counter+1
#         }
#       }
#     }
#   }
# }
napcall<-function(bedfiles,famfiles,mapfiles,parfiles,
                   es=c(0.9,1,1.1),mods=c(1,2),loci=100,
                   iter=100,
                   k=5,
                   dryrun=T,
                   override=F
                  ){
  counter<-1
  for(i in 1:length(bedfiles)){
    for(e in es){
      for(mod in mods){
        for(l in loci){
        command<-paste("sbatch napcall.sh ",
                       "--param", parfiles[i],
                       "--fam", famfiles[i],
                       "--map", mapfiles[i],
                       "--bed", bedfiles[i],
                       "--e", e,
                       "--m", mod,
                       "--l", l,
                       "--o", override,
                       "--i", iter ,
                       "--k", k
                       )
        
        
        dname<-dirname(parfiles[i])
        bname<-gsub(x=basename(parfiles[i]),pattern=".param.txt",replacement="")
        finalbase=paste0(dname,"/",bname,".results.",
                    ifelse(mod==1,"a_","m_"),
                    paste0("e",e,"_"),
                    paste0("l",l)
                  )
        # finalfile<-paste0(finalbase,".tsv")
        # finallog<-paste0(finalbase,".log")
        # finalrda<-paste0(finalbase,".rda")
        # finaltsv<-paste0(finalbase,".tsv")
        successfile<-paste0(finalbase,".success")
        if(all(file.exists(c(successfile)) ) ){
          prevAIC<-fn(head(read.table(successfile),1))
          if(abs(prevAIC) < 1e+200){
            if(override){
              message("Overwriting mode!")
              if(!dryrun){system(command)}  
            }else{
              message("Result files already exist and run AIC does not indicate failed optimization! Stopping now")
              next
            }
          }else{
              message("Result files already exist but AIC indicates a failed optimization run")
              if(!dryrun){system(command)}      
          }
        }else{
          message(command)
          if(!dryrun){system(command)}  
        }
        
        }
      }
    }
  }
}




