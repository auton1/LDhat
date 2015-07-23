###############################################################
# Various functions to run pop genetic analyses and practicals
###############################################################

#library(popgen);

#########################################################
# To read in genotype data (0,1,2,? as data types
# could also implement PHASE type. 1 is heterozygote)
#########################################################

read.genotype.data<-function(file, data.format="unphased", names=FALSE, header=TRUE, positions=TRUE) {

	if (header == TRUE) {data <- read.table(file, fill=TRUE, na.strings="?", skip=1, header=FALSE);}
	else {data <- read.table(file, fill=TRUE, na.strings="?", header=FALSE);}

	if (positions == TRUE) {n.seq <- nrow(data)-1;}
	else {n.seq <- nrow(data);}
	if (names == TRUE) {l.seq <- ncol(data)-1;}
	else {l.seq <- ncol(data);}

	if (positions == TRUE) {
		pos <- as.numeric(as.vector(data[1,]));
		pos <- pos[1:l.seq];
		data <- data[2:nrow(data),];
	}
	else {pos <- c(1:l.seq);}

	if (names == TRUE) {
		nm<- as.vector(data[,1]);
		data <- data[,2:ncol(data)];
	}
	else {nm <- c(1:n.seq);}
	
	seq <- as.matrix(data);
;
	return(list(type="genotype", n.seq=n.seq, l.seq=l.seq, pos=pos, names=nm, seq=seq));
}


#Reads in haplotype data in PHASE-like format

read.haplotype.data<-function(file, data.format="phased", names=FALSE, header=TRUE, positions=TRUE) {

	if (header == TRUE) {data <- read.table(file, fill=TRUE, na.strings="?", skip=1, header=FALSE);}
	else {data <- read.table(file, fill=TRUE, na.strings="?", header=FALSE);}

	if (positions == TRUE) {n.seq <- nrow(data)-1;}
	else {n.seq <- nrow(data);}
	if (names == TRUE) {l.seq <- ncol(data)-1;}
	else {l.seq <- ncol(data);}

	if (positions == TRUE) {
		pos <- as.numeric(as.vector(data[1,]));
		pos <- pos[1:l.seq];
		data <- data[2:nrow(data),];
	}
	else {pos <- c(1:l.seq);}

	if (names == TRUE) {
		nm<- as.vector(data[,1]);
		data <- data[,2:ncol(data)];
	}
	else {nm <- c(1:n.seq);}
	
	seq <- as.matrix(data);

	return(list(type="haplotype", n.seq=n.seq, l.seq=l.seq, pos=pos, names=nm, seq=seq));


}

#Reads the output from Selsim

 read.selsim.data<-function(seq.file="seq.txt", tree.file="SStree.txt", plot=TRUE) {

  x=read.table(seq.file,skip=4,fill=T, as.is=T)
  pos<-as.matrix(x[1,]);
  x<-as.matrix(x[2:(nrow(x)-2),]);
  data <- list(type="haplotype", n.seq=nrow(x), l.seq=ncol(x), pos=pos, names=c(1:nrow(x)), seq=x);
  data <- subset.data(data, maf=0)
  if (plot==TRUE) {
  tree<-read.table(tree.file, as.is=T, header=T);
  x11();
  o<-plot.tree(tree)
  haplotype.plot(data, order=o);
  }
  return(data);
}


#Generic wrapper for reading data files

read.data <- function(file, format="haplotype", names=FALSE, positions=TRUE, header=TRUE, aligned=TRUE, gap.as.state=TRUE, biallelic.only=TRUE) {

	if (format=="selsim")           {data<-(read.selsim.data(file));}
	else if (format == "haplotype") {data<-(read.haplotype.data(file, names=names, positions=positions, header=header));}
	else if (format == "genotype")  {data<-(read.genotype.data(file, names=names, positions=positions, header=header));}
	else if (format == "fasta") {data <- read.fasta(file, aligned=aligned, gap.as.state=gap.as.state, biallelic.only=biallelic.only);}
	
	cat(paste("\n\nRead ",data$n.seq, " sequences of length ",data$l.seq, "\n\n",sep=""));

	return(data);
}



#snp.frequency

snp.frequency <- function(data, maf=TRUE) {

	ct<-apply(data$seq, 2, sum, na.rm=T)
	norm<-apply(!is.na(data$seq), 2, sum)
	if (data$type == "genotype") norm <- norm*2;

	freq<-ct/norm;

	if (maf==TRUE) {
		f.inv<-1-freq;
		freq<-apply(cbind(freq, f.inv), 1, min);
	}
	return(freq)
}


#Selects data subset on basis of position and MAF


subset.data <- function(data, maf=-1, pos.lim=c(), select.snp=c(), select.ind=c(), convert.to.minor=FALSE, missing.data.thresh=0.2) {

	#Just keep those individuals desired
	if (length(select.ind)>0) {
		data$seq<-data$seq[select.ind,];
		data$n.seq<-length(select.ind);
		data$names<-data$names[select.ind];
	}


	#Remove SNPs with >thresh missing data
	na<-apply(is.na(data$seq), 2, mean);
	data$seq<-data$seq[,na<missing.data.thresh];
	data$pos<-data$pos[na<missing.data.thresh];
	if (length(data$snp.id>0)) data$snp.id<-data$snp.id[na<missing.data.thresh];
	data$l.seq<-sum(na<missing.data.thresh);

	#Remove inds with >thresh missing data
	na<-apply(is.na(data$seq), 1, mean);
	data$seq<-data$seq[na<missing.data.thresh,];
	data$names<-data$names[na<missing.data.thresh];
	data$n.seq<-sum(na<missing.data.thresh);

	select <- c(1:data$l.seq)
	
	#On allele frequency
	if (maf>=0) {
		f <- snp.frequency(data);
		select <- select[f>maf];
	}

	#On limits	
	if (length(pos.lim) == 2) {
		select <- select[data$pos[select]>pos.lim[1] & data$pos[select]<pos.lim[2]];
	}

	#On defined snps - identified by relative position - problematic if using filters from above
	if (length(select.snp)>0) {
		select<-select[select %in% select.snp];
	}

	if (length(select) == 0) {
		cat("\nNo data to return\n");
		return(NULL);
	}

	if (convert.to.minor==TRUE) {
		f<-apply(data$seq[,select], 2, mean);
		data$seq[,select[f>0.5]]<-1-data$seq[,select[f>0.5]]
	}

	cat(paste("\n\nNumber of SNPs returned = ", length(select), "\n\n", sep=""));
	
	if (length(data$snp.id>0)) {
		return(list(type=data$type, n.seq=data$n.seq, l.seq=length(select), pos=data$pos[select], snp.id=data$snp.id[select],
			names=data$names, seq=as.matrix(data$seq[,select])));
	}
	else {
		return(list(type=data$type, n.seq=data$n.seq, l.seq=length(select), pos=data$pos[select], 
			names=data$names, seq=as.matrix(data$seq[,select])));
	}


}


#Summarises data

summarise.data <- function(data, single.plot=TRUE, ld.plot=TRUE, maf=0.0, n.int=10, new.plot=TRUE) {

	if (data$l.seq<2) {
		cat("\n\n*** Error: too little data to analyse (S<2) ***\n\n");
	}

	data<-subset.data(data, maf=0.0)

	if (data$l.seq<2) {
		cat("\n\n*** Error: too little data to analyse (S<2) ***\n\n");
	}

	cat("\n\n#####################################################\n");
	cat(" Coalescent summary of data\n");
	cat("#####################################################\n\n");
	
	if (data$type=="haplotype") cat(paste("Number of sequences = ",data$n.seq,"\n",sep=""));
	if (data$type=="genotype") cat(paste("Number of individuals = ",data$n.seq,"\n",sep=""));
	cat(paste("Number of polymorphic sites = ",data$l.seq,"\n",sep=""));
       

        if (new.plot) x11();
        if (single.plot == TRUE) par(mfrow=c(2,2));

	freq <- snp.frequency(data);

	#Count singletons
	fq<-apply(data$seq, 2, sum, na.rm=T);
	n.na<-apply(is.na(data$seq), 2, sum);
	if (data$type=="haplotype") n.sing <- sum(fq==1)+sum(fq==(data$n.seq-n.na-1));
	if (data$type=="genotype") n.sing <- sum(fq==1)+sum(fq==(2*data$n.seq-2*n.na-1));
	if (n.sing<1) n.sing<-0;

	d <- dist(data$seq, method="manhattan");
	if (data$type=="genotype") d<-d/2;

	pwd = mean(d);
        hl=hclust(d, method="average");
	ns <- data$n.seq;
	if (data$type=="genotype") ns<-2*ns;
	s <- data$l.seq;

        con=vector(length=10);
        con[1]=sum(1/c(1:(ns-1)));
        con[2]=sum((1/c(1:(ns-1)))^2);
        con[3]=(ns+1)/(3*ns-3);
        con[4]=2*(ns*ns+ns+3)/(9*ns*ns-9*ns);
        con[5]=con[3]-1/con[1];
        con[6]=con[4]-(ns+2)/(ns*con[1])+con[2]/(con[1]*con[1]);
        con[7]=con[5]/con[1];
        con[8]=con[6]/(con[1]*con[1]+con[2]);
        con[9]=2*(ns*con[1]-2*(ns-1))/((ns-1)*(ns-2));
        con[10]=con[9]+(ns-2)/((ns-1)*(ns-1))+((2/(ns-1))*(1.5-(2*(con[1]+1/ns)-3)/(ns-2)-1/ns));
        con[11]=(ns*ns*con[2]/((ns-1)*(ns-1))+con[1]*con[1]*con[10]-2*ns*con[1]*(con[1]+1)/((ns-1)*(ns-1)))/(con[1]*con[1]+con[2]);
        con[12]=ns/(ns-1)*(con[1]-ns/(ns-1))-con[11];

        cat(paste("\nWatterson's estimate of theta = ", round(data$l.seq/con[1],2), "\n", sep=""));
        cat(paste("Average paiwise differences = ", round(pwd,2), "\n", sep=""));
        cat(paste("Number of singletons = ", n.sing, "\n", sep=""));
        cat(paste("Tajima D statistic = ", round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2), "\n", sep=""));
        cat(paste("Fu and Li D statistic = ", round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2), "\n", sep=""));

	breaks = c(-1.5,-0.5,0.5,1.5,2.5)
        image(x=data$pos, z=t(data$seq[hl$order,]), col=c("white", "blue", "yellow","red"), breaks=breaks, xlab="Position", ylab="Sample - same order as cluster", main="Haplotype plot", yaxt="n");

	if (single.plot == FALSE) {x11();}
	breaks<-c(0:n.int)/(2*n.int);
	breaks[length(breaks)]<-breaks[length(breaks)]*1.1;
	exp<-rep(0, n.int);
	for (i in 1:(data$n.seq-1)) {
		wt<-1/i;
		f.i<-findInterval(min(i, data$n.seq-i)/data$n.seq, breaks);
		exp[f.i]<-exp[f.i]+wt;
	}
	obs<-hist(freq, breaks=breaks, plot=F)$counts;
	exp<-exp*length(freq)/sum(exp);
	mpt<-(breaks[1:n.int]+breaks[2:(n.int+1)])/2;
	barplot(rbind(obs, exp), xlab="Minor Allele Frequency", ylab="Count", beside=T,	
		col=c("blue", "cyan"), main="Minor Allele Frequency", names=mpt, space=c(0,0.25),
		cex.names=0.5);
	


	if (single.plot == FALSE) {x11();}
        hist(d, main="Pairwise differences", xlab="Pairwise differences", ylab="Count", col="blue");

	if (single.plot == FALSE) x11();
	plot(x=hl, labels=data$names, main="UPGMA tree of pairwise differences", ylab="Distance", hang=-1);

	if (ld.plot == TRUE) {
		x11();
		par(mfrow=c(1,1));
          	plotld(data, cutoff=0.05);
	}
	
}


#Makes LD matrix from data

makeld<-function(data){

  cor=array(0, dim=c(data$l.seq, data$l.seq))
  

    for(i in 1:(nrow(cor)-1)){

      for(j in (i+1):ncol(cor)){

	which.use<-which(!is.na(data$seq[,i]) & !is.na(data$seq[,j]));

	if (length(which.use>0)) {

        	if (data$type=="haplotype") {
			d<-(length(which.use)-1)/length(which.use)*cov(data$seq[which.use,i],data$seq[which.use,j], use="pairwise");
			fi<-mean(data$seq[which.use,i]);
			fj<-mean(data$seq[which.use,j]);
		}
		if (data$type=="genotype") {
			d = estimateD.genotype.EM(data$seq[which.use,c(i,j)]);
			fi<-d$fi;
			fj<-d$fj;
			d<-d$D;
		}
		r<-d/sqrt(fi*fj*(1-fi)*(1-fj));
		cor[i,j]=r^2;
        	if(!is.na(d) & d>=0) cor[j,i]=d/min(fi*(1-fj),(1-fi)*fj);
        	if(!is.na(d) & d<0) cor[j,i]=-d/min(fi*fj,(1-fi)*(1-fj));
	}
	if (length(which.use)==0) {
		is.na(cor[j,i])<-TRUE;
		is.na(cor[i,j])<-TRUE;
	}
    }
  }

  if (data$type=="haplotype") {
    min.rec = vector(length=ncol(cor));
    min.rec[1]=0;
    if (!is.na(cor[2,1]) & cor[2,1]<0.999) {min.rec[2]=1;}
    else {min.rec[2]=0;}

    for (k in 3:ncol(cor)) {
	min.rec[k]=0;
	for (i in 2:(k-1)) {
		if (!is.na(cor[k,i]) & cor[k,i]<0.999) {r.add=1;}
		else {r.add=0;}
		if (min.rec[i]+r.add>min.rec[k]) {min.rec[k]=min.rec[i]+r.add;}
	}
    }
  }

    if (data$type=="genotype") min.rec<-min.rec.gt(data);
    cat(paste("\n\nMinimum number of recombination events (HK85) = ",min.rec[ncol(cor)],"\n\n", sep=""));
  


  return(cor);


}

plotld=function(data,cutoff=0,left=-1,right=-1,cor=c(), col=c(), breaks=c()){


  data<-subset.data(data, maf=cutoff);

  cat("\n\nCalculating LD matrix");
  
  cor=makeld(data)

  if (length(col)==0) {
  	blues=rgb(0,0,(20:0)/20);
  	reds=rgb((0:20)/20,0,0);
	col=c(blues, reds);
  }

  if (length(breaks)==0) {
	image(data$pos,data$pos,cor,zlim=c(0,1.02),xlab="Position",ylab="Position", 
		main="LD matrix: r2 - |D'|", col=col)
  }
  else {
	image(data$pos,data$pos,cor,zlim=c(0,1.5),xlab="Position",ylab="Position", 
		main="LD matrix: r2 - |D'|", col=col, breaks=breaks)
  }
}




#########################################################
# Haplotype plot of data
#########################################################

haplotype.plot<-function(data, maf=-1, plot.lims=c(), plot.new=TRUE, col=c(), order=c(), spacing="distance", 
	type="block", return=FALSE, add.legend=TRUE, cex=1, sort.hap.position=0, f.half=0.2, method="normal") {

	if (data$l.seq<2) {
		cat("\n\n*** Error: too little data to analyse (S<2) ***\n\n");
	}

	if (maf>=0 | length(plot.lims)>0) data<-subset.data(data, maf=maf, pos=plot.lims);


        if (plot.new==TRUE) {x11();}
	if (length(order)==0) {
 		if (sort.hap.position==0) {d <- dist(data$seq, method="manhattan");}
		else {d<-make.weighted.dist(data, position=sort.hap.position, f.half<-f.half, method=method);}
        	hl<-hclust(d, method="average");
		order<-hl$order
	}

	if (is.na(pmatch(spacing, "distance"))) {pos<-c(1:data$l.seq);}
	else {pos<-data$pos;}

	#Block like picture
	if (is.na(pmatch(type, "block"))==FALSE) {
		breaks = c(-1.5,-0.5,0.5,1.5,2.5)
		if (length(col)==0) col=c("white", "blue", "yellow","red");
        	image(x=pos, z=t(data$seq[order,]), col=col, breaks=breaks, 
			xlab="Position", ylab="Chromosome", main="Haplotype plot", yaxt="n");
	}

	#Stick and ball picture
	else {
		plot(0,0,type="n", xlim=c(min(pos), max(pos)), ylim=c(0,(data$n.seq*1.05)), yaxt="n", xlab="Position", bty="n", ylab="");
		if (length(col)==0) {col=rep("blue", data$l.seq);}
		for (i in 1:data$n.seq) segments(x0=pos[1], y0=i, x1=pos[length(pos)], y1=i, lwd=1, col=grey(0.85));
		for (i in 1:data$n.seq) {
			snps<-pos[data$seq[order[i],]==1]
			points(x=snps, y=rep(i, length(snps)), pch=19, col=col[data$seq[order[i],]==1], cex=cex);
		}
		if (is.na(pmatch(spacing, "distance")) & add.legend==TRUE) {
			segments(x0=pos[1], y0=data$n.seq*1.05, x1=pos[length(pos)], y1=data$n.seq*1.05, lwd=1, col=grey(0.75));
			low<-data$pos[1];
			high<-data$pos[data$l.seq];
			snp.pos<-1+(data$l.seq-1)*(data$pos-low)/(high-low);
			segments(snp.pos, data$n.seq*1.04, snp.pos, data$n.seq*1.05, lwd=1, col=col);
			segments(pos, data$n.seq+1, snp.pos, data$n.seq*1.03, lwd=1, col="black");
		}
	}

	if (return==TRUE) return(order);

}



#############################################################
# Make distance matrix weighting exponentially from location
#############################################################

make.weighted.dist<-function(data, position=0, f.half=1, method="exponential") {

	if (method=="exponential") {
		lambda<-log(2)/(f.half*(data$pos[length(data$pos)]-data$pos[1]));
		wts<-exp(-lambda*abs(data$pos-position));
	}
	else {
		sig<-f.half*(data$pos[length(data$pos)]-data$pos[1])/sqrt(2*log(2));
		wts<-exp(-(data$pos-position)^2/(2*sig^2));
	}

	d<-matrix(nrow=data$n.seq, ncol=data$n.seq);
	for (i in 1:(data$n.seq-1)) for (j in (i+1):data$n.seq) {
		d[j,i]<-sum(((data$seq[i,]-data$seq[j,])^2)*wts);
	}
	return(as.dist(d));

}


#########################################################
#Summarises output from interval program in LDhat2.1
#########################################################


summarise.interval=function(rates.file="rates.txt", burn.in=30, locs.file=FALSE) {

  x = read.table(rates.file, skip=1, fill=T);
  x = as.matrix(x);

  low = as.integer(nrow(x)*burn.in/100);

  cat("\n\nSummarise output from MCMC estimation of recombination rates in INTERVAL (LDhat 2.1)\n\n");
  cat(paste("Number of SNPs = ", ncol(x), "\n", sep=""));
  cat(paste("Number of samples = ", nrow(x), "\n", sep=""));
  cat(paste("Burn-in period = ", low, " samples\n", sep=""));

  
  x11();
  par(mfrow=c(1,2));
  plot(x[,1], type="s", col=rgb(0,0,0.5), xlab="Sample", ylab="Total map length", 
	main="Mixing of total map length");
  image(x=c(1:nrow(x)), y=c(1:(ncol(x)-1)), z=log(x[,2:ncol(x)]), xlab="Sample", 
	ylab="log(rate) at SNP", main="Mixing of rates");

  means<-apply(x[low:nrow(x),], 2, mean, na.rm=T);
  q.95<-apply(x[low:nrow(x),], 2, quantile, probs=c(0.025, 0.5, 0.975), na.rm=T);

  cat(paste("\nMean posterior total map length (4Ner) = ", signif(means[1], 4), "\n", sep=""));
  
  x11();
  if (locs.file==FALSE) {pos<-c(1:ncol(x)); xlab<-"Position (SNP)";}
  else {pos<-as.vector(as.matrix(read.table(locs.file, as.is=T, skip=1))); xlab<-"Position"}

  plot(pos[1:(length(pos)-1)], y=means[2:length(means)], type="s", col=rgb(0,0,0.5), 	
	xlab=xlab, ylab="Posterior mean rate", main="Posterior mean rates");
  lines(pos[1:(length(pos)-1)], y=q.95[1,2:length(means)], type="s", col=grey(0.75), lty="dotted");
  lines(pos[1:(length(pos)-1)], y=q.95[3,2:length(means)], type="s", col=grey(0.75), lty="dotted");

  op<-cbind(means, t(q.95));
  colnames(op)<-c("Mean", "q2.5", "Median", "q97.5");

  return(op);


}


#########################################################
#summarises output from rhomap program in LDhat2.1
#########################################################

summarise.rhomap<-function(rates.file="rates.txt", burn.in=30, locs.file=FALSE) {

	return(summarise.interval(rates.file=rates.file, burn.in=burn.in, locs.file=locs.file));

}



#########################################################
#summarises output from pairwise program in LDhat2.1
#########################################################


summarise.pairwise<-function(surf=TRUE, window=FALSE, rm=FALSE, test=FALSE, ci=FALSE, locs.file=FALSE) {

  cat("\n\nSummarising output from PAIRWISE program in LDhat 2.1\n\n");

  if (locs.file==TRUE) pos<-as.vector(as.matrix(read.table(locs.file, as.is=T, skip=1)));
  
  if (surf) {
    surface = read.table("outfile.txt", skip=10, fill=T);
    x11();
    plot(surface[,1:2], type="l", col=rgb(0,0,0.5), xlab="4Ner", 
		ylab="Composite likelihood", main="Likelihood surface for 4Ner over region");

    xmax=surface[1,1];
    ymax = surface[1,2];

    for (i in 2:nrow(surface)) {
      if (surface[i,2]>ymax) {ymax=surface[i,2]; xmax=surface[i,1];}
    }
    cat(paste("Maximum composite likelihood estimate of 4Ner = ", xmax, "\n", sep=""));
  }

  if (window) {
    win = read.table("window_out.txt", skip=5);
    x11();
    par(mfrow=c(3,1));
    plot(x=(win[,1]+win[,2])/2, win[,3], type="n", xlab="Window", ylab="Number SNPs", main="SNP density");
    segments(x0=win[,1], y0=win[,3], x1=win[,2], y1=win[,3], col="blue");
    plot(x=(win[,1]+win[,2])/2, win[,4], type="n", xlab="Window", ylab="4Ner per kb/bp", main="Local recombination rate");
    segments(x0=win[,1], y0=win[,4], x1=win[,2], y1=win[,4], col="red");

    av = vector(length=nrow(win));
    for(i in 1:nrow(win)) av[i]=xmax/(win[nrow(win),2]-win[1,1]);
    lines(x=(win[,1]+win[,2])/2, y=av, col="green");
    plot(x=(win[,1]+win[,2])/2, y=win[,5], type="n", xlab="Window midpoint", ylab="CLR", main="Composite likelihood ratio");
    segments(x0=win[,1], y0=win[,5], x1=win[,2], y1=win[,5], col="black");
  }

  if (rm) {
    rmin = read.table("rmin.txt", skip=5, fill=T);
    rmin= as.matrix(rmin[1:nrow(rmin),2:ncol(rmin)]);

    f<-as.matrix(read.table("freqs.txt", as.is=T, skip=5)[,2:6]);
    n.chr<-sum(f[1,]);
    anc<-apply(f[,2:5], 1, sum);
    mx<-apply(f[,2:5], 1, max);
    seg<-mx<anc;
    
    rmin<-rmin[seg[1:(length(seg)-1)],];
    rmin<-rmin[,seg[2:length(seg)]];

    maxrm=0;
    maxscrm=0;
    for (i in 1:nrow(rmin)) {
      if (i<ncol(rmin)) {for (j in (i):ncol(rmin)) {if (rmin[i,j]>maxrm) maxrm=rmin[i,j];}}
      if (i>1) {for (j in 1:(i-1)) {if(rmin[i,j]>maxscrm) maxscrm=rmin[i,j];}}
    }

    xu=matrix(0,nrow=nrow(rmin)+1, ncol=nrow(rmin)+1);
    xl=matrix(0,nrow=nrow(rmin)+1, ncol=nrow(rmin)+1);
    
    for (i in 1:nrow(rmin)) {
      if (i<ncol(rmin)) {
        for (j in (i):ncol(rmin)) {xu[i,j+1]=rmin[i,j]; xl[i,j+1]=0;}
      }
      if (i>1) {
        for (j in 1:(i-1)) {xl[i,j]=rmin[i,j]; xu[i,j+1]=0;}
      }
      xu[i,i]=0;
      xl[i,i]=0;
    }

    blues = c(rgb(0,0,(15:6)/15));
    reds = c(rgb((6:15)/15,0,0));

    x11();
    par(mfrow=c(2,1));
    
    image(xu,  xaxt="n", yaxt="n", xlab="Position", ylab="Position", main="Total recombination", 
	col=c("white", blues, reds, "yellow"));
    image(log(xl+0.001), main="Recombination scaled by distance", xaxt="n", yaxt="n", 
	xlab="Position", ylab="Position", col=c("white", blues, reds, "yellow"), 
	breaks=as.vector(quantile(xl, probs=seq(0,1,length.out=23))));
    


    rm(xu);
    rm(xl);
    rm(rmin);
  }

  if (test) {
    xx<-max(surface[,2]);
    rd<-as.matrix(read.table("rdist.txt", skip=4)[,2:5]);

    x11();
    xlim=range(rd[,2]);
    hist(rd[,2], col="blue", breaks=seq(xlim[1], xlim[2], length.out=20), , 
	main="Likelihood permutation test", xlab="Composite Likelihood");
    points(xx,1,pch=25, col="red", bg="red");
  }

  if (ci) {
	cat("\n\nNothing implemented yet to display sampling distribution\n\n");
  }

}



#########################################################
#Function to plot the results of Jonathan's ps program
#########################################################

plot.ps.results=function(ps.output) {

  x11();
  pal = c(rgb(0,0,(15:0)/15), rgb((0:15)/15,0,0));
  image(x=c(1:nrow(ps.output$Q)), y=c(1:ncol(ps.output$Q)), z=ps.output$Q, xlab="Individual", ylab="Posterior probability from K",col=pal, main="Posterior probs for individual assignment");

  x11();
  pal = c("red", "blue", "green", "cyan", "black", "purple", "orange");
  cmax = max(ps.output$c);
  plot(ps.output$c[,1], type="n", xlab="Sample", ylab="c-parameter", ylim=c(0,cmax), main="Samples of the population c parameters");
  for (i in 1:ncol(ps.output$c)) {
	lines(ps.output$c[,i], col=pal[i]);
  }
}


#########################################################
#Writes LDhat format file from data - need to cope with missing data still
#currently outputs as NA rather than ?/N/n/-
#########################################################

write.ldhat.format = function(x, seq.file="seqs.txt", loc.file="locs.txt", in.kb=TRUE, maxseq=-1, na.char="?") {
	
	if (x$type == "haplotype") {hd=1;}
	else {hd=2;}

	if ("names" %in% names(x)) {
		names<-x$names;
	}
	else {
		names<-paste("seq", 1:x$n.seq, sep="");
	}

	if (in.kb) {
		write(paste(x$l.seq, as.character(max(x$pos)/1000), "L", sep=" "), loc.file);
	}
	else {
		write(paste(x$l.seq, max(x$pos), "L", sep=" "), loc.file);
	}
	if (in.kb) {write(as.character(x$pos/1000), file=loc.file, append=TRUE, ncolumns=1);}
	else {
		write(x$pos, loc.file, append=TRUE, ncolumns=1);
	}

	if (maxseq>0) {
		n.seq <- min(maxseq, x$n.seq);
	}
	else {n.seq <- x$n.seq;}

	cat(paste(n.seq, x$l.seq, hd, "\n", sep=" "), file=seq.file);
	for (i in 1:n.seq) {
		
		s<-as.character(x$seq[i,]);
		s[is.na(s)]<-na.char;
		#s[s=="2"]<--1;
		#s[s=="1"]<-2;
		#s[s=="-1"]<-1;
		
		#cat (paste(">", names[i], "\n", s, "\n", sep="", coll=""), file=seq.file, append=TRUE);
		cat (paste(">", names[i], "\n", sep=""), file=seq.file, append=TRUE);
		cat (paste(s, sep="", collapse=""), file=seq.file, append=TRUE);
		cat ("\n", file=seq.file, append=TRUE);
	}

}

#########################################################
#Writes PHASE format file from data - need to cope with missing data still
#currently outputs as NA rather than ?/N/n/-
#########################################################

write.phase.format = function(x, file="data.phase") {
	
	cat(x$pos, sep="\t", file=file);
	for (i in 1:x$n.seq) {
		cat("\n", file=file, append=T);
		cat(x$seq[i,], file=file, append=T);
	}
}



#########################################################
#Convert from standard data format to ps format
#########################################################


convert2ps <- function(data) {

	y = array(dim=c(data$n.seq, 2, data$l.seq));

	if (data$type == "haplotype") {
		for (i in 1:data$n.seq) {
			for (j in 1:data$l.seq) {
				y[i,1,j]=data$seq[i,j];
				y[i,2,j]=-1;
			}
		}
	}

	else if (data$type == "genotype") {
		for (i in 1:data$n.seq) {
			for (j in 1:data$l.seq) {
				if (data$seq[i,j]==0) {y[i,1,j]=0; y[i,2,j]=0;}
				else if (data$seq[i,j]==1) {y[i,1,j]=0; y[i,2,j]=1;}
				else if (data$seq[i,j]==2) {y[i,1,j]=1; y[i,2,j]=1;}
				else {y[i,1,j]=-1; y[i,2,j]=-1;}
			}
		}
	}
  return(y);


}



#########################################################
#To plot haplotype structure around a SNP
#########################################################

plot.snp.hap<-function(data, snp=1, pos=c(), maf=0.0, gap=3, cex=2, scale=10, min.r2=0.2, show=0) {

	snppos<-data$pos[snp];
	d<-subset.data(data, pos=pos, maf=maf);
	snp<-c(1:length(d$pos))[d$pos==snppos];
	d.l<-d$seq[,1:(snp-1)];
	d.r<-d$seq[,(snp+1):d$l.seq];

	f1pos<-sum(d$seq[,snp]);
	
	haps.l<-unique(d.l);
	haps.r<-unique(d.r);
	nhaps.l<-nrow(haps.l);
	nhaps.r<-nrow(haps.r);

	dist.l<-dist(haps.l, method="manhattan");
	dist.r<-dist(haps.r, method="manhattan");
	h.l<-hclust(dist.l, method="average");
	h.r<-hclust(dist.r, method="average");

	hap.count.l<-matrix(nrow=nhaps.l, ncol=6);
	hap.id.l<-matrix(nrow=nhaps.l, ncol=1);
	hap.count.r<-matrix(nrow=nhaps.r, ncol=6);
	hap.id.r<-matrix(nrow=nhaps.r, ncol=1);
	colnames(hap.count.l)<-c("Haplotype", "Count0", "Count1", "Freq", "D", "r2");
	colnames(hap.count.r)<-c("Haplotype", "Count0", "Count1", "Freq", "D", "r2");

	for (hap in 1:nhaps.l) {
		hap.id.l[hap,1]<-paste(haps.l[hap,], sep="", collapse="");
		hap.count.l[hap,1]<-hap;
		m<-t(d.l)-haps.l[hap,];
		mn<-apply(m, 2, min);
		mx<-apply(m, 2, max);
		id<-c(1:ncol(m))[mn==0 & mx==0];
		allele<-d$seq[id,snp];
		f1hap<-sum(allele);
		fhap<-length(id);
		hap.count.l[hap,3]<-f1hap;
		hap.count.l[hap,2]<-fhap-f1hap;
		hap.count.l[hap,4]<-fhap/d$n.seq;
		dnum<-f1hap*d$n.seq-fhap*f1pos;
		dden<-f1pos*(d$n.seq-f1pos)*fhap*(d$n.seq-fhap);
		hap.count.l[hap,5]<-dnum;
		hap.count.l[hap,6]<-dnum*dnum/dden;
	}

	for (hap in 1:nhaps.r) {
		hap.id.r[hap,1]<-paste(haps.r[hap,], sep="", collapse="");
		hap.count.r[hap,1]<-hap;
		m<-t(d.r)-haps.r[hap,];
		mn<-apply(m, 2, min);
		mx<-apply(m, 2, max);
		id<-c(1:ncol(m))[mn==0 & mx==0];
		allele<-d$seq[id,snp];
		f1hap<-sum(allele);
		fhap<-length(id);
		hap.count.r[hap,3]<-f1hap;
		hap.count.r[hap,2]<-fhap-f1hap;
		hap.count.r[hap,4]<-fhap/d$n.seq;
		dnum<-f1hap*d$n.seq-fhap*f1pos;
		dden<-f1pos*(d$n.seq-f1pos)*fhap*(d$n.seq-fhap);
		hap.count.r[hap,5]<-dnum;
		hap.count.r[hap,6]<-dnum*dnum/dden;
	}

	plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", 
		xlim=c(0,d$l.seq+1+2*gap), ylim=c(0,max(nhaps.l, nhaps.r)), bty="n");
	pal<-c("white", "black");
	mpt<-max(nhaps.l, nhaps.r)/2;

	for (hap in 1:nhaps.l) if (hap.count.l[hap,6]>min.r2) {
		col<-"red";
		if (show==0 & hap.count.l[h.l$order[hap],5]>0) {col<-"blue"};
		if (show==1 & hap.count.l[h.l$order[hap],5]<0) {col<-"blue"};
		segments(x0=snp+gap, y0=mpt, x1=snp-1, y1=hap, lwd=scale*hap.count.l[h.l$order[hap],6], col=col);
	}
	for (hap in 1:nhaps.r) if (hap.count.r[hap,6]>min.r2) {
		col<-"red";
		if (show==0 & hap.count.r[h.r$order[hap],5]>0) {col<-"blue"};
		if (show==1 & hap.count.r[h.r$order[hap],5]<0) {col<-"blue"};
		segments(x0=snp+gap, y0=mpt, x1=snp+2*gap, y1=hap, lwd=scale*hap.count.r[h.r$order[hap],6], col=col);
	}

	for (hap in 1:nhaps.l) {
		lines(x=c(1, ncol(haps.l)), y=c(hap, hap));
		points(x=c(1:ncol(haps.l)), y=rep(hap, ncol(haps.l)), pch=21, bg=pal[haps.l[h.l$order[hap],]+1], cex=cex);
	}
	for (hap in 1:nhaps.r) {
		lines(x=c(snp+2*gap+1, ncol(haps.r)+snp+2*gap), y=c(hap, hap));
		points(x=c(1:ncol(haps.r))+snp+2*gap, y=rep(hap, ncol(haps.r)), pch=21, bg=pal[haps.r[h.r$order[hap],]+1], cex=cex);
	}
	points(snp+gap, mpt, pch=21, bg=pal[show+1], cex=cex*2);


	
	
	return(list(hapID.L=hap.id.l, hapCT.L=hap.count.l, hapID.R=hap.id.r, hapCT.R=hap.count.r));

}



#########################################################
# Make a list of all descenddant of a node
#########################################################

find.descendants<-function(node, nodes, descendant.list) {
	if (nodes[node,3]==0) {descendant.list<-c(descendant.list, node);}
	else {
		descendant.list<-find.descendants(nodes[node,3], nodes, descendant.list);
		descendant.list<-find.descendants(nodes[node,4], nodes, descendant.list);
	}

	return(descendant.list);
}



#########################################################
# Find time underneath a node
#########################################################

node.time<-function(nodes, node, time) {

	if (nodes[node,3]>0) {
		time<-time+nodes[node,2]-nodes[nodes[node,3],2];
		time<-node.time(nodes, nodes[node,3], time);
		time<-time+nodes[node,2]-nodes[nodes[node,4],2];
		time<-node.time(nodes, nodes[node,4], time);
	}
	return(time);
}


#########################################################
# Choose a node to mutate on the basis of a cumulative time
#########################################################

choose.node<-function(nodes, ctime, mrca) {

	parent<-match(c(1:(mrca-1)), nodes[,3:4])%%mrca;
	parent[parent==0]=mrca;
	dtime<-nodes[parent,2]-nodes[1:(mrca-1),2];

	for (node in 1:(mrca-1)) {
		ctime<-ctime-dtime[node];
		if (ctime<=0) return(node);
	}

	cat("\n\n*** Error: gone beyond end of nodes ***\n\n");

}




#########################################################
# Define order
#########################################################

order.seqs<-function(data, node=nrow(data), ord=c()) {

	if (data[node,3]==0) {
		ord[length(ord)+1]<-node;
	}
	else {
		ord<-order.seqs(data=data, node=data[node, 3], ord);
		ord<-order.seqs(data=data, node=data[node, 4], ord);
	}
	return(ord);
}



#########################################################
# To plot a tree vertically
# Data is 
# Col1 = name
# Col2 = Time
# COl3 = D1
# Col4 = No. mutations
#########################################################


plot.tree<-function(data, maxdepth=-1, add.mutations=FALSE, y.low=0, names.plot=TRUE, lwd=1, cex.mtn=0.5, plot.order=c(),cex.txt=0.5, srt.txt=90, col=c()) {

	nseq<-(nrow(data)+1)/2;
	mrca = nrow(data);
	if (length(plot.order)==0) plot.order<-order.seqs(data);

	if (maxdepth<0) maxdepth=max(data[,2]);
	plot(0,0, type="n", xaxt="n", xlab="", ylab="Time", bty="n", ylim=c(y.low,maxdepth), xlim=c(0,nseq+1));

	branch.pos<-order(plot.order);
	for (node in (nseq+1):mrca) {
		branch.pos[node]<-(branch.pos[data[node,3]]+branch.pos[data[node,4]])/2;
#DO horizontal line
		segments(x0=branch.pos[data[node,3]], y0=data[node,2], x1=branch.pos[data[node,4]], y1=data[node, 2], lwd=lwd);
#DO LH vertical
		segments(branch.pos[data[node,3]], data[data[node,3],2], branch.pos[data[node,3]], data[node, 2], lwd=lwd);
		if (add.mutations==T && data[data[node,3],5]>0) {
			time.span<-c(data[data[node,3],2], data[node,2]);
			y.vals<-time.span[1]+runif(data[data[node,3],5])*(time.span[2]-time.span[1]);
			points(rep(branch.pos[data[node,3]], data[data[node,3],5]), y.vals, pch=19, col=rgb(0,0,0.5), cex=cex.mtn);
		}

#DO RH vertical
		segments(branch.pos[data[node,4]], data[data[node,4],2], branch.pos[data[node,4]], data[node, 2], lwd=lwd);
		if (add.mutations==T && data[data[node,4],5]>0) {
			time.span<-c(data[data[node,4],2], data[node,2]);
			y.vals<-time.span[1]+runif(data[data[node,4],5])*(time.span[2]-time.span[1]);
			points(rep(branch.pos[data[node,4]], data[data[node,4],5]), y.vals, pch=19, col=rgb(0,0,0.5), cex=cex.mtn);
		}
	}

	if (names.plot==TRUE & length(col)==0) {
		text(x=c(1:nseq), y=rep(0, nseq), labels=data[plot.order,1], srt=srt.txt, cex=cex.txt, pos=1)
	}
	if (length(col)==nseq) {
		points(x=c(1:nseq), y=rep(0, nseq), pch=19, col=col[plot.order]);
	}
	return(plot.order)
}






#########################################################
# To plot a tree horizontally
# Data is 
# Col1 = name
# Col2 = Time
# COl3 = D1
# Col4 = No. mutations
#########################################################


plot.tree.horizontal<-function(data, maxdepth=-1, add.mutations=FALSE, y.low=0, names.plot=TRUE, lwd=1, cex.mtn=0.5) {

	nseq<-(nrow(data)+1)/2;
	mrca = nrow(data);
	plot.order<-order.seqs(data);

	if (maxdepth<0) maxdepth=max(data[,2]);
	plot(0,0, type="n", yaxt="n", xlab="Time", ylab="", bty="n", xlim=c(y.low,maxdepth), ylim=c(0,nseq+1));

	branch.pos<-order(plot.order);
	for (node in (nseq+1):mrca) {
		branch.pos[node]<-(branch.pos[data[node,3]]+branch.pos[data[node,4]])/2;
#DO horizontal line
		segments(y0=branch.pos[data[node,3]], x0=data[node,2], y1=branch.pos[data[node,4]], x1=data[node, 2], lwd=lwd);
#DO LH vertical
		segments(data[data[node,3],2], branch.pos[data[node,3]], data[node, 2], branch.pos[data[node,3]], lwd=lwd);
		if (add.mutations==T && data[data[node,3],5]>0) {
			time.span<-c(data[data[node,3],2], data[node,2]);
			x.vals<-time.span[1]+runif(data[data[node,3],5])*(time.span[2]-time.span[1]);
			points(y=rep(branch.pos[data[node,3]], data[data[node,3],5]), x=x.vals, pch=19, col=rgb(0,0,0.5), cex=cex.mtn);
		}

#DO RH vertical
		segments(data[data[node,4],2], branch.pos[data[node,4]], data[node, 2], branch.pos[data[node,4]], lwd=lwd);
		if (add.mutations==T && data[data[node,4],5]>0) {
			time.span<-c(data[data[node,4],2], data[node,2]);
			x.vals<-time.span[1]+runif(data[data[node,4],5])*(time.span[2]-time.span[1]);
			points(y=rep(branch.pos[data[node,4]], data[data[node,4],5]), x=x.vals, pch=19, col=rgb(0,0,0.5), cex=cex.mtn);
		}
	}

	if (names.plot==TRUE) {
		text(y=c(1:nseq), x=rep(0, nseq), labels=data[plot.order,1], srt=0, pos=2, cex=0.5)
	}
	return(plot.order)
}







#########################################################
# To simulate data under the coalescent
#########################################################

simulate.coalescent<-function(sample=10, theta=10, sites=10, rho=0, rmap=c(), nrun=1, 
	plot.tree=TRUE, return=TRUE, n.tree=5, all.seg=TRUE, inf.sites=TRUE) {

	cat("\n\n#####################################################\n");
	cat(" Coalescent simulation of data\n");
	cat(paste(" samples = ",sample,"\n",sep=""));
	cat("#####################################################\n\n");



	if (sample>1000) {
		cat("\n\nToo many samples: please use n<=100\n\n");
	}
	if (rho>100) {
		cat("\n\nToo much recombination: please use rho<=100\n\n");
	}
	if (length(rmap)>0) if (rmap[length(rmap)]>100) {
		cat("\n\nToo much recombination: please use rho<=100\n\n");
	}

	#No recombination
	if (length(rmap)==0 && rho<1e-6) {

	cat(paste("\n\nNo recombination: infinite-sites theta = ",theta,"\n",sep=""));
	nodes<-matrix(nrow=2*sample-1, ncol=5);
	colnames(nodes)<-c("Name", "Time", "D1", "D2", "Mutations");

	nodes[,1]<-c(1:nrow(nodes));
	nodes[,2:5]<-0;

	k<-sample;
	klist<-c(1:sample);
	time<-0;
	current.node<-k+1;

	while(k>1) {
		rate<-k*(k-1)/2;
		dt<--log(runif(1))/rate;
		time<-time+dt;
		l1<-ceiling(runif(1)*k);
		tmp<-klist[l1];
		klist[l1]<-klist[k];
		klist[k]<-tmp;
		l2<-ceiling(runif(1)*(k-1));
		nodes[current.node,2]<-time;
		nodes[current.node,3]<-klist[k];
		nodes[current.node,4]<-klist[l2];
		nodes[nodes[current.node,3],5]<-rpois(1, theta*(time-nodes[nodes[current.node,3],2])/2);
		nodes[nodes[current.node,4],5]<-rpois(1, theta*(time-nodes[nodes[current.node,4],2])/2);
		klist[l2]<-current.node;

		current.node<-current.node+1;
		k<-(k-1);
		klist<-klist[1:k];
	}

	n.mtns<-sum(nodes[,5]);
	mut.list<-runif(n.mtns);
	seqs<-matrix(0, nrow=sample, ncol=n.mtns);
	cum.mut<-1;
	for (node in 1:(nrow(nodes)-1)) if (nodes[node,5]>0) {
		who.list<-find.descendants(node, nodes, c());
		for (mut in 1:nodes[node,5]) {
			seqs[who.list,cum.mut]=1;
			cum.mut<-cum.mut+1;
		}
	}
	seqs<-seqs[,order(mut.list)];
	mut.list<-mut.list[order(mut.list)];
	

	if (plot.tree == TRUE) plot.tree(nodes, add.mutations=T);

	if (return==TRUE) return(list(tree=nodes, 
		data=list(type="haplotype", n.seq=sample, l.seq=n.mtns, pos=mut.list, names=c(1:sample), seq=seqs)
		
		));
	}


	#With recombination, use finfast type
	else {
	cat(paste("\n\nWith recombination: sites = ",sites,"\n",sep=""));
	if (length(rmap)>0) {cat(paste("Using Rmap: length = ",rmap[length(rmap)],"\n",sep=""));}
	else {cat(paste("Assuming constant rec rate: rho = ",rho,"\n",sep=""));}

	nodes<-array(dim=c(sites, 2*sample-1, 5)); 
	dimnames(nodes)<-list(c(1:sites),c(1:(2*sample-1)), 
		c("Name", "Time", "D1", "D2", "Mutations"));
	for (site in 1:sites) {
		nodes[site,,1]<-c(1:(2*sample-1));
		nodes[site,,2:5]<-0.0;
	}
	
	k<-sample;
	klist<-c(1:sample);
	anc<-matrix(ncol=sample, nrow=sites);
	for (i in 1:sites) anc[i,]<-c(1:sample);
	time<-0;
	if (length(rmap)==0) rmap<-c(0:(sites-1))*rho/(sites-1);
	current.nodes<-rep(sample+1,sites); #Ref to current node at each site
	rho.lin<-cbind(rep(1,sample), rep(sites, sample), rep(rmap[sites]-rmap[1],sample)); #Rho for each lineage
	colnames(rho.lin)=c("Lower", "Upper", "Rho");

	n.rec<-0;
	mrca<-2*sample-1;
	new.lin<-sample+1;

	while(k>1) {
		rho.tot<-sum(rho.lin[,3]);
		rate<-k*(k-1)/2+rho.tot/2;
		dt<--log(runif(1))/rate;
		time<-time+dt;

		#Recombination
		if (runif(1)<rho.tot/(rho.tot+k*(k-1))) {
			n.rec<-n.rec+1;
			cum.rec<-cumsum(rho.lin[,3]);
			tmp<-runif(1)*rho.tot;
			l1<-c(1:k)[cum.rec>tmp][1];

			cum.rec<-rmap[anc[,l1]>0]
			tmp<-runif(1)*rho.lin[l1,3];
			pos<-rho.lin[l1,1]+c(1:sites)[cum.rec-rmap[rho.lin[l1,1]]>tmp][1]-1;

#			cat(paste("\nRecombination on Lin = ", l1, " start = ", rho.lin[l1,1], 
#				" end = ", rho.lin[l1,2], " pos = ", pos, sep="")); 
			
			klist<-c(klist,new.lin);
			new.lin<-new.lin+1;
			anc<-cbind(anc, anc[,l1]);
			anc[(1:sites)>=pos,l1]=0
			anc[(1:sites)<pos,k+1]=0
			
			rho.lin<-rbind(rho.lin, rho.lin[l1,]);
			rho.lin[k+1,2]=rho.lin[l1,2];
			rho.lin[k+1,1]=min(c(1:sites)[anc[,k+1]>0]);
			rho.lin[k+1,3]=rmap[rho.lin[k+1,2]]-rmap[rho.lin[k+1,1]]
			rho.lin[l1,2]=max(c(1:sites)[anc[,l1]>0]);
			rho.lin[l1,3]=rmap[rho.lin[l1,2]]-rmap[rho.lin[l1,1]]

			k<-k+1;		
		}

		#Coalescent
		else {
			l1<-ceiling(runif(1)*k);
			#First move chosen lineage to end of klist, anc and rho.lin
			tmp<-klist[l1];
			klist[l1]<-klist[k];
			klist[k]<-tmp;
			tmp<-anc[,l1];
			anc[,l1]<-anc[,k];
			anc[,k]<-tmp;
			tmp<-rho.lin[l1,];
			rho.lin[l1,]<-rho.lin[k,];
			rho.lin[k,]<-tmp;
			l2<-ceiling(runif(1)*(k-1));

			new.nodes<-apply(anc[,c(k,l2)], 1, max);
			is.co<-c(1:sites)[apply(anc[,c(k,l2)], 1, min)>0];
			
			co.nodes<-current.nodes[is.co];
			for (co in 1:length(is.co)) {
				nodes[is.co[co],co.nodes[co],2]<-time;
				nodes[is.co[co],co.nodes[co],3]<-anc[is.co[co],k];
				nodes[is.co[co],co.nodes[co],4]<-anc[is.co[co],l2];
			}
			new.nodes[is.co]<-co.nodes;
			anc[,l2]<-new.nodes;
			current.nodes[is.co]<-co.nodes+1;

#			cat(paste("\nCoalescent between lineages ", l1, " and ", l2,sep=""));

			k<-k-1;
			klist<-klist[1:k];
			anc<-anc[,1:k];
			rho.lin<-rho.lin[1:k,];

			#Remove lineages at MRCA
			if(k>1) {
			anc[anc[,l2]==mrca,l2]<-0;
			if (sum(anc[,l2])>0) {
				rho.lin[l2,1]<-min(c(1:sites)[anc[,l2]>0]);
				rho.lin[l2,2]<-max(c(1:sites)[anc[,l2]>0]);
				if (rho.lin[l2,1]<rho.lin[l2,2]) {
					rho.lin[l2,3]<-rmap[rho.lin[l2,2]]-rmap[rho.lin[l2,1]];
				}
				else {rho.lin[l2,3]<-0;}
			}
			else {
				klist[k]<-klist[l2];
				anc[,l2]<-anc[,k];
				rho.lin[l2,]<-rho.lin[k,];
				k<-k-1;
				klist<-klist[1:k];
				anc<-anc[,1:k];
				rho.lin<-rho.lin[1:k,];
			}
			}
		}
	}

	cat(paste("\nTotal of ",n.rec," recombination events\n\n",sep=""));

	if (plot.tree == TRUE) {
		par(mfrow=c(1,n.tree));
		pos<-ceiling(c(1:n.tree)*sites/(n.tree+1));
		mx<-max(nodes[pos,,2]);
		for (tree in 1:n.tree) {
			plot.tree(nodes[pos[tree],,], maxdepth=mx, add.m=TRUE);
			title(main=paste("Position ", pos[tree], sep=""));
		}
		cat(paste("\nPrinting trees at position", pos, "\n", sep=" "));
		par(mfrow=c(1,1));
	}

	seqs<-matrix(0, nrow=sample, ncol=sites);
	is.seg<-c(1:sites);
	times.sites<-rep(0,sites);
	for (site in 1:sites) times.sites[site]<-node.time(nodes[site,,], mrca, 0);
	if (all.seg == FALSE) {
		for (site in 1:sites) {
			if (runif(1)<exp(-theta*times.sites[site]/(2*sites))) is.seg[site]<-0;
		}
		cat(paste("Total of ", sum(is.seg>0), " sites segregating\n", sep=""));
	}
	else {
		cat(paste("All sites (", sites, ") segregating\n", sep=""));
	}

	#Note will only throw down a max of 1 mtn per site
	for (site in 1:sites) if (is.seg[site]>0) {
		node<-choose.node(nodes[site,,], runif(1)*times.sites[site], mrca);
		who.list<-find.descendants(node, nodes[site,,], c());
		seqs[who.list,site]=1;
	}

	
	if (return==TRUE) return(list(tree=nodes, 
		data=list(type="haplotype", n.seq=sample, l.seq=sum(is.seg>0), 
		pos=c(1:sites)[is.seg>0], names=c(1:sample), seq=seqs[,is.seg>0])));
	}
	

}


#########################################################################
#To simulate new haplotypes from existing ones using LS method (approx)
#########################################################################

simulate.seqs.LS<-function(data, rho=0, rmap=c(), n.seq=1, theta=0.001) {

	h<-data$seq;

	if (length(rmap)==0) rmap=data$pos/(max(data$pos)-min(data$pos))*rho;
	rmap=1-exp(-diff(rmap)/data$n.seq);
	p.ident<-exp(-theta/data$n.seq);

	seqs.new<-matrix(nrow=n.seq, ncol=data$l.seq);
	h.new<-vector(length=data$l.seq);

	for (i in 1:n.seq) {
		states<-as.integer(runif(data$l.seq)*data$n.seq)+1;
		switch<-runif(data$l.seq-1)<=rmap;
		mutate<-runif(data$l.seq)>p.ident;
		h.new[1]<-states[1];
		for (j in 1:length(switch)) {
			if (switch[j]==TRUE) {h.new[j+1]<-states[j+1];}
			else {h.new[j+1]<-h.new[j];}
		}
		for (j in 1:data$l.seq) {
			if (mutate[j]==FALSE) {seqs.new[i,j]<-h[h.new[j],j];}
			else {seqs.new[i,j]<-1-h[h.new[j],j];}
		}
	}
	

	return(list(type="haplotype", n.seq=n.seq, l.seq=data$l.seq, 
		pos=data$pos, names=c(1:n.seq), seq=seqs.new));	

}




#########################################################
# Some example coalescent simulations
#########################################################

example.sims<-function(value=1, seed=235654) {

	set.seed(seed);

	if (value==1) {
		cat("\n\nData simulated with no recombination\n\n");
		s<-simulate.coalescent(sample=50, sites=50);
		summarise.data(s$data);
	}

	else if (value==2) {
		cat("\n\nData simulated with intermediate recombination\n\n");
		s<-simulate.coalescent(rho=10, sample=50, sites=50);
		summarise.data(s$data);
	}

	else if (value==3) {
		cat("\n\nData simulated with high recombination\n\n");
		s<-simulate.coalescent(rho=100, sample=50, sites=50);
		summarise.data(s$data);
	}

	else if (value==4) {
		cat("\n\nData simulated with a hotspot model\n\n");
		s<-simulate.coalescent(rmap=c(rep(0,25), rep(50,25)), sample=50, sites=50);
		summarise.data(s$data);
	}

	return(s);
}



###########################################################
# Count the number of haplotypes and haplotype homozygosity
###########################################################



count.haps<-function(data) {

	s<-data$seq;
	d<-dist(s, method="manhattan");
	hc<-hclust(d, method="average")$order;
	d<-as.matrix(d);
	d<-d[hc,hc];
	n.hap<-1;
	n.hap.i<-1;
	pairs.ident<-0;
	last.hap<-1;
	for (j in 2:nrow(d)) {
		if (d[last.hap,j]==0) {n.hap.i<-n.hap.i+1;}
		else {
			n.hap<-n.hap+1;
			pairs.ident<-pairs.ident+n.hap.i*(n.hap.i-1)/2;
			n.hap.i<-1;
			last.hap<-j;
		}
	}

	cat(paste("\n\nNumber of sequences = ", data$n.seq, sep=""));
	cat(paste("\nNumber of haplotypes = ", n.hap, sep=""));
	cat(paste("\n\nHaplotype homozygosity = ", pairs.ident*2/(data$n.seq*(data$n.seq-1)),
		"\n\n", sep=""));
}




#############################################################
# To simulate 2-population data from beta-binomial model
#############################################################

simulate.bb<-function(n.seq=c(10,10), l.seq=10, fst=0.01, f.min=0.01) {

	#Location of SNPs
	pos<-c(1:l.seq);

	#Allele frequencies in ancestor: drawn from 1/x distribution with min = f.min
	anc.freq<-exp(runif(l.seq)*log(f.min));
	alpha<-anc.freq/(fst*(1+anc.freq));
	beta<-alpha*(1-anc.freq)/anc.freq;

	#Allele frequencies in each population
	pop1.freq<-rbeta(l.seq, alpha, beta);
	pop2.freq<-rbeta(l.seq, alpha, beta);
	
	#Haplotypes for each population
	pop1.haps<-matrix(nrow=n.seq[1], ncol=l.seq);
	for (i in 1:n.seq[1]) pop1.haps[i,]<-rbinom(l.seq, 1, pop1.freq);
	pop2.haps<-matrix(nrow=n.seq[2], ncol=l.seq);
	for (i in 1:n.seq[2]) pop2.haps[i,]<-rbinom(l.seq, 1, pop2.freq);

	return(list(type="haplotype", n.seq=sum(n.seq), l.seq=l.seq, pos=pos, 
		names=c(1:sum(n.seq)), seq=rbind(pop1.haps, pop2.haps)));
}



############################################################################
# To simulate data under the coalescent with n popns and single split time
############################################################################

simulate.coalescent.split<-function(sample=10, theta=10, sites=10, rho=0, rmap=c(), nrun=1, 
	plot.tree=TRUE, return=TRUE, n.tree=5, all.seg=TRUE, inf.sites=TRUE,
	pop.assign=c(), split.time=0.0) {

	cat("\n\n#####################################################\n");
	cat(" Coalescent simulation of data\n");
	cat(paste(" samples = ",sample,"\n",sep=""));
	cat("#####################################################\n\n");



	if (sample>100) {
		cat("\n\nToo many samples: please use n<=100\n\n");
	}
	if (rho>100) {
		cat("\n\nToo much recombination: please use rho<=100\n\n");
	}
	if (length(rmap)>0) if (rmap[length(rmap)]>100) {
		cat("\n\nToo much recombination: please use rho<=100\n\n");
	}
	if (length(pop.assign)==0) {
		pop.assign<-rep(1, sample);
	}

	#No recombination
	if (length(rmap)==0 && rho<1e-6) {

	cat(paste("\n\nNo recombination: infinite-sites theta = ",theta,"\n",sep=""));
	nodes<-matrix(nrow=2*sample-1, ncol=5);
	colnames(nodes)<-c("Name", "Time", "D1", "D2", "Mutations");

	nodes[,1]<-c(1:nrow(nodes));
	nodes[1:sample,1]<-pop.assign;
	nodes[,2:5]<-0;

	k<-sample;
	klist<-c(1:sample);
	time<-0;
	current.node<-k+1;

	while(k>1) {
		rate<-k*(k-1)/2;
		dt<--log(runif(1))/rate;
		time<-time+dt;
		if (time>split.time) pop.assign<-rep(1, k);

		#Choose a first lineage and move to end of klist
		l1<-ceiling(runif(1)*k);
		tmp<-klist[l1];
		klist[l1]<-klist[k];
		klist[k]<-tmp;
		tmp<-pop.assign[l1];
		pop.assign[l1]<-pop.assign[k];
		pop.assign[k]<-tmp;

		#Choose a second lineage
		l2<-ceiling(runif(1)*(k-1));
		if (pop.assign[k]==pop.assign[l2]) {
			nodes[current.node,2]<-time;
			nodes[current.node,3]<-klist[k];
			nodes[current.node,4]<-klist[l2];
			nodes[nodes[current.node,3],5]<-rpois(1, theta*dt/2);
			nodes[nodes[current.node,4],5]<-rpois(1, theta*dt/2);
			klist[l2]<-current.node;

			current.node<-current.node+1;
			k<-(k-1);
			klist<-klist[1:k];
			pop.assign<-pop.assign[1:k];
		}
	}

	n.mtns<-sum(nodes[,5]);
	mut.list<-runif(n.mtns);
	seqs<-matrix(0, nrow=sample, ncol=n.mtns);
	cum.mut<-1;
	for (node in 1:(nrow(nodes)-1)) if (nodes[node,5]>0) {
		who.list<-find.descendants(node, nodes, c());
		for (mut in 1:nodes[node,5]) {
			seqs[who.list,cum.mut]=1;
			cum.mut<-cum.mut+1;
		}
	}
	seqs<-seqs[,order(mut.list)];
	mut.list<-mut.list[order(mut.list)];


	if (plot.tree == TRUE) plot.tree(nodes, add.mutations=T);

	if (return==TRUE) return(list(tree=nodes, 
		data=list(type="haplotype", n.seq=sample, l.seq=n.mtns, pos=mut.list, names=c(1:sample), seq=seqs)
		
		));
	}


	else {
		cat("\n\nNot yet imlpemented recombination and pop division\n\n");
	}


}



######################################################################
# To read in data dumped from HapMap web-site
######################################################################


read.hapmap.data<-function(file="dumped_region.txt", maf=FALSE, convert.to.binary=TRUE) {

	data<-read.table(file, as.is=T, fill=T, sep="-");
	l<-nrow(data);
	snp<-c(1:l)[data[,1]=="snps:"];
	hap<-c(1:l)[data[,1]=="phased_haplotypes:"];
	snps<-data[(snp+1):(hap-1),2];
	n.snp<-length(snps);
	pos<-c();
	snp.id<-c();
	snps<-strsplit(snps, ":");
	for (i in 1:n.snp) {
		pos[i]<-as.numeric(snps[[i]][2]);
		snp.id[i]<-strsplit(snps[[i]][1], " ")[[1]][2];
	}
	haps<-data[(hap+1):nrow(data),2];
	n.haps<-length(haps);
	haps<-strsplit(haps, ":");
	seqs.char<-c();
	names<-c();
	for (i in 1:n.haps) {
		names[i]<-haps[[i]][1];
		seqs.char[i]<-haps[[i]][2];
	}
	
	#Convert characters to values A=1, C=2, G=3, T=4
	seqs<-matrix(nrow=n.haps, ncol=n.snp);
	for (i in 1:n.haps) {
		hap<-strsplit(seqs.char[i], "")[[1]];
		hap<-hap[2:length(hap)];
		base<-c(1:n.snp)[hap=="A"];
		seqs[i,base]<-1;
		base<-c(1:n.snp)[hap=="C"];
		seqs[i,base]<-2;
		base<-c(1:n.snp)[hap=="G"];
		seqs[i,base]<-3;
		base<-c(1:n.snp)[hap=="T"];
		seqs[i,base]<-4;
	}

	#Convert to 0/1
	if (convert.to.binary==TRUE) {
		mn<-apply(seqs, 2, min);
		for (i in 1:n.snp) {
			seqs[,i]<-seqs[,i]==mn[i]
		}
	}

	#Convert 0 to major allele - best not to really
	if (maf==TRUE) {
		sm<-apply(seqs, 2, sum);
		for (i in 1:n.snp) {
			if (sm[i]>n.haps/2) seqs[,i]<-1-seqs[,i];
		}
	}


	data<-list(type="haplotype", n.seq=n.haps, l.seq=n.snp, pos=pos, snp.id=snp.id, names=names, seq=seqs);
	cat(paste("\n\nRead ", n.haps, " sequences of length ", n.snp, "\n\n", sep=""));
	
	return(data);

}




######################################################################
# To read in data from multiple populations dumped from HapMap web-site
######################################################################


read.combined.hapmap.data<-function(file.list=c("dumped_region.ceu.txt", "dumped_region.yri.txt","dumped_region.jpt+chb.txt"), convert.to.binary=TRUE) {

	d<-list();
	for (i in 1:length(file.list)) d[[i]]<-read.hapmap.data(file.list[i], convert=F);		
	n.snp<-d[[1]]$l.seq;
	for (i in 2:length(file.list)) if (d[[i]]$l.seq != n.snp) {
		cat("\n\nDifferent numbers of snps in data sets\n\n");
		return();
	}
	
	#Combine data
	seqs<-d[[1]]$seq;
	for (i in 2:length(file.list)) seqs<-rbind(seqs, d[[i]]$seq);

	#Convert to 0/1
	if (convert.to.binary==TRUE) {
		mn<-apply(seqs, 2, min);
		for (i in 1:n.snp) {
			seqs[,i]<-seqs[,i]==mn[i]
		}
	}

	#Combine names;
	names<-d[[1]]$names;
	tot.haps<-d[[1]]$n.seq;
	for (i in 2:length(file.list)) {
		names<-c(names, d[[i]]$names);
		tot.haps<-tot.haps+d[[i]]$n.seq;
	}

	return(list(type="haplotype", n.seq=tot.haps, l.seq=n.snp, pos=d[[1]]$pos, snp.id=d[[1]]$snp.id, names=names, seq=seqs));

}



######################################################################
# To calculate various selection statistics for a data set
######################################################################

selection.statistics<-function(data, derived=TRUE, rho=-1, fixed=FALSE, ehh.position=-1, f.cut=0.05, f.half=0.1) {

	ns <-data$n.seq;
	s <- data$l.seq;
	d <- dist(data$seq, method="manhattan");
	pwd = 2*sum(d)/(data$n.seq*(data$n.seq-1));
	f<-apply(data$seq, 2, sum);
	n.sing<-sum(f==1);
	if (derived == FALSE) n.sing<-n.sing + sum(f==(data$n.seq-1));

      con=vector(length=10);
      con[1]=sum(1/c(1:(ns-1)));
      con[2]=sum((1/c(1:(ns-1)))^2);
      con[3]=(ns+1)/(3*ns-3);
      con[4]=2*(ns*ns+ns+3)/(9*ns*ns-9*ns);
      con[5]=con[3]-1/con[1];
      con[6]=con[4]-(ns+2)/(ns*con[1])+con[2]/(con[1]*con[1]);
      con[7]=con[5]/con[1];
      con[8]=con[6]/(con[1]*con[1]+con[2]);
      con[9]=2*(ns*con[1]-2*(ns-1))/((ns-1)*(ns-2));
      con[10]=con[9]+(ns-2)/((ns-1)*(ns-1))+((2/(ns-1))*(1.5-(2*(con[1]+1/ns)-3)/(ns-2)-1/ns));
      con[11]=(ns*ns*con[2]/((ns-1)*(ns-1))+con[1]*con[1]*con[10]-2*ns*con[1]*(con[1]+1)/((ns-1)*(ns-1)))/(con[1]*con[1]+con[2]);       con[12]=ns/(ns-1)*(con[1]-ns/(ns-1))-con[11];

	if (derived==TRUE) {
		f<-snp.frequency(data, maf=FALSE)*data$n.seq;
		h<-hist(f[f>0 & f<data$n.seq], breaks=c(0:(data$n.seq-1)), plot=F)$counts;
		theta.der<-2*sum(h*c(1:(data$n.seq-1))^2)/(data$n.seq*(data$n.seq-1));
	}

	if (derived == FALSE) {
		return(list(n.seq=ns, s=s, pi=round(pwd,2), theta.Watterson=round(data$l.seq/con[1],2), 
			theta.singleton = n.sing,
			Tajima.D= round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2), 
			FuLi.D = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2)));
	}

	if (ehh.position>0 & fixed==FALSE) {
		pos<-which(data$pos==ehh.position);
		iehh<-calculate.iehh(data=data, snp.pos=pos, rho=rho, f.cut=f.cut, f.half=f.half);

		return(list(n.seq=ns, s=s, pi=round(pwd,2), theta.Watterson=round(data$l.seq/con[1],2),
			theta.derived = round(theta.der,2), theta.singleton = n.sing,
			Tajima.D= round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2), 
			FuLi.D = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2),
			FayWu.H = round(pwd,2)-round(theta.der, 2),
			ehh.position=ehh.position, ehh.freq = sum(data$seq[,pos])/data$n.seq, iehh=iehh));
	}

	if (ehh.position<1 | fixed==TRUE) {

		return(list(n.seq=ns, s=s, pi=round(pwd,2), theta.Watterson=round(data$l.seq/con[1],2), 
			theta.derived = round(theta.der,2), theta.singleton = n.sing,
			Tajima.D= round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2), 
			FuLi.D = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2),
			FayWu.H = round(pwd,2)-round(theta.der, 2)));

		
	}

}


#Calulate Integrated EHH and plot

calculate.iehh<-function(data, snp.pos=1, rho=-1, der.allele=1, anc.allele=0, f.cut=0.05, plot.order=c(), f.half=0.1) {

	if (length(rho)==1) {
		rmap=(data$pos-data$pos[1])/(data$pos[data$l.seq]-data$pos[1]);
		if (rho>0) rmap<-rmap*rho;
	}

	list.der<-which(data$seq[,snp.pos]==der.allele);
	seq.der<-data$seq[data$seq[,snp.pos]==der.allele,];
	f.der<-apply(seq.der, 2, mean);
	seq.anc<-data$seq[data$seq[,snp.pos]==anc.allele,];
	f.anc<-apply(seq.anc, 2, mean);

	#Set all derived mutations polymoprhic on der background only to 0
	#seq.der[,f.der>0 & f.der<1 & f.anc==0]<-0;

	#Indetify core set left and right from selected SNP
	core<-rep(0, data$l.seq);
	core[snp.pos]<-1;

	#Set up output sequences
	seq2<-data$seq;
	seq2[list.der, snp.pos]<-seq2[list.der, snp.pos]+1.5;
	iehh<-0;

	#Left
	which.core<-c(1:nrow(seq.der));
	pos<-snp.pos-1;
	while (length(which.core)>1 & pos>0) {
		f<-mean(seq.der[which.core,pos]);
		if (f<0.5) core[pos]<-0;
		if (f>=0.5) core[pos]<-1;
		which.core<-which.core[seq.der[which.core,pos]==core[pos]];
		if (length(which.core)/data$n.seq < f.cut) which.core<-c();
		seq2[list.der[which.core], pos]<-seq2[list.der[which.core], pos]+0.5;
		iehh<-iehh+length(which.core)*(rmap[pos+1]-rmap[pos]);
		pos<-pos-1;
	}

	#Right
	which.core<-c(1:nrow(seq.der));
	pos<-snp.pos+1;
	while (length(which.core)>1 & pos<=data$l.seq) {
		f<-mean(seq.der[which.core,pos]);
		if (f<0.5) core[pos]<-0;
		if (f>=0.5) core[pos]<-1;
		which.core<-which.core[seq.der[which.core,pos]==core[pos]];
		if (length(which.core)/data$n.seq < f.cut) which.core<-c();
		seq2[list.der[which.core], pos]<-seq2[list.der[which.core], pos]+0.5;
		iehh<-iehh+length(which.core)*(rmap[pos]-rmap[pos-1]);
		pos<-pos+1;
	}

	pal<-c(hsv(2/3,0.3,1), hsv(2/3,1,1), hsv(1/6,0.3,1), hsv(1/6,1,1), "red");
	breaks = c(-1.5,0.2,0.75,1.25,2, 3);
	if (length(plot.order)==0) {
		d<-make.weighted.dist(data, method="normal", f.half=f.half, position=data$pos[snp.pos]);
		plot.order<-hclust(d, method="average")$order
	}
      image(x=data$pos, z=t(seq2[plot.order,]), col=pal, breaks=breaks, 
			xlab="Position", ylab="Chromosome", main="EHH Haplotype plot", yaxt="n");

	return(iehh);	

}






#Summarise edited output from recmin file

summarise.recmin = function(file) {

  x = read.table(file);
  x = as.matrix(x);
  pos<-x[1,];
  x<-x[2:nrow(x),];
  di<-cbind(c(1:nrow(x)), c(1:nrow(x)));
  x[di]<-0;

  pal = c(rgb(0,0,(0:15)/15), heat.colors(12));

  cat(paste("Summarising data for ", file, "\n\n", sep=""));
  cat(paste("Minimum number of recombination events = ", x[1,ncol(x)],"\n",sep=""));

  x11();
  par(mfrow=c(2,1));
  breaks<-c(-10,0.5,1000000);
  image(x=pos, y=pos, z=x, col=c("white", "black"), breaks=breaks,
	main="Incompatibility matrix", xlab="Position", ylab="Position");
  image(x=pos, y=pos, z=x, main="Rh matrix", xlab="Position", ylab="Position", col=pal);

  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (j>i) {x[i,j] = x[i,j]/(pos[j]-pos[i]);}
      else {x[i,j]=0;}
    }
  }

  x11();
  par(mfrow=c(1,1));
  image(x=pos, y=pos, z=x, main="Scaled Rh matrix", xlab="Position", ylab="Position", col=pal, zlim=c(0,0.003));
  
  
}



###############################################################################
# To read the OXSTAT human genetic map from the web
# Can return the map for the whole chromosome (default), 
# a segment (using lims) or selected positions (using pos)
###############################################################################


read.genetic.map<-function(lims=c(), pos=c(), chromosome=1, ncbi.build=35) {

	if (ncbi.build==35) fname<-paste("http://www.stats.ox.ac.uk/~mcvean/OXSTAT/GeneticMap_b", 
		ncbi.build, "/genetic_map_chr", chromosome, ".txt.gz", sep="");
	if (ncbi.build==36) fname<-paste("http://www.stats.ox.ac.uk/~mcvean/OXSTAT/GeneticMap_b", 
		ncbi.build, "/genetic_map_chr", chromosome, "_b36.txt.gz", sep="");

	map<-scan(gzcon(url(fname)), sep=" ", skip=1);

	if (ncbi.build==35) {
		map<-matrix(map, nrow=length(map)/6, ncol=6, byrow=T);
		colnames(map)<-c("position_b35", "CEU_rate(cM/Mb)", "YRI_rate(cM/Mb)",
			"JPT+CHB_rate(cM/Mb)", "COMBINED_rate(cM/Mb)", "Genetic_Map(cM)");
		map<-map[,c(1,5,6)];
	}
	
	if (ncbi.build==36) {
		map<-matrix(map, nrow=length(map)/3, ncol=3, byrow=T);
		colnames(map)<-c("position_b36", "COMBINED_rate(cM/Mb)", "Genetic_Map(cM)");
	}
	

	if (length(lims)>0) {
		map<-map[map[,1]>=lims[1] & map[,1]<=lims[2],];
		if(nrow(map)==0) {
			cat("\n\n*** Warning: not data to return - check limits ***\n\n");
			return();
		}
		return(map);
	}

	if (length(pos)>0) {
		op<-approx(map[,1], map[,3], xout=pos);
		rts<-diff(op$y)*1000000/diff(op$x);
		rts<-c(rts, rts[length(rts)]);
		map<-cbind(pos, rts, op$y);
		if (ncbi.build==35) colnames(map)<-c("position_b35", "COMBINED_rate(cM/Mb)", "Genetic_Map(cM)");
		if (ncbi.build==36) colnames(map)<-c("position_b36", "COMBINED_rate(cM/Mb)", "Genetic_Map(cM)");
		return(map);
	}

	return(map);

}







###########################################################
# Routine to read in data from ms/msHOT
###########################################################


read.ms<-function(file) {

	data<-list();

	input<-scan(file, what="character");

	n.samp<-as.numeric(input[2]);
	n.rep<-as.numeric(input[3]);
	rep<-0;

	for (i in 1:length(input)) {

		if (input[i]=="//") {

			rep<-rep+1;
			l.seq<-as.numeric(input[i+2]);

			new.data<-list(type="haplotype", n.seq=n.samp, l.seq=l.seq, 
			pos=as.numeric(input[(i+4):(i+3+l.seq)]), names=c(1:n.samp), seq=matrix(0, nrow=n.samp, ncol=as.numeric(input[i+2])));

			i<-i+3+l.seq;
			for (j in 1:n.samp) {
				new.data$seq[j,]<-as.numeric(strsplit(input[i+j], split="")[[1]]=="1")
			}

			data[[rep]]<-new.data;
		}
	}

	return(data);
}



###########################################################
# Routine to convert data from ms/msHOT to LDhat format
###########################################################

ms2ldhat<-function(file, total.length=1) {

	input<-read.ms(file);

	for (i in 1:length(input)) {

		input[[i]]$pos<-round(input[[i]]$pos*total.length, 3);

		locs.file<-paste("locs_", i, ".txt", sep="");
		sites.file<-paste("sites_", i, ".txt", sep="");
		write.ldhat.format(input[[i]], seq.file=sites.file, loc.file=locs.file);
	}

}






###########################################################
#Convert hclust output into my standard tree format
###########################################################

hc2tree<-function(hc) {

	n<-length(hc$height)+1;
	tree<-matrix(0, nrow=2*n-1, ncol=5);
	colnames(tree)<-c("Node", "Height", "D1", "D2", "Mutations");
	tree[,1]<-c(1:nrow(tree));
	
	for (i in 1:nrow(hc$merge)) {
		which.node<-n+i;
		tree[which.node,2]<-hc$height[i];
		if (hc$merge[i,1]<0) {tree[which.node,3]<--hc$merge[i,1];}
		else {tree[which.node,3]<-hc$merge[i,1]+n;}
		if (hc$merge[i,2]<0) {tree[which.node,4]<--hc$merge[i,2];}
		else {tree[which.node,4]<-hc$merge[i,2]+n;}
	}

	return(tree);

}



###########################################################
#Script to read in data to R in FASTA format
###########################################################

read.fasta<-function(file, aligned=TRUE, gap.as.state=TRUE, biallelic.only=TRUE, recode.biallelic=TRUE) {

	allowed.bases<-c("T", "C", "A", "G", "t", "c", "a", "g", "U", "u", "?", "N", "n", "W", "w", "R", "r", "Y", "y", "M", "m", "K", "k", "S", "s", "B", "b", "D", "d", "H", "h", "V", "v", "-");
	bases.coding <-c(1, 2, 3, 4, 1, 2, 3, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5);

	input<-scan(file, what="char", comment="#", sep="\n");
	
	seq.starts<-grep(">", input);

	seqs<-list();
	for (i in 1:length(seq.starts)) {
		seqs[[i]]<-list();
		seqs[[i]]$name <-strsplit(input[seq.starts[i]], ">")[[1]][2];
		s<-c();
		if (i != length(seq.starts)) up<-seq.starts[i+1]-1;
		if (i==length(seq.starts)) up<-length(input);
		for (j in (seq.starts[i]+1):up) {
			s<-c(s, input[j]);
		}
		s<-unlist(strsplit(s, ""));
		s<-s[s %in% allowed.bases];
		s<-bases.coding[match(s, allowed.bases)];
		seqs[[i]]$seq<-s;
		seqs[[i]]$length <- length(seqs[[i]]$seq);
	}
	
	#If not aligned, just return;
	if (aligned==FALSE) return (seqs);

	#If aligned, check all same length
	for (i in 2:length(seqs)) {
		if (seqs[[i]]$length != seqs[[1]]$length) {
			cat("\n\n*** Error in file: sequences must be same length ***\n\n");
			return();
		}
	}

	#Now create numeric representation of data - first convert to 
	#0=?/N/ambig, 1=T/t, 2=C/c, 3=A/a, 4=G/g, 5=gap/-

	
	seqs.numeric<-matrix(0, nrow=length(seqs), ncol=seqs[[1]]$length);
	for (i in 1:length(seqs)) seqs.numeric[i,]<-seqs[[i]]$seq;

	#Now set things to missing data
	for (i in 1:ncol(seqs.numeric)) {
		is.na(seqs.numeric[seqs.numeric[,i]==0,i])<-TRUE;
		if (gap.as.state==FALSE) is.na(seqs.numeric[which(seqs.numeric[,i]==5),i])<-TRUE;
	}

	#Now count number of alleles
	allele.ct<-rep(0, ncol(seqs.numeric));
	for (i in 1:ncol(seqs.numeric)) {
		a<-unique(seqs.numeric[,i]);
		a<-a[is.na(a)==F]
		allele.ct[i]<-length(a);
	}

	if (biallelic.only==TRUE) which.keep<-which(allele.ct==2);
	if (biallelic.only==FALSE) which.keep<-which(allele.ct>1);
	pos<-which.keep;
	seqs.numeric<-seqs.numeric[,which.keep, drop=F];
	allele.ct<-allele.ct[which.keep];

	#select poly sites and re-code biallelic ones as 0/1 with 0=major allele;
	for (i in 1:length(pos)) if (allele.ct[i]==2 & recode.biallelic) {
		a.base<-max(seqs.numeric[,i], na.rm=T);
		a.alt<-min(seqs.numeric[,i], na.rm=T);
		f<-mean(seqs.numeric[,i]==a.base, na.rm=T);
		if (f>0.5) {
			seqs.numeric[,i]<-as.numeric(seqs.numeric[,i]==a.alt);
		}
		if (f<=0.5) {
			seqs.numeric[,i]<-as.numeric(seqs.numeric[,i]==a.base);
		}
	}
		
	names<-c();
	for (i in 1:length(seqs)) names[i]<-seqs[[i]]$name;

	
	

	#Finally make a big data type
	data<-list(type="haplotype", n.seq=nrow(seqs.numeric), l.seq=ncol(seqs.numeric), names=names, pos = pos, seqs=seqs.numeric);

	return(data);
	



	
}

################################################################
#Script to estimate LD statistic D from genotype data using EM
################################################################

#Work out if evidence for recombination for 2 SNPs with GT data

is.rec.gt<-function(gts) {

	gts<-gts[apply(is.na(gts), 1, sum)<1,];
	val<-gts[,1]*3+gts[,2]+1;
	counts<-hist(val, plot=F, breaks=c(0:9))$counts

	is.present<-rep(0,4);
	is.present[1]<-counts[1]+counts[2];
	is.present[2]<-counts[2]+counts[3]+counts[6];
	is.present[3]<-counts[4]+counts[7]+counts[8];
	is.present[4]<-counts[8]+counts[9];

	return(as.numeric(min(is.present)>0));
}

min.rec.gt<-function(data) {

	min.r<-0;
	min.i<-rep(0, data$l.seq);
	min.i[2]<-is.rec.gt(data$seq[,1:2]);
	
	for (i in 3:data$l.seq) {
		curr.max<-min.i[i-1];
		min.i[i]<-curr.max;
		j<-i-1;
		while(j>0 & min.i[j]==curr.max) {
			new.rec<-is.rec.gt(data$seq[,c(j,i)]);
			if (new.rec) {
				min.i[i]<-curr.max+1;
				j<-1;
			}
			j<-j-1;
			if (j==0) break();
		}
	}
	
	return(min.i);
}


#GTs in order 00, 01, 02, 10, 11, 12, 20, 21, 22
#HTs in order 00, 01, 10, 11

llk.gt.2locus<-function(gt.counts, ht.freq){

	llk<-1;

	gt.freq<-rep(0,9);
	gt.freq[1]<-ht.freq[1]^2;
	gt.freq[2]<-2*ht.freq[1]*ht.freq[2];
	gt.freq[3]<-ht.freq[2]^2;
	gt.freq[4]<-2*ht.freq[1]*ht.freq[3];
	gt.freq[5]<-2*ht.freq[1]*ht.freq[4]+2*ht.freq[2]*ht.freq[3];
	gt.freq[6]<-2*ht.freq[2]*ht.freq[4];
	gt.freq[7]<-ht.freq[3]^2;
	gt.freq[8]<-2*ht.freq[3]*ht.freq[4];
	gt.freq[9]<-ht.freq[4]^2;

	which.present<-which(gt.counts>0);
	if (min(gt.freq[which.present])<=0) is.na(llk)<-TRUE;
	
	if (!is.na(llk)) llk<-sum(gt.counts[which.present]*log(gt.freq[which.present]));

	return(llk);

}


estimateD.genotype.EM<-function(gts, tol=1e-4, max.iterations=100, return.freq=FALSE, verbose=FALSE){

	#Input genotypes are columns of 0/1/2/NA

	#Counts contains 00/01/02/10/11/12/20/21/22 genotypes
	#Freqs contains 00/01/10/11 haplotype frequencies

	gts<-gts[apply(is.na(gts), 1, sum)<1,];
	val<-gts[,1]*3+gts[,2]+1;
	counts<-hist(val, plot=F, breaks=c(0:9))$counts
	freq<-vector(length=4);
	freq[1]<-2*counts[1]+counts[2]+counts[4]+counts[5]*0.5;
	freq[2]<-2*counts[3]+counts[2]+counts[6]+counts[5]*0.5;
	freq[3]<-2*counts[7]+counts[4]+counts[8]+counts[5]*0.5;
	freq[4]<-2*counts[9]+counts[6]+counts[8]+counts[5]*0.5;
	freq.dh<-counts[5]/sum(counts);
	freq<-freq/sum(freq);
	D0<-freq[1]*freq[4]-freq[2]*freq[3];

	if(verbose) {
		cat(paste("\nInitial haplotype freqs = ", paste(freq, sep="\t", coll=""), sep=""));
		cat(paste("\nInitial log likelihood  = ", llk.gt.2locus(counts, freq), sep=""));
	}

	#Only perform EM if double hets exist
	if (freq.dh>1e-6) {
		del<-1;
		its<-0;
		while(del>tol & its<max.iterations) {
			rr<-freq[1]*freq[4]/(freq[1]*freq[4]+freq[2]*freq[3]);
			freq[1]<-2*counts[1]+counts[2]+counts[4]+counts[5]*rr;
			freq[2]<-2*counts[3]+counts[2]+counts[6]+counts[5]*(1-rr);
			freq[3]<-2*counts[7]+counts[4]+counts[8]+counts[5]*(1-rr);
			freq[4]<-2*counts[9]+counts[6]+counts[8]+counts[5]*rr;
			freq<-freq/sum(freq);
			D1<-freq[1]*freq[4]-freq[2]*freq[3];
			del<-abs(D1-D0);
			D0<-D1;
			its<-its+1;

			if(verbose) {
				cat(paste("\nIt ", its, " haplotype freqs = ", paste(freq, sep="\t", coll=""), sep=""));
				cat(paste("\nIt ", its, " log likelihood  = ", llk.gt.2locus(counts, freq), sep=""));
			}
		}
	}

	if (return.freq==TRUE) return(freq);
	

	return(list(D=D0, fi=(freq[3]+freq[4])/sum(freq), fj=(freq[2]+freq[4])/sum(freq)));
	

}





######################################################################
# To read in data directly from HapMap web-site
######################################################################


get.hapmap.data<-function(pop=c("CEU", "YRI"), chr=17, start=38349000, stop=38579000, maf=FALSE, convert.to.binary=TRUE, rel="phased_hapmap3", md.thresh=0.2) {

	url<-"http://www.hapmap.org/cgi-perl/";
	tmp<-pop[1];
	if (length(pop)>1) for (i in 2:length(pop)) {
		tmp<-paste(tmp, "+", pop[i], sep="");
	}
	file<-paste(url, rel, "?chr=chr", chr, "&pop=", tmp, "&start=", as.integer(start), "&stop=", as.integer(stop), "&ds=r2&out=txt", sep="");

	data<-scan(file, what="character");
	pos<-as.integer(data[grep("rs", data)+1]);
	snp.id<-data[grep("rs", data)];
	snp.id<-unlist(strsplit(snp.id, ":"));
	n.snp<-length(pos);
	names<-data[grep("NA", data)];
	names<-unlist(strsplit(names, ":"));
	seqs.char<-data[grep("NA", data)+1];
	n.haps<-length(seqs.char);

	#Convert characters to values A=1, C=2, G=3, T=4
	seqs<-matrix(nrow=n.haps, ncol=n.snp);
	for (i in 1:n.haps) {
		hap<-strsplit(seqs.char[i], "")[[1]];
		seqs[i,hap=="A"]<-1;
		seqs[i,hap=="C"]<-2;
		seqs[i,hap=="G"]<-3;
		seqs[i,hap=="T"]<-4;
		is.na(seqs[i,hap=="N"])<-TRUE;
	}

	#Remove inds with high rates NA
	md<-apply(is.na(seqs), 1, mean);
	which.keep<-which(md<=md.thresh);
	seqs<-seqs[which.keep,];
	n.haps<-length(which.keep);
	names<-names[which.keep];


	#Convert to 0/1/NA
	if (convert.to.binary==TRUE) {
		mn<-apply(seqs, 2, min, na.rm=T);
		for (i in 1:n.snp) {
			seqs[,i]<-seqs[,i]==mn[i]
		}
	}

	#Convert 0 to major allele - best not to really
	if (maf==TRUE) {
		sm<-apply(seqs, 2, sum, na.rm=T);
		for (i in 1:n.snp) {
			if (sm[i]>n.haps/2) seqs[,i]<-1-seqs[,i];
		}
	}


	data<-list(type="haplotype", n.seq=n.haps, l.seq=n.snp, pos=pos, snp.id=snp.id, names=names, seq=seqs);
	cat(paste("\n\nRead ", n.haps, " sequences of length ", n.snp, "\n\n", sep=""));
	
	return(data);

}

######################################################################
#Routine to compute HWE p-value for vector (L=3) of GT counts
######################################################################

hwe<-function(cts, verbose=TRUE, return=FALSE) {

	n<-sum(cts);

	if (n<2) {
		cat("\n\n*** Error: sample size < 2 ***\n\n");
		return();
	}

	f1<-(cts[2]+2*cts[3])/(2*n);
	g1<-cts/sum(cts);
	which.nz<-which(cts>0);
	l1<-sum(cts[which.nz]*log(g1[which.nz]));
	
	g0<-c((1-f1)^2, 2*f1*(1-f1), f1^2);
	l0<-sum(cts[which.nz]*log(g0[which.nz]));

	dl<-2*(l1-l0);
	p<-1-pchisq(dl, 1);

	if (verbose) cat("\nP-value for HWE = ", p, "\n\n");
	if (return) return(p);

}


#Function to take a set of sites that represent the switch points (pos) in a set of sites and return those 
#e.g. 1, 3, 9 in 1....10 gives c(1,4,5,6,7,8,9)

switch.seq<-function(pos, sites) {
	pos.n<-c(0, pos, sites);
	val<-(-1)^(1:(length(pos.n)-1));
	nn<-rep(val, times=diff(pos.n));
	return(which(nn<0));
}


######################################################################
# To simulate data under the coalescent - template switching version
# If data is serial need to put # samples at each time as a vector and
# Indicate times in 'times.sample'
# Most recent samples in first entry, times must be positive/incrementing.
######################################################################

simulate.coalescent.ts<-function(sample=c(10), theta=10, sites=20, rho=0, rmap=c(), 
	plot.tree=TRUE, return=TRUE, n.tree=5, all.seg=TRUE, inf.sites=TRUE, cx1=FALSE, finite=FALSE,
	times.sample=rep(0,length(sample)), verbose=FALSE, e.switch=5) {


	if (verbose) {
		cat("\n\n#####################################################\n");
		cat(" Coalescent simulation with template switching\n");
		cat(paste(" samples = ",sample,"\n",sep=""));
		cat("#####################################################\n\n");
	}
	
	#Checks
	fl<-1;
	if (sum(sample)>1000) {
		cat("\n\nToo many samples: please use n<=100\n\n");
		fl<-0;
	}
	if (rho>100) {
		cat("\n\nToo much recombination: please use rho<=100\n\n");
		fl<-0;
	}
	if (length(rmap)>0) if (rmap[length(rmap)]>100) {
		cat("\n\nToo much recombination: please use rho<=100\n\n");
		fl<-0;
	}
	if (sites<20) {
		cat("\n\nToo few sites: please use sites>=20\n\n");
		fl<-0;
	}
	if (length(sample) != length(times.sample)) {
		cat("\n\nIf use serial sample need to have same length vectors for sample and times.sample\n\n");
		fl<-0;
		if (times.sample[1] !=0 ) fl<-0;
	}


	if (fl==0) return();

	tot.sample<-sum(sample);

	nodes<-array(dim=c(sites, 2*tot.sample-1, 5)); #make array to store node info
	dimnames(nodes)<-list(c(1:sites),c(1:(2*tot.sample-1)), 
		c("Name", "Time", "D1", "D2", "Mutations"));
	ss<-1:tot.sample;
	ts.init<-rep(times.sample, times=sample);
	for (site in 1:sites) {
		nodes[site,,1]<-c(1:(2*tot.sample-1));
		nodes[site,,2:5]<-0.0;
		nodes[site,ss,2]<-ts.init;
	}	#put values in names, times and zeroes in rest	
	
	#Initialise
	k<-tot.sample;			#k = number of active lineages
	klist<-c(1:k);			# klist = list of lineages eg: [1] 1 2 3 4
	active<-rep(0, tot.sample);
	active[1:sample[1]]<-1;		#Define which lineages are active
	anc<-matrix(ncol=tot.sample, nrow=sites);	# make anc matrix - who is ancestral to whom
	for (i in 1:sites) anc[i,]<-c(1:tot.sample);	#fill anc matrix
	time<-0;				#initialise time=0
	if (length(rmap)==0) rmap<-c(0:(sites-1))*rho/(sites-1);  # make up rmap if none specified 
	r.prob<-diff(rmap);  id.sites<-1:(sites-1);
	current.nodes<-rep(tot.sample+1,sites); #Ref to current node at each site
	rho.lin<-cbind(rep(1,tot.sample), rep(sites, tot.sample), rep(rmap[sites]-rmap[1],tot.sample)); #Rho for each lineage
	colnames(rho.lin)=c("Lower", "Upper", "Rho");#name rho.lin cols

	n.rec<-0;			#Cumulative number of recombination events
	mrca<-2*tot.sample-1;	#Identifier for MRCA at each site
	new.lin<-tot.sample+1;	#keeps record of latest lineage
	epoch<-1;
	next.time<- -1;
	if (length(sample)>1) next.time<-times.sample[2];

	while(k>1) {

		w.active<-which(active==1);

		rho.tot<-sum(rho.lin[w.active,3]);	#total rho in sample. i.e k*rho!
		ke<-length(w.active);			#Allows 'inactive' lineages to occur
		rate<-ke*(ke-1)/2+rho.tot/2;		#rate of time to next event k(K – 1 + rho)/2
		if (rate>0) {
			dt<- -log(runif(1))/rate;		#inverse transform sampling from ~exp(rate)
			time<-time+dt;				#update time by adding time+dt
			fl<-1;
			ne<-as.integer(runif(1)<rho.tot/(rho.tot+ke*(ke-1)));  #'1' if REC event, 0 if CO event
		} else {
			time<-next.time;	#Allows for cases where nothing can happen in lineages left but there are ones yet to add
		}
		if (verbose) cat("\nCurrently ", ke, " active lineages and event rate ", rate, ": NE at ", time);

		#If serial and time exceeds next addition of samples - add in next tranche;
		if ((next.time>0) & (time>=next.time)) {

			epoch<-epoch+1;
			time<-next.time;
			to.add<-(sum(sample[1:(epoch-1)])+1):sum(sample[1:epoch]);
			active[klist %in% to.add]<-1;
			
			if (epoch >= length(sample)) {
				next.time<- -1;
			} else {
				next.time<-times.sample[epoch+1];
			}
			fl<-0;

			if (verbose) cat("\nNew epoch - adding samples ", to.add, " at time ", time);
		}

		#Recombination event happens with prob (rho/(k-1+rho))
		if (fl==1 & ne==1) {

			# *** Choose sequence to recombine ***
			n.rec<-n.rec+1;					# number of recombination events??
			l1<-sample(w.active, 1, prob=rho.lin[w.active,3]);

			# *** Choose position to recombine ***
			if(cx1==FALSE){
			pos=sort(sample(id.sites,rpois(1,e.switch-1)+1, prob=r.prob)) #this was the old implementation
			#pos=sort(sample(id.sites,rpois(1,(4*sites/9719))+1), prob=r.prob) #this is the new implementation... as #sites increase, the rate comes closer to 5/sequence
			sw=switch.seq(pos,sites)
			} else {
				pos=sample(id.sites,1, prob=r.prob);
				sw=switch.seq(pos,sites)		#Function that determines which sites will go into the recombinant
			}
				
			# *** Make recombinant****			
			klist<-c(klist,new.lin);		#add new lineage to klist
			new.lin<-new.lin+1;			#update new.lin for next round
			anc<-cbind(anc,anc[,l1]*0)		#add new hap to anc
			anc[sw,(k+1)]=anc[sw,l1]
			anc[sw,l1]=0	#populate anc with correct ancestral material from switching 
			
			#check that ancestral lineages actually have material on them - if not just remove
			r.flag<-1;
			#First - if L1 no longer has ANC material, put it in k+1 position
			if(sum(anc[,l1])==0){
				tmp=anc[,l1]
				anc[,l1]=anc[,k+1]
				anc[,k+1]=tmp
			}

			#Now check if k+1 position has ancestral material and if not remove			
			if(sum(anc[,k+1])==0){
				anc<-anc[,-(k+1)];		#remove faulty anc lineage
				new.lin<-new.lin-1;	#correct new.lin
				klist<-klist[-(k+1)];	#remove last klist entry
				n.rec<-n.rec-1;
				r.flag<-0;
			}
	
			# *** update lists *** - ONLY if k+1th element HAS ancestral material
			if (r.flag==1) {
				k<-k+1;
				rho.lin<-rbind(rho.lin, rho.lin[l1,]);
				anc.tmp<-which(anc[,k]>0);
				rho.lin[k,2]<-max(anc.tmp);#ensure correct upper value
				rho.lin[k,1]<-min(anc.tmp);
				rho.lin[k,3]<-rmap[rho.lin[k,2]]-rmap[rho.lin[k,1]];
				anc.tmp<-which(anc[,l1]>0);
				rho.lin[l1,2]<-max(anc.tmp);
				rho.lin[l1,1]<-min(anc.tmp);
				rho.lin[l1,3]<-rmap[rho.lin[l1,2]]-rmap[rho.lin[l1,1]];
				active<-c(active, 1);
			}
							
		}

		#Coalescent event happens with prob (k-1)/(k-1+rho)
		else if (fl==1) {

			l1<-sample(w.active, 2);	#choose lineages to coalesce
			l2<-min(l1);
			l1<-max(l1);

			#First move chosen lineage to end of klist, anc and rho.lin
			tmp<-klist[l1];
			klist[l1]<-klist[k];
			klist[k]<-tmp;
			tmp<-anc[,l1];
			anc[,l1]<-anc[,k];
			anc[,k]<-tmp;
			tmp<-rho.lin[l1,];
			rho.lin[l1,]<-rho.lin[k,];
			rho.lin[k,]<-tmp;
			tmp<-active[l1];
			active[l1]<-active[k];
			active[k]<-tmp;

			new.nodes<-apply(anc[,c(k,l2)], 1, max);#merge 
			is.co<-which(apply(anc[,c(k,l2)], 1, min)>0);
			
			co.nodes<-current.nodes[is.co];
			for (co in 1:length(is.co)) {
				nodes[is.co[co],co.nodes[co],2]<-time;
				nodes[is.co[co],co.nodes[co],3]<-anc[is.co[co],k];
				nodes[is.co[co],co.nodes[co],4]<-anc[is.co[co],l2];
			}
			new.nodes[is.co]<-co.nodes;
			anc[,l2]<-new.nodes;
			current.nodes[is.co]<-co.nodes+1;

			if (verbose) cat(paste("\nCoalescent between lineages ", l1, " and ", l2,sep=""));

			k<-k-1;
			keep<-1:k;
			klist<-klist[keep];
			anc<-anc[,keep];
			rho.lin<-rho.lin[keep,];
			active<-active[keep];

			#Remove lineages at MRCA
			if(k>1) {
				anc[anc[,l2]==mrca,l2]<-0;
				if (sum(anc[,l2])>0) {
					anc.tmp<-which(anc[,l2]>0);
					rho.lin[l2,1]<-min(anc.tmp);
					rho.lin[l2,2]<-max(anc.tmp);
					rho.lin[l2,3]<-rmap[rho.lin[l2,2]]-rmap[rho.lin[l2,1]];
				}
				else {
					klist[k]<-klist[l2];
					anc[,l2]<-anc[,k];
					rho.lin[l2,]<-rho.lin[k,];
					active[k]<-active[l2];
					k<-k-1;
					keep<-1:k;
					klist<-klist[keep];
					anc<-anc[,keep];
					rho.lin<-rho.lin[keep,];
					active<-active[keep];
				}
			}
		}
	}

	if (verbose) cat(paste("\n\nTotal of ",n.rec," recombination events\n\n",sep=""));

	if (plot.tree == TRUE) {
		par(mfrow=c(1,n.tree));
		pos<-ceiling(c(1:n.tree)*sites/(n.tree+1));
		mx<-max(nodes[pos,,2]);
		for (tree in 1:n.tree) {
			plot.tree(nodes[pos[tree],,], maxdepth=mx, add.m=TRUE);
			title(main=paste("Position ", pos[tree], sep=""));
		}
		if (verbose) cat(paste("\nPrinting trees at position", pos, "\n", sep=" "));
		par(mfrow=c(1,1));
	}

	seqs<-matrix(0, nrow=tot.sample, ncol=sites);
	is.seg<-c(1:sites);
	times.sites<-rep(0,sites);
	for (site in 1:sites) times.sites[site]<-node.time(nodes[site,,], mrca, 0);
	if (all.seg == FALSE) {
		for (site in 1:sites) {
			if (runif(1)<exp(-theta*times.sites[site]/(2*sites))) is.seg[site]<-0;
		}
		if (verbose) cat(paste("Total of ", sum(is.seg>0), " sites segregating\n", sep=""));
	} else {
		if (verbose) cat(paste("All sites (", sites, ") segregating\n", sep=""));
	}

	#infinite sites option
   if(!finite){
	for (site in 1:sites) if (is.seg[site]>0) {
		node<-choose.node(nodes[site,,], runif(1)*times.sites[site], mrca);
		who.list<-find.descendants(node, nodes[site,,], c());
		seqs[who.list,site]=1;
	}
   }	
   
   #else -> i.e if finite==true, then do finite sites option
   else{
   	for (site in 1:sites) if (is.seg[site]>0) {
        nmut<-rpois(1, times.sites[site]*theta/(2*sites));
        if (nmut>0) for (mut in 1:nmut) {
           node<-choose.node(nodes[site,,], runif(1)*times.sites[site], mrca);
           who.list<-find.descendants(node, nodes[site,,], c());
           seqs[who.list,site]=1-seqs[who.list,site];
        }
	if (nmut==0) is.seg[site]<-0;
    }
   	}
	
	w.seg<-which(is.seg>0);
	if (return==TRUE) return(list(tree=nodes, 
		data=list(type="haplotype", n.seq=tot.sample, l.seq=length(w.seg), 
		pos=w.seg, names=c(1:tot.sample), seq=seqs[,w.seg],n.rec=n.rec)));
	

}








#End of scripts





















