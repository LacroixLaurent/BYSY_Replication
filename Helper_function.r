### This file contains home-made function used for the BYSY replication project.

### like sapply with mclapply
smclapply <- function(X, FUN, ...,
	mc.preschedule = TRUE, mc.set.seed = TRUE,
	mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
	mc.cleanup = TRUE, mc.allow.recursive = TRUE)
{
	require(parallel)
	simplify2array(mclapply(X, FUN, ...,
			mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
			mc.silent = mc.silent, mc.cores = mc.cores,
			mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive))}

### a function to change seqinf of a GRanges
NewSeqinfo <- function(GR,seqin) {
	seqlevels(GR,pruning.mode="coarse") <- seqlevels(seqin)
	seqinfo(GR) <- seqin
	return(GR)
}

### a function to plot forks from NFS data
plotforks <- function(toto,b2a.thr=0.02,plot.raw=F,Exp="EXP")
{

	suppressMessages(require(tidyverse))
	require(ggprism)
	require(gridExtra)
	theme_set(theme_bw())

	mypal=RColorBrewer::brewer.pal(12,"Paired")
	pl=list()
	for (i in 1:nrow(toto))
	{
		test <- toto %>% dplyr::slice(i)
		minpos <- min(test$signalr[[1]]$positions)
		maxpos <- max(test$signalr[[1]]$positions)
		if (plot.raw) {
			pl[[i]] <- ggplot(test$signalr[[1]]) +
				geom_point(aes(x=positions,y=Bprob,col="data.raw"),size=0.5,alpha=0.5,shape=16)+
#				geom_text(data=test$sl2[[1]],aes(x=sl.x,y=0,col="RDP_seg_type",label=sl.pat2,fontface="bold"))+
				geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
				geom_line(data=test$RDP[[1]],aes(x=x,y=y,col="RDP_segment"))+
				geom_hline(yintercept=b2a.thr,linetype="dashed") +
				geom_segment(data=test$forks[[1]],aes(x=X1,xend=X2,y=(0.7+sign(d.Y)/40),yend=(0.7+sign(d.Y)/40),col="NFS_fork_chase"),arrow=arrow(length = unit(0.1,"cm")))+
				geom_segment(data=test$forks[[1]],aes(x=X0,xend=X1,y=(0.7+sign(d.Y)/40),yend=(0.7+sign(d.Y)/40),col="NFS_fork_pulse"),arrow=arrow(length = unit(0.2,"cm")))+
				geom_text(data=test$forks[[1]],aes(x=(X0+X1)/2,y=(0.8+sign(d.Y)/20),col="NFS_fork_pulse",label=speed),size=4,show.legend=F)+
				xlab(paste0(test$chrom,"_",Exp," (kb)"))+
				scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(minpos,maxpos),
		breaks=seq(floor(minpos/10000)*10000,maxpos,10000),
		minor_breaks=seq(floor(minpos/5000)*5000,maxpos,5000),
		expand=c(0,0)
		)+
				guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,size=4,labels="")))+
				theme(legend.position = "right")+
				scale_color_manual(breaks = c("data.smoothed","data.raw","RDP_segment","RDP_seg_type","NFS_fork_pulse","NFS_fork_chase","NFS_speed","data.gap"),values = mypal[c(2,1,4,3,6,5,8,10)])+
				coord_cartesian(ylim=c(0,1))
		}else{
			pl[[i]] <- ggplot(test$signalr[[1]]) +
				geom_text(data=test$sl2[[1]],aes(x=sl.x,y=0,col="RDP_seg_type",label=sl.pat2,fontface="bold"), show.legend = F)+
				geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
				geom_line(data=test$RDP[[1]],aes(x=x,y=y,col="RDP_segment"))+
				geom_hline(yintercept=b2a.thr,linetype="dashed") +
				geom_segment(data=test$forks[[1]],aes(x=X1,xend=X2,y=(0.5+sign(d.Y)/40),yend=(0.5+sign(d.Y)/40),col="NFS_fork_chase"),arrow=arrow(length = unit(0.2,"cm")), show.legend = F)+
				geom_segment(data=test$forks[[1]],aes(x=X0,xend=X1,y=(0.5+sign(d.Y)/40),yend=(0.5+sign(d.Y)/40),col="NFS_fork_pulse"),arrow=arrow(length = unit(0.1,"cm")), show.legend = F)+
				geom_text(data=test$forks[[1]],aes(x=(X0+X1)/2,y=(0.8+sign(d.Y)/20),fontface="bold",col="NFS_speed",label=speed),size=2, show.legend = F)+
				xlab(paste0(test$chrom,"_",Exp," (kb)"))+
				scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(minpos,maxpos),
		breaks=seq(floor(minpos/100000)*100000,maxpos,100000),
		minor_breaks=seq(floor(minpos/50000)*50000,maxpos,50000),
		expand=c(0,0))+
				guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
				theme(legend.position = "right")+
				scale_color_manual(breaks = c("data.smoothed","data.raw","RDP_segment","RDP_seg_type","NFS_fork_pulse","NFS_fork_chase","NFS_speed","data.gap"),values = mypal[c(2,1,4,3,6,5,8,10)])+
				coord_cartesian(ylim=c(0,1))
		}
		if (test$gap_pos[[1]]$gap_start[1]>0)
		{
			pl[[i]] <- pl[[i]]+
				geom_segment(data=test$gap_pos[[1]],aes(x=gap_start,xend=gap_end,y=1,yend=1,col="data.gap"),size=4, show.legend = F)
		}
	}
	return(pl)
}

### a function to compute mean trace with padding for NFS data

padandtrace <- function(EXPforks,pad.max=500)
{

	suppressMessages(require(tidyverse))

	position <- (-10:(pad.max-1))*100
	EXPtrace0f <- EXPforks %>%
		select(trac,exp) %>%
		mutate(sig=map(trac, function(x) x%>% pull(signal))) %>%
		mutate(sig2=map(sig, function(x) {if(length(x)<pad.max){res=c(x,rep(NA,(pad.max+10-length(x))))}else{res=x[1:(pad.max+10)]};return(res)}))
	EXPtrace1mean <- colMeans(do.call(rbind,EXPtrace0f$sig2),na.rm=T)
	out <- suppressMessages(bind_cols(position,EXPtrace1mean))
	names(out) <- c("position","mean_trace")
	return(out)
}

compute_meantrace <- function(PLStib,trac.xmax=50000)
{
	toplot <- PLStib
	pad.max0 <- trac.xmax/100

	totrace <- split(toplot, toplot$exp)
	totrace2 <- lapply(totrace, function(x) padandtrace(x,pad.max0))
	totrace3 <- lapply(totrace2, function(x) do.call(bind_cols,x))
	totrace4 <- tibble(exp=factor(names(totrace3),levels=levels(toplot$exp)),trace=totrace3)
	return(totrace4)
}

plot_meantrace <- function(totrace4,explist,ymax0=0.8,xmax=50000,normalise=F,root_title="",pathout="",nameout="EXP_meantrace",expor=T)
{
	totrace1 <- totrace4 %>% filter(exp %in% explist)
	if (normalise) {
		totrace1 <- totrace1 %>%
			mutate(trace=map(trace, function(x) {x %>% mutate(mean_trace=mean_trace/max(mean_trace,na.rm=T))}))
		ymax=1
		lab.y="%BrdU (norm)"
	}else{
		ymax=ymax0
		lab.y="%BrdU"
	}

	totrace5 <- totrace1 %>% unnest(cols=c(trace))
	p1 <- ggplot(totrace5) +
		geom_line(aes(x=position, y=mean_trace,col=exp))+
		coord_cartesian(xlim=c(-1000,xmax),ylim=c(0,ymax))+
		paletteer::scale_color_paletteer_d("ggthemes::Classic_20")+
		ylab(lab.y)+
		xlab("Position (relative to X0, bp)")+
		ggtitle(paste0(root_title,"_",nameout))
	if(expor)
	{
		ggsave(paste0(pathout,root_title,"_",nameout,".pdf"), h=6,w=8)
	}else{p1}
}

### rescaling function
# this rescaling function put the signal to 0 for the quantile corresponding to infq and to 1 for the signal corresponding to supq.
myscaling0 <- function(x,infq=0.005,supq=0.995,...)
	{
		upper <- quantile(x,supq)
		lower <- quantile(x,infq)
		output <- (x-lower)/(upper-lower)
		return(output)
	}

