# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(jlutils))

indata=readRDS(infile)
str(indata)

peakdata=indata@metadata[['PeakCallSummary']]
print('===> peakdata')
str(peakdata)
print(peakdata)
write.txt(peakdata, file=sprintf('%s/%s_peakcallsummary.txt.gz', outdir, bname))

if (length(group) > 0) {
	if (invert) {
		peakdata=peakdata[!(peakdata$Group %in% group), ]
	} else {
		peakdata=peakdata[peakdata$Group %in% group, ]
	}
	print('===> peakdata after group subsetting')
	str(peakdata)
	print(peakdata)
}

peakdata$Group=sub("\\s*\\([^)]+\\)", "", peakdata$Group)

print('===> trim () from peakdata$Group')
str(peakdata)
print(peakdata)

freq_sums=peakdata %>%
group_by(Group) %>%
summarise(Total_Freq=sum(Freq)*1000) %>%
ungroup() %>%
arrange(Total_Freq) # not arrange(desc(Total_Freq))

print('===> freq_sums')
str(freq_sums)
print(freq_sums)

peakdata=peakdata %>%
group_by(Group) %>%
mutate(Percentage=Freq/sum(Freq)) %>%
ungroup()

print('===> factor peakdata$Group, percentage')
str(peakdata)
print(peakdata)

peakdata=peakdata %>%
left_join(
	freq_sums %>% select(Group, Total_Freq),
	by="Group"
	) %>%
mutate(Group=paste0(Group, " (n=", format(Total_Freq, big.mark = ","), ")"))

freq_sums=freq_sums %>%
mutate(Group=paste0(Group, " (n=", format(Total_Freq, big.mark = ","), ")"))
peakdata$Group=factor(peakdata$Group, levels=freq_sums$Group)

print('===> add (n=total) to peakdata$Group')
str(peakdata)
print(peakdata)
write.txt(peakdata, file=sprintf('%s/%s_peakdata.txt.gz', outdir, bname))

colornames=c(
	brewer.pal(9, "Set1")
	, brewer.pal(8, "Set2")
	, brewer.pal(12, "Set3")
	)

if (mode=='proportion') {
	p=ggplot(data=peakdata, aes(x=Group, y=Percentage, fill=Var1, color=Var1)) +
	geom_bar(stat='identity', position='stack', show.legend=T, width=0.5) +
	xlab(xlab) +
	ylab(ylab) +
	ggtitle(main) +
	scale_colour_manual(values=colornames[1:nlevels(peakdata$Group)], name='') +
	scale_fill_manual(values=colornames[1:nlevels(peakdata$Group)], name='') +
	scale_y_continuous(label=scales::percent, breaks=scales::pretty_breaks(n=5)) +
	theme_bw() +
	theme(
		plot.background=element_blank()
		, panel.grid.major=element_blank()
		, panel.grid.minor=element_blank()
		, panel.border=element_blank()
		, plot.title=element_text(hjust=0.5)
		, axis.line=element_line(color='black')
		# , axis.text.x=element_text(vjust=0.5, hjust=0.5)
		, axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)
		, axis.text=element_text(color='black', size=size)
		, legend.position=legendposition
		, legend.key.size=unit(0.1, 'in')
		)
} else if (mode=='number') {
	p=ggplot(data=peakdata, aes(x=Group, y=Freq*1000, fill=Var1, color=Var1)) +
	geom_bar(stat='identity', position='stack', show.legend=T, width=0.5) +
	xlab(xlab) +
	ylab(ylab) +
	ggtitle(main) +
	scale_colour_manual(values=colornames[1:nlevels(peakdata$Group)], name='') +
	scale_fill_manual(values=colornames[1:nlevels(peakdata$Group)], name='') +
	scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
	theme_bw() +
	theme(
		plot.background=element_blank()
		, panel.grid.major=element_blank()
		, panel.grid.minor=element_blank()
		, panel.border=element_blank()
		, plot.title=element_text(hjust=0.5)
		, axis.line=element_line(color='black')
		# , axis.text.x=element_text(vjust=0.5, hjust=0.5)
		, axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)
		, axis.text=element_text(color='black', size=size)
		, legend.position=legendposition
		, legend.key.size=unit(0.1, 'in')
		)
}

outfile=sprintf('%s/%s.pdf', outdir, bname)
ggsave(p, file=outfile, width=width, height=height, useDingbats=F)
