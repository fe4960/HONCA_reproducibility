# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

indata=read.txt('stdin', header=F)
names(indata)=c('x', 'y')

if (fixorder) {
	indata$x=factor(indata$x, levels=unique(indata$x))
}
print('===> indata')
print(indata)

colornames=c(
	brewer.pal(9, "Set1")
	, brewer.pal(8, "Set2")
	, brewer.pal(12, "Set3")
	)
labelx=levels(indata$x)
color_map=setNames(
  colornames[seq_along(labelx)],
  if (alphaorder) sort(labelx) else labelx
)
print('===> color_map')
print(color_map)

p=ggplot(data=indata, aes('', y, fill=x)) +
geom_bar(stat='identity', width=1) +
coord_polar('y', start=0) +
geom_text(aes(label=y), position=position_stack(vjust=0.5)) +
labs(
	x=NULL
	, y=NULL
	, fill=NULL
	, title=main
	) +
scale_fill_manual(values=color_map, name='') +
theme_classic() +
theme(
	axis.line=element_blank()
	, axis.text=element_blank()
	, axis.ticks=element_blank()
	, plot.title=element_text(hjust=0.5)
	)

# if (ncolor<=8) {
# 	p=p+scale_fill_brewer(palette='Set2')
# } else {
# 	p=p+scale_fill_manual(values=colornames[1:ncolor], name='')
# }

if (endsWith(outfile, '.png')) {
	ggsave(p, file=outfile, width=width, height=height, units='in', dpi=500)
} else {
	ggsave(p, file=outfile, width=width, height=height, useDingbats=F)
}
