PlotTheme = theme_bw() +
            theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              text = element_text(size = 16, color = 'black'))

theme_set(PlotTheme)

orderit = c('nasal', 'oropharynx', 'trachea','lung','cardiovascular',
 'lymphoid', 'gastrointestinal', 'kidney', 'ocular','rib','CNS', 'ref')


aminoacids = c('G','A','L','M','F','W','K',
                'Q','E','S','P','V','I','C',
                'Y','H','R','N','D','T','*')


aacolors = as.vector(pals::alphabet(21))
names(aacolors) = aminoacids
aa_colScale_fill <- scale_fill_manual(name = "aa",values = aacolors)
aa_colScale <- scale_colour_manual(name = "aa",values = aacolors)


voc_cols = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2","#D55E00","#CC79A7","#999999","#000000")

haplotype_cols = c("#9e0142","#3288bd","#ffffbf","#cbc3e3","#5e4fa2",
"#d53e4f","#fdae61","#66c2a5","#fee08b","#abdda4","#f46d43")

ntlist = c('A','T','G','C')
