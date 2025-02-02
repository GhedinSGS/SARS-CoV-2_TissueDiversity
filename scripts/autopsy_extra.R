# Functions and plot themes for autopsy SARS-CoV-2 analyses:

# lists: 
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

utrs = c('5\'UTR','3\'UTR')


nsps = c('nsp1', 'nsp2',
        'nsp3','nsp4', 'nsp5',
        'nsp6','nsp7','nsp8',
        'nsp9','nsp10','nsp11',
        'nsp12a','nsp12b','nsp13',
        'nsp14','nsp15','nsp16')

column_pull = c('name', 'general_location', 'detailed_location', 
                'lower_resp', 'predicted_time', 'resp', 
                'grouping', 'predicted_date', 'segment', 
                'gene_id', 'ntpos',  'totalcount', 
                'aapos',  'ay119_nt', 'ay119_aa','major')

# plots and colors:                 
PlotTheme1 = theme_bw() +
            theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              text = element_text(size = 16, color = 'black'))


col_list = c('black','#bdcee0','#dac8b9','#84bf3b',
                    '#66bace',
                    '#4372d6',
                    '#949cf3',
                    '#6233e3',
                    '#ba5fe5',
                    '#b03766',
                    '#db7032',
                    '#eece50',
                    '#cacaca',
            'black')
short_gene_list = c('5\'UTR','ORF1a','ORF1b','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10','3\'UTR')
names(col_list) = short_gene_list
gene_colScale_fill <- scale_fill_manual(name = "grp",values = col_list)
gene_colScale <- scale_colour_manual(name = "grp",values = col_list)



loc_colors = as.vector(c("#5a5156", "#939393", '#4475b4', "#fe7968",
                        "#a3cc7a", "#f1a5c9",'#ef9820', '#b3823e', "#ace7ea",
                        '#cd9fcb', '#ffd480', '#afe4c8'))
names(loc_colors) = orderit
loc_colScale_fill = scale_fill_manual(name = "location", values = loc_colors)
loc_colScale <- scale_colour_manual(name = "location", values = loc_colors)

aacolors = as.vector(pals::alphabet(21))
names(aacolors) = aminoacids
aa_colScale_fill <- scale_fill_manual(name = "aa",values = aacolors)
aa_colScale <- scale_colour_manual(name = "aa",values = aacolors)


gene_ord_list = c('5\'UTR',
                  'nsp1',
                  'nsp2',
                  'nsp3',
                  'nsp4',
                  'nsp5',
                  'nsp6',
                  'nsp7',
                  'nsp8',
                  'nsp9',
                  'nsp10',
                  'nsp11',

                  'nsp12a',
                  'nsp12b',
                  'nsp13',
                  'nsp14',
                  'nsp15',
                  'nsp16',

                  'S',
                  'ORF3a',
                  'E',
                  'M',
                  'ORF6',
                  'ORF7a',
                  'ORF7b',
                  'ORF8',
                  'N',
                  'ORF10',
                  '3\'UTR',
                'INTERGENIC')

cds_regions = c('nsp1','nsp2','nsp3',
                  'nsp4','nsp5','nsp6','nsp7',
                  'nsp8','nsp9','nsp10','nsp11',
                  'nsp12a','nsp12b','nsp13','nsp14',
                  'nsp15','nsp16','S','ORF3a',
                  'E','M','ORF6','ORF7a',
                  'ORF7b','ORF8','N','ORF10')

gene_color_list = c('#cacaca',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',
                    '#bdcee0',


                    '#dac8b9',
                    '#dac8b9',
                    '#dac8b9',
                    '#dac8b9',
                    '#dac8b9',
                    '#dac8b9',


                    '#84bf3b',
                    '#66bace',
                    '#4372d6',
                    '#949cf3',
                    '#6233e3',
                    '#ba5fe5',
                    '#b03766',
                    '#db7032',
                    '#eece50',
                    '#cacaca',
                    '#cacaca',
                  '#cacaca')


detailed_orderlist = c(
"appendix", 
"proximal-trachea",
"jejunum",
"L-bronchus",
"L-ventricle",
"L-inferior-lobe",
"R-bronchus",
"pericardium",
"L-superior-lobe",
"R-middle-lobe",
"rib",
"R-inferior-lobe",
"L-eye-optic-nerve",
"R-eye-retina",
"R-eye-choroid-sclera",
"nasal-placode",    
"sinus-turbinate",
"R-eye-cornea",
"tongue",
"R-kidney")


names(gene_color_list) = gene_ord_list

gene_colScale_fill <- scale_fill_manual(name = "grp",values = gene_color_list)

gene_colScale <- scale_colour_manual(name = "grp",values = gene_color_list)



# functions: 
buildAA = function(fasta, var_colname){
    
    # read aa fasta and set as a df
    myProtein = readAAStringSet(file = fasta) %>% as.data.frame() 
    
    # grab name
    n = gsub(" ", "", as.character(rownames(myProtein)[1])) 

    # make into a df
    df = cbind(as.vector(str_split_fixed(myProtein, pattern = "", n = nchar(myProtein))),
              seq(1, nchar(myProtein), by = 1),
               n) %>%
          as.data.frame()

    # adjust colnames
    colnames(df) = c(var_colname,'aapos','gene_id')
    
    return(df)
}

buildNT = function(fasta, var_colname){
    # read in sequence
    mySeq = as.character(readDNAStringSet(file = fasta)[[1]])

    # bind to df w/ ntpos info
    df = cbind(as.vector(str_split_fixed(mySeq, pattern = "", n = nchar(mySeq))),
      seq(1, nchar(mySeq), by = 1)) %>%
      as.data.frame()

    # adjust colnames
    colnames(df) = c(var_colname,'ntpos')

    # make sure correct case
    df[var_colname] = toupper(df[[var_colname]])

    return(df)
}

appendVoc = function(voc_aa, voc_nt, varfile){
    vardf = read.csv(varfile, header = TRUE, na.strings=c("","NA")) 
        
    vardf = merge(vardf, voc_nt, by = c('ntpos'), all.x = TRUE)
        
    vardf = merge(vardf, voc_aa, by = c('gene_id','aapos'), all.x = TRUE)
    
    return(vardf)
}


returnConList = function(vardf, var_column, coverage_cutoff=5, frequency_cutoff=0.01){
    # pull positions that differ from delta - return the list for each sample
    
    t1 = nrow(vardf)
    
    vardf = vardf %>%
                filter(.data[[var_column]] != major &
                          totalcount >= coverage_cutoff &
                          majorfreq >= frequency_cutoff | # added on 20230613
                          .data[[var_column]] != refnt &
                          totalcount >= coverage_cutoff &
                          majorfreq >= frequency_cutoff ) %>%
            unique() %>%
            as.data.frame()
    
    varlist = c(levels(factor(vardf$ntpos)))
    
    return(varlist)
    
}

returnMinList = function(df, ntlist, coverage_cutoff = 50, frequency_cutoff = 0.01){
    # not yet filtering for differences to the refnt we want to compare to
    df = df %>%
                filter(minor %in% ntlist & 
                       totalcount >= coverage_cutoff &
                       minorfreq >= frequency_cutoff) %>%
            unique() %>%
            as.data.frame()
    
    varlist = c(levels(factor(df$ntpos)))
    
    return(varlist)
}


GrabCoverage = function(df, cov_requirement, percent_requirement){
    total_nt = df %>% select(ntpos) %>% unique() %>% nrow()
    pass_nt = df %>% filter(totalcount >= cov_requirement) %>% select(ntpos) %>% unique() %>% nrow()
    temp = df %>% select(name) %>% unique() %>%
                mutate(genome_size = total_nt, 
                       pass_size = pass_nt,
                       percentage = pass_nt/total_nt,
                       rd_requirement = cov_requirement,
                       pass = ifelse(percentage >= percent_requirement, 'PASS','FAIL')
                       )
    return(temp)
}

grabPositions = function(df, ntpos_list){
    
     df = df %>%
            filter(ntpos %in% ntpos_list) %>%
            unique() %>%
            as.data.frame()
        
    return(df)
}


    