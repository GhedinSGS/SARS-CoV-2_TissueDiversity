# plot prep and filter lists:
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

d = c('delta1','delta2','delta3','delta4','delta5','mock')
isolate_colors = c('#5F50A2', '#9E1B45', '#F68D64','#e14d7d','#efa2bb','gray')
names(isolate_colors) = d
delta_colScale_fill = scale_fill_manual(name = "isolate",values = isolate_colors)
delta_colScale = scale_colour_manual(name = "isolate",values = isolate_colors)


ntlist = c('A','G','T','C','-')

aminoacids = c('G','A','L','M','F','W','K',
                'Q','E','S','P','V','I','C',
                'Y','H','R','N','D','T','*')


aacolors = as.vector(pals::alphabet(21))
names(aacolors) = aminoacids
aa_colScale_fill <- scale_fill_manual(name = "aa",values = aacolors)
aa_colScale <- scale_colour_manual(name = "aa",values = aacolors)


nsps = c('nsp1','nsp2', 'nsp3','nsp4','nsp5','nsp6',
                  'nsp7','nsp8','nsp9','nsp10','nsp11','nsp12a',
                  'nsp12b','nsp13','nsp14','nsp15','nsp16')

select_list = c('name', 'gene_id', 'ntpos','major','majorfreq','totalcount','aapos',
               'majoraa','majorcodon','refnt','refcodon','refaa','ay119_nt','ay119_aa',
                'sample_number','sample_details','plaque_pick','location')

merge_list = c('name', 'gene_id', 'ntpos','aapos',
               'refnt','refcodon','refaa','ay119_nt','ay119_aa',
               'sample_number','sample_details','plaque_pick','location')


# FUNCTIONS:
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


grabConsensus = function(vardf, ntpos_list, aminoacids, varaa_column){
    # grab the ntpos across all positions were a consensus change was identified
    vardf = vardf %>%
                filter(ntpos %in% ntpos_list) %>%
            unique() %>%
            as.data.frame() %>%
            rowwise() %>%
            mutate(check_syn = ifelse(majoraa %in% aminoacids & 
                                      .data[[varaa_column]] %in% aminoacids &
                                      .data[[varaa_column]] == majoraa,'syn',                      
                                 ifelse(majoraa %in% aminoacids & 
                                      .data[[varaa_column]] %in% aminoacids &
                                      .data[[varaa_column]] != majoraa,'nonsyn', NA))) %>%
            ungroup()
    
    return(vardf)
    
}

grabMin = function(df, ntpos_list){
    
     df = df %>%
            filter(ntpos %in% ntpos_list) %>%
            unique() %>%
            as.data.frame()
        
    return(df)
    
}


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


