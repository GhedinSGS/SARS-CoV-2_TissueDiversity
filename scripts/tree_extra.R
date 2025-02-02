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

locations = c('NY','FL','OTHER','p45')
loc_colors = c('#57c2ff','#ef9820','#a2a2a2','#83c03b')
names(loc_colors) = locations
loc_colScale_fill = scale_fill_manual(name = "location", values = loc_colors)
loc_colScale <- scale_colour_manual(name = "location", values = loc_colors)