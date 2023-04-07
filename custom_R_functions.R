##### ============================================= #####
#####				Custom R functions				#####
#####		(c)	Grigorios Georgolopoulos 2018-2022	#####
#####				georgolog@gmail.com				#####
##### ============================================= #####

#!/usr/bin/env Rscript

##### Split vector into specified chunks or lines; equivalent to unix split #####
splt = function(x, chunks = NA, length = NA) {
    if (is.na(chunks) & is.na(length) | chunks < 1 & length < 1 | is.na(chunks) & length < 1 | is.na(length) & chunks < 1 ) {
        stop('Number of chunks or lines to split to not specified')
    } else if (!is.na(chunks) & !is.na(length)) {
        stop('Specify either number of chunks or number of lines')
    } else if (is.numeric(chunks) & is.na(length)) {
        return(split(x, sort((1:length(x))%%chunks)))
    } else if (is.na(chunks) & is.numeric(length)) {
        return(split(x, ceiling(seq_along(x)/length)))
    }
}

##### Calculate Cosine Distance #####
cos.sim <- function(X, index) 
{
    ix = index
    A = X[ix[1],]
    B = X[ix[2],]
    return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}

##### Get the last element of a vector #####
last <- function(x) {return(x[length(x)])} 

##### Compute Jaccard Similarity Index #####
jacc = function(x, y) {
    S = (x + y)
    if(sum(x) == 0 & sum(y) == 0) {
        return(NA)
    } else {
        n00 = sum(S==0)
        n10 = sum(x==1 & y == 0)
        n01 = sum(y==1 & x == 0)
        n11 = sum(S == 2)
        return(n11/sum(n01,n10, n11))
    }
}

##### Calculate Counts Per Million (CPM) #####
cpm <- function (expr_mat) 
{
    norm_factor <- colSums(expr_mat)
    return(t(t(expr_mat)/norm_factor) * 10^6)
}

##### "ggridis" lol! custom color palette #####
ggiridis <- function(x) {
    colorRampPalette(c('#92BD9B','#219066','#3E5D7B','#420555','#BE5337','#DCC22F','#FDEE98'))(x)
}

##### Calculate Geometric Mean #####
geo_mean = function(x, na.rm=TRUE, zero.propagate = FALSE) {
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

##### Get density of points in 2 dimensions.
# param x A numeric vector.
# param y A numeric vector.
# param n Create a square n by n grid to compute density.
# return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

##### Calculate the midpoints of groups of values #####
# This function is helpful to plot cluster labels
# in heatmaps in the middle of each cluster

mids = function(x) {
    s = sort(x)
    w = which(!duplicated(s))
    t = unlist(as.numeric(table(s)))
    m = floor(t/2)
    return(w+m)
}

##### Emulate ggplot Color Palette #####

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
} 

##### Shuffle a vector or a matrix row-wise or column-wise #####
# options:
#	mar = Margin; 1 for shuffling rows-wise, 2 for column-wise
#	fixed = if (mar > 1) to fix each margin or shuffle rows (or columns) independently
shuff = function(x, mar = 1, fixed = T) {
    if (length(dim(x)) > 2) {
        stop('x must be vector or matrix of 2 dimensions')
    } else if (length(dim(x)) == 1) {
        l = lenth(x)
        return(sample(x))
    } else if (length(dim(x)) == 2) {
        if (fixed) {
            if (!mar %in% 1:2) {
                stop('Please enter 1 for rows, 2 for columns')
            } else if (mar == 1) {
                v = 1:(dim(x)[[1]])
                return(x[sample(v),])
            } else if (mar == 2) {
                v = 1:(dim(x)[[2]])
                return(x[,sample(v)])
            } 
        } else if (!fixed) {
            if (!mar %in% 1:2) {
                stop('Please enter 1 for rows, 2 for columns')
                } else if (mar == 1) {
                    return(apply(x, MARGIN = 2, FUN = sample))
                } else if (mar == 2) {
                    return(t(apply(x, MARGIN = 1, FUN = sample)))
                } 
        } 
    }
}



##### Map values to colors #####

val2col <- function(x,pal,limits=NULL){

    if(is.null(limits)) limits=range(x)

    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]

}



##### Variance Normalization for sequencing data based on Satija et al. 2015 #####

normalize_variance = function(x, cut.along = NULL, nbin = 20, log = TRUE) {
    
    if (min(x) < 0) {
        
        stop('Matrix must contain non-negative values')
        
    } else {
        
        if (log != TRUE) {
            
            c = min(x[x!=0])
            
            X = log10(x + c*1.1)
        
        } else {
            
            X = x
            
        }
        
        if (is.null(cut.along)) {
            
            mu = rowMeans(X)
            
        } else if (length(cut.along) != dim(x)[[1]]) {
        
        stop('Values provided must be the same length as dim(x)[1]') 
            
            } else {
            
            mu = cut.along
        }
        
        names(mu) = 1:nrow(X)
        
        brk = hist(mu, nbin, plot = F)$breaks
        
        int = findInterval(mu, brk)
        
        bins = unique(sort((int)))
        
        if (any(table(int) == 1)) {
            
            tab = table(int)
            
            nm = names(tab[tab == 1])
            
            int[int == nm] = bins[length(bins)-1]
            
            bins = unique(sort((int)))
            
        } else {
            
            NULL
        }
        
        z.var = lapply(bins, FUN = function(i) {
            
            return(scale(rowVars(X[int %in% i,]), center = T, scale = T))
            
        })
        
        return(unlist(z.var))
    }
       
}

##### Revrse-compliment a sequence #####

rev_com = function(x, rc = TRUE, r = FALSE) {
tab = c('A','T','G','C','N'); names(tab)=c('T','A','C','G','N')

if(!is.character(x)) {
	
	stop('Please enter a valid sequence')

}

x = toupper(x)

x_splt = unlist(strsplit(x, split = ''))

if (!all(x_splt %in% tab)) {

	stop('Please enter a valid sequence')

} else {

rv_splt = x_splt[length(x_splt):1]

rc_splt = names(tab)[match(rv_splt, tab)]

rc_out = paste(rc_splt, collapse='')

if (rc == TRUE) {

	return(rc_out)

	} else if (r == TRUE & rc == FALSE)

		return(paste(rv_splt, collapse = ''))

	}
}

##### Write BED formatted tables #####
write.bed = function(x, file) {

write.table(x, file = file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

}

##### Shortcut to Standard Error of the Mean function #####

stderr = function(x) {

sd(x)/sqrt(length(x))

}

##### Matplotlib colormaps #####

#matplotlib_camp = function(pal, n) {
#	
#	pals = list( tab10 = c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'))
#
#	return(pals[[pal]][1:n])
#
#}
matplotlib_cmap = function(n, pal) {

pals = list(tab10 = c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'),
			tab20 = c('#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2','#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5'),
			tab20b = c('#393b79','#5254a3','#6b6ecf','#9c9ede','#637939','#8ca252','#b5cf6b','#cedb9c','#8c6d31','#bd9e39','#e7ba52','#e7cb94','#843c39','#ad494a','#d6616b','#e7969c','#7b4173','#a55194','#ce6dbd','#de9ed6'),
			tab20c = c('#3182bd','#6baed6','#9ecae1','#c6dbef','#e6550d','#fd8d3c','#fdae6b','#fdd0a2','#31a354','#74c476','#a1d99b','#c7e9c0','#756bb1','#9e9ac8','#bcbddc','#dadaeb','#636363','#969696','#bdbdbd','#d9d9d9'),
			tab40=c('#3182bd','#6baed6','#9ecae1','#c6dbef','#e6550d','#fd8d3c','#fdae6b','#fdd0a2','#31a354','#74c476','#a1d99b','#c7e9c0','#756bb1','#9e9ac8','#bcbddc','#dadaeb','#636363','#969696','#bdbdbd','#d9d9d9',
'#393b79','#5254a3','#6b6ecf','#9c9ede','#637939','#8ca252','#b5cf6b','#cedb9c','#8c6d31','#bd9e39','#e7ba52','#e7cb94','#843c39','#ad494a','#d6616b','#e7969c','#7b4173','#a55194','#ce6dbd','#de9ed6'),
			gnuplot = c('#000000','#2f0037','#43006b','#52009a','#5f01c2','#6a01e1','#7402f6','#7d04ff','#8605fc','#8e08ed','#960ad3','#9d0eaf','#a41283','#ab1751','#b11d1c','#b72300','#bd2b00','#c33300','#c93d00','#ce4800','#d45400','#d96100','#de6f00','#e37f00','#e89100','#eda300','#f1b800','#f6ce00','#fbe600','#ffff00'))
 if (!pal %in% names(pals)) {

 stop(paste('Please choose one of the valid color palettes:', paste(names(pals), collapse = ', ')))

 }

 if (pal == 'gnuplot') {

 	return(colorRampPalette(pals[[pal]])(n))
 } else {

 	if (n > length(pals[[pal]])) {

		n = length(pals[[pal]])
		
		print(paste('Returning maximum of', length(pals[[pal]]), 'colors'))
		
		return(pals[[pal]][1:n])
		
		}
		
		return(pals[[pal]][1:n])
	}

}

##### Calculate Robust Z-score #####
zscore.rob = function(x) {
    
    (x-median(x))/mad(x)
    
}

##### Convert between Gene ID formats #####

get_gene_names= function(keys, from_type, to_type, format = NA, organism,fill.na = TRUE) {

organism.str <- organism
organism <- eval(sym(organism.str))

library(organism.str, character.only = TRUE)

conversion = mapIds(organism, keys = keys, keytype = from_type, column=to_type)

if (fill.na) {
	conversion[is.na(conversion)] = names(conversion[is.na(conversion)])
}

if (format == 'data.frame') {
	conversion_df = data.frame(V1=names(conversion), V2=as.character(conversion), stringsAsFactors=FALSE)
	colnames(conversion_df) = c(from_type, to_type)
	return(conversion_df)
} else {

	return(conversion)
}
}

##### Round a value to the next decade or decimal, up or down #####

roundup = function(n, scale = 10, round = 'up') {
  if (!is.numeric(n)) {
    stop('n must be a number')
  } else {
    nn = n/scale
    if (round == 'up') {
      nn = ceiling(nn)*scale
      return(nn)
    } else if (round == 'down') {
      nn = floor(nn)*scale
      return(nn)
    } else {
      stop('round must be one of "up" ord "down"')
    }
  }
}

##### Expanded brewer.pal "Spectral" color palette ######

Spectral2 = c('#000000','#181346','#20417f','#3288BD','#66C2A5','#ABDDA4','#E6F598','#FFFFBF','#FEE08B','#FDAE61','#F46D43','#D53E4F','#481211','#000000','#481211','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#ffffff')
Clima = c("#68C3A5",
"#ADD8A4",
"#3289BD",
"#FBF7C1",
"#FCE08A",
"#E4EB99",
"#F36D43",
"#D63F4F",
"#20417F",
"#272363",
"#3E2E75",
"#683996",
"#91449A",
"#BD5BA3",
"#DA8EBD",
"#ECC3DC",
"#FFFFFF",
"#A1242E",
"#701111",
"#FBAE60")

#### Clean column names -- Borrowed from `janitor` package probably ####
clean_names <- function(.data, unique = FALSE) {
  n <- if (is.data.frame(.data)) colnames(.data) else .data
  n <- gsub("%+", "_pct_", n)
  n <- gsub("\\$+", "_dollars_", n)
  n <- gsub("\\++", "_plus_", n)
  n <- gsub("\\*+", "_star_", n)
  n <- gsub("\\-+", ".", n)
  n <- gsub("#+", "_n_", n)
  n <- gsub("&+", "_and_", n)
  n <- gsub("@+", "_at_", n)
  n <- gsub("[^a-zA-Z0-9_]+", "_", n)
  n <- gsub("([A-Z][a-z])", "_\\1", n)
  n <- tolower(trimws(n))
  
  n <- gsub("(^_+|_+$)", "", n)
  
  n <- gsub("_+", "_", n)
  
  if (unique) n <- make.unique(n, sep = "_")
  
  if (is.data.frame(.data)) {
    colnames(.data) <- n
    .data
  } else {
    n
  }
}

#### Custom `hist` function to plot "step" line histogram ####

hist_step <- function(x, breaks = 30, freq = TRUE, lwd = 1, lty = 1, col = 'black', add = FALSE, xlab = NULL, main = NA, xaxt = NULL, yaxt = NULL, ...) {
    hst <- hist(x, breaks = breaks, plot = FALSE)
    
    if (is.null(xlab)) {
        xlab <- deparse(substitute(x))
    }

    binwidth <- 2*(hst$breaks[1] - hst$mids[1])
    
    if (add) {
        if (freq) {
            lines(sort(c(hst$breaks[-1], hst$breaks[-1] + binwidth) ), rep(hst$counts, each = 2), 
                  col = col, 
                  lwd = lwd, 
                  lty = lty)
        } else if (!freq) {
            lines(sort(c(hst$breaks[-1],hst$breaks[-1] + binwidth) ), rep(hst$density, each = 2), 
                  col = col, 
                  lwd = lwd, 
                  lty = lty)
        }
    } else {
        if (freq) {
            plot(sort(c(hst$breaks[-1], hst$breaks[-1] + binwidth) ), rep(hst$counts, each = 2), 
                 ylab = 'Frequency', 
                 xlab = xlab, 
                 type = 'l', 
                 col = col, 
                 lwd = lwd, 
                 lty = lty,
		 main = main,
		xaxt = xaxt,
		yaxt = yaxt)
        } else if (!freq) {
            plot(sort(c(hst$breaks[-1], hst$breaks[-1] + binwidth) ), rep(hst$density, each = 2), 
                 ylab = 'Density', 
                 xlab = xlab, 
                 type = 'l', 
                 col = col, 
                 lwd = lwd, 
                 lty = lty,
		 main = main,
		xaxt = xaxt,
		yaxt = yaxt)
        }
    }
}

#' @description
#' `plotHomerKnown` plots known motif enrichment results from HOMER 
#' findMotifs tool and returns a formatted table which includes caclulated 
#' log2 fold change of target sequences over background
#' Percent columns are converted to numeric values for downstream processing
#' 
#' @details
#' Depends on `dplyr`, `ggplot2`, `RColorBrewer`, `forcats`, `grDevices`
#' and `data.table`
#' 
#' @path Path to knownResults.txt file
#' @ntop Plot the n top results
#' @by Filter the top results by this column
#' @sort Sort the plot labes by this column
#' @col (optional) A vector of colors
#' @title (optional) A plot title
#' 
#' @return prints a dotplot with results and quitely returns a list
#' with the plot and the results table

plotHomerKnown <- function(path, ntop=10, by='q_value_benjamini', sort='l2f_enrichment', col=NULL, title=NULL) {
    require(dplyr)
    require(ggplot2)
    require(RColorBrewer)
    require(forcats)
    require(grDevices)
    require(data.table)
    
    results <- fread(path, data.table = FALSE)
    results <- clean_names(results)
    colnames(results)[6] <- 'n_of_target_sequences_with_motif'
    colnames(results)[8] <- 'n_of_background_sequences_with_motif'
    
    results <- results %>%
    mutate(pct_of_target_sequences_with_motif = as.numeric(gsub('%','',pct_of_target_sequences_with_motif)) / 100,
          pct_of_background_sequences_with_motif = as.numeric(gsub('%','',pct_of_background_sequences_with_motif)) / 100) %>%
    mutate(l2f_enrichment = log2(pct_of_target_sequences_with_motif) - log2(pct_of_background_sequences_with_motif)) %>%
    mutate(label = gsub('\\-ChIP.*', '', motif_name)) %>%
    mutate(label = ifelse(nchar(label) > 50, paste0(substr(label, 1, 50),'...'), label)) %>%
    mutate(q_value_benjamini = ifelse(q_value_benjamini == 0, 1e-04, q_value_benjamini))
    
    if (by %in% colnames(results)) {
        if (by == 'q_value_benjamini') {
            to.plt <- results %>%
            slice_max(n = ntop, order_by = -log(q_value_benjamini), with_ties = FALSE) %>%
            mutate(label = fct_reorder(.f = as.factor(label), .x = get(sort)))
        
        } else {
                    
            to.plt <- results %>%
            slice_max(n = ntop, order_by = get(by), with_ties = FALSE) %>%
            mutate(label = fct_reorder(.f = as.factor(label), .x = get(sort)))
        }
    } else {
        stop(sprintf("%s is not a valid field to filter by. Chose one of: 'p_value', 'log_p_value', 'q_value_benjamini', 'l2f_enrichment'"))
    }

    if (is.null(col)) {
        plt <- to.plt %>%
        mutate(q_value_benjamini = ifelse(q_value_benjamini == 0, 1e-04, q_value_benjamini)) %>%
        ggplot(aes(x = l2f_enrichment, y = label, size = pct_of_target_sequences_with_motif, fill = -log10(q_value_benjamini))) +
        geom_point(shape = 21) +
        scale_fill_distiller(palette = 'PuBuGn', 
                             direction = 1, 
                             name = expression('-'*italic('log')[10]*' q-value'), 
                             limits = c(0,4), 
                             breaks = seq(0,4, length.out = 5), 
                             labels = c(0,(10^-seq(0,3,length.out = 4)[-1]), '≤1e-04')
                            ) +
        scale_size_continuous(name = 'Motif frequency') +
        scale_y_discrete(limits = to.plt %>% arrange(l2f_enrichment) %>% pull(label)) +
        theme_linedraw() +
        xlab(expression(italic(log)[2]*' fold enrichment over background')) +
        ylab('Motif name') +
        ggtitle(ifelse(is.null(title),'Top 10 motifs by q-value (Benjamini)', title)) +
        theme(plot.title = element_text(hjust = .5, face = 'bold'))
        
        } else {
            
        pal <- colorRampPalette(col)(30)
        plt <- to.plt %>%
        mutate(q_value_benjamini = ifelse(q_value_benjamini == 0, 1e-04, q_value_benjamini)) %>%
        ggplot(aes(x = l2f_enrichment, y = label, size = pct_of_target_sequences_with_motif, fill = -log10(q_value_benjamini))) +
        geom_point(shape = 21) +
        scale_fill_gradientn(colours = pal,
                             name = expression('-'*italic('log')[10]*' q-value'), 
                             limits = c(0,4), 
                             breaks = seq(0,4, length.out = 5), 
                             labels = c(0,(10^-seq(0,3,length.out = 4)[-1]), '≤1e-04')
                            ) +
        scale_size_continuous(name = 'Motif frequency') +
        scale_y_discrete(limits = to.plt %>% arrange(l2f_enrichment) %>% pull(label)) +
        theme_linedraw() +
        xlab(expression(italic(log)[2]*' fold enrichment over background')) +
        ylab('Motif name') +
        ggtitle(ifelse(is.null(title),sprintf('Top %d motifs by %s',ntop, by), title)) +
        theme(plot.title = element_text(hjust = .5, face = 'bold'))
    }
    print(plt)
    invisible(list('plot' = plt, 'results' = results))    
}


#' Save workspace progress to a file
#'
#' This function saves the current workspace to a file with a filename that
#' includes the current date and time, in the specified directory. If the
#' overwrite argument is set to FALSE and a progress file with the same name
#' already exists in the directory, the function will not save the workspace
#' again.
#'
#' @param dir A character string specifying the directory to save the progress file in. 
#' If not specified, defaults to the current working directory.
#' @param overwrite A logical value indicating whether to overwrite previous progress 
#' files that already exists in the directory.
#'
#' @return Prints a message indicating whether the image was saved successfully.
#'
#' @examples
#' # Save workspace progress in the current working directory
#' save.progress()
#'
#' # Save workspace progress in a specific directory
#' save.progress(dir = '~/my-project/progress/')
#'
#' # Save workspace progress but do not overwrite previous progress file
#' save.progress(overwrite = FALSE)
#'
#' @export
save.progress <- function(dir, overwrite = TRUE) {
    if (rlang::is_missing(dir)) {
        dir <- getwd()
    }
    date <- format(Sys.time(), format = "%Y%m%d_%H%M%S")
    savename <- paste0('save_progress_',date,'.Rdata')

    if (overwrite) {
        old.file <- list.files(dir, pattern = 'save_progress_')
        if (length(old.file) > 0) {
            message(sprintf('Previous files found: %s. Overwriting...', paste(old.file, collapse = ',')))
            file.remove(old.file)
            save.image(savename)
            message('Image saved successfully.')
        } else {
            message('No previous progress files found! Proceeding with saving a new.')
            save.image(savename)
            message('Image saved successfully.')
        }

    } else {
        save.image(savename)
        message('Image saved successfully.\n')
    }
}
