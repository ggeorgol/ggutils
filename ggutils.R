##### ============================================= #####
#####               Custom R functions              #####
#####    (c) Grigorios Georgolopoulos 2018-2024     #####
#####  		   georgolog@gmail.com		    #####
##### ============================================= #####

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
x_ordered <- x[order(x[,1],x[,2],x[,3]),]
write.table(x_ordered, file = file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

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

get_gene_names= function(keys, from_type, to_type, format = NA, fill.na = TRUE) {

library("org.Hs.eg.db")

conversion = mapIds(org.Hs.eg.db, keys = keys, keytype = from_type, column=to_type)

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
  n <- gsub("^=", "equals", n)
  n <- gsub("=", "_equals_", n)
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
