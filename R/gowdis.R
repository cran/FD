`gowdis` <-
function(x, w, asym.bin = NULL, ord = c("podani", "metric") ){

if (length(dx <- dim(x)) != 2 || !(is.data.frame(x) || is.numeric(x))) stop("x is not a dataframe or a numeric matrix.","\n")

# n = number of rows, p = number of variables
n <- dx[1]
p <- dx[2]

ord <- match.arg(ord)

varnames <- dimnames(x)[[2]]

# check for weight vector, add equal weights if missing
if (!missing(w)){

	# check if correct class and length
	if (length(w) != p | !is.numeric(w) ) stop("w needs to be a numeric vector of length = number of variables in x.","\n")
	w <- w / sum(w)
	}
else w <- rep(1, p) / sum(rep(1, p))


if (is.data.frame(x)) {
        type <- sapply(x, data.class)
        }
    else {
        type <- rep("numeric", p)
        names(type) <- colnames(x)
	}

# replace character variables by factors
for (i in 1:p) if (type[i] == "character") x[,i] <- as.factor(x[,i]) else x[,i] <- x[,i]

# check for binary variables
is.bin <- function(k) all(k[!is.na(k)] %in% c(0,1))

bin.var <- rep(NA,p); names(bin.var) <- varnames
for (i in 1:p) bin.var[i] <- is.bin(x[,i])

if (any(type[bin.var] != "numeric")) stop("Binary variables should be of class 'numeric'.")

type[type %in% c("numeric", "integer")] <- "C"
type[type == "ordered"] <- "O"
type[type %in% c("factor", "character")] <- "N"
type[bin.var] <- "B"

# convert asymmetric binary variables, if present
if (!is.null(asym.bin) ){
	if (!all(bin.var[asym.bin])) stop("Asymetric binary variables must only contain 0 or 1.")
	else type[asym.bin] <- "A"
	}

# convert factors to their internal numeric codes
x <- data.matrix(x)

# convert ordinal variables to ranks, following Podani (1999)
for (i in 1:p) if (type[i] == "O") x[,i] <- rank(x[,i], na.last = "keep") else x[,i] <- x[,i]


# create a vector of the coordinates (position) of each value in the distance matrix
full.col <- rep(1:n, n)
full.row <- rep(1:n, rep(n, n))
full.pos <- cbind(full.col, full.row)


# compute the range of each variable (this will only be used for numeric and ordinal variables)
range.Data <- function(v){
	r.Data <- range(v, na.rm = T)
	res <- r.Data[2]-r.Data[1]
	return(res)
	}

range2<- apply(x, 2, range.Data)

# compute Timax, and Timin for each variable (these will only apply to ordinal variables, see Podani [1999], eq. 2b.)
comp.Timax <- function(v){
	Ti.max <- max(v, na.rm = T)
	no.na <- v[!is.na(v)]
	res <- length(no.na[no.na == Ti.max])
	return(res)
	}

Timax <- apply(x, 2, comp.Timax)
	
comp.Timin <- function(v){
	Ti.min <- min(v, na.rm = T)
	no.na <- v[!is.na(v)]
	res <- length(no.na[no.na == Ti.min])
	return(res)
	}

Timin <- apply(x, 2, comp.Timin)


# function to compute partial similarities between two rows of x
part.sim <- function(pos, Data, type, weight, ord, range2, Timax, Timin){

	j <- pos[1]
	k <- pos[2]

	# compute Tij, and Tik for each variable (even though these will only apply to ordinal variables, see Podani [1999], eq. 2b.)
	comp.Tij <- function(v, pos){
			j <- pos[1]
			Tj <- v[j]
			if (is.na(Tj)) res <- NA
			else{
			    no.na <- v[!is.na(v)]
			    temp1 <- no.na[no.na == Tj]
		 	    res <- length(temp1)
			  }
			return(res)	
			}
	Tij <- apply(Data, 2, comp.Tij, pos=pos)	

	comp.Tik <- function(v, pos){
			k <- pos[2]
			Tk <- v[k]
			if (is.na(Tk)) res <- NA
			else{
			    no.na <- v[!is.na(v)]
			    temp2 <- no.na[no.na == Tk]
		 	    res <- length(temp2)
			  }
			return(res)	
			}
	Tik <- apply(Data, 2, comp.Tik, pos=pos)

	temp <- rbind(Data[j,], Data[k,], weight, range2, Tij, Tik, Timax, Timin)
	
	# put weights to 0 for j-k comparisons with NAs
	put.na <- function(v){
		if (is.na(v[1]) || is.na(v[2])) v[3] <- 0
		else v[3] <- v[3]
		return(v)
		}

	temp <- apply(temp, 2, put.na)
	temp <- t(temp)
	temp <- data.frame(temp, type)

	# put weights to 0 for j-k comparisons involving double-zeros for asymmetric binary variables
	temp[temp$type == "A" & !is.na(temp[1]) & !is.na(temp[2]) & temp[1] == 0 & temp[2] == 0,3] <- 0

	# convert type to numbers
	temp[,9] <- as.character(temp[,9])
	temp[temp$type == "C",9] <- 1
	temp[temp$type == "N",9] <- 2
	temp[temp$type == "A",9] <- 3
	temp[temp$type == "O",9] <- 4
	temp[temp$type == "B",9] <- 5
	temp[,9] <- as.numeric(temp[,9])	

	part.sim2 <- function(v, ord){
		
		# if one value is NA, sim = 0 (note that weight = 0)
		if (is.na(v[1]) | is.na(v[2])) sim <- 0
		
		else{
			# for type C (1)
			if (v[9] == 1) sim <- 1- ( abs(v[1] - v[2]) / v[4] )
			
			# note that the weights for double zeros for type A have already been set to 0, so we can treat them like nominal variables here
		
			# for types N (2) or A (3) or B (5)
			if (v[9] == 2 | v[9] == 3 | v[9] == 5) sim <- ifelse(v[1] == v[2], 1 ,0)
		
			# for type O (4)
			# The following is Eqs. 2a-b of Podani (1999)
			if (v[9] == 4 & ord == "podani"){
				
				if (v[1] != v[2]) sim <- 1 - ( (abs(v[1] - v[2]) - (v[5] - 1)/2 - (v[6] - 1)/2 ) / (v[4] - (v[7] - 1)/2 - (v[8] - 1)/2) )
				else sim <- 1
				}
			# The following is Eq. 3 of Podani (1999)
			if (v[9] == 4 & ord == "metric") sim <- 1 - ( abs(v[1] - v[2]) / v[4] )
			}
		return(sim)
		} # end of part.sim2


	p.sim <- apply(temp, 1, part.sim2, ord=ord)
	dis.jk <- 1- (sum(temp[,"weight"] * p.sim) / sum(temp[,"weight"]) )
	return(dis.jk)
	} # end of part.sim function

# compute dissimilarities to lower half of dissimilarity matrix
full.pos <- full.pos[full.pos[,1] < full.pos[,2],]
full.pos<-full.pos[order(full.pos[,1],full.pos[,2]),]
full.dis <- apply(full.pos, 1, part.sim, Data=x, type=type, weight=w, ord=ord, range2=range2, Timax=Timax, Timin=Timin)

# create the dissimilarity matrix

dummy.dis <- matrix(NA, n, n)
dummy.dis[lower.tri(dummy.dis)]<-full.dis
full.dis <- as.dist(dummy.dis)

if (any(is.na(full.dis))) attr(full.dis, "NA.message") <- "NA-values in the dissimilarity matrix !"
attr(full.dis, "Labels") <- dimnames(x)[[1]]
attr(full.dis, "Size") <- n
attr(full.dis, "Metric") <- "Gower"
attr(full.dis, "Types") <- type

return(full.dis)

}

