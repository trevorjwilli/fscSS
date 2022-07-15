

#' Read in arp files of single DNA Loci
#'
#' This function reads in .arp files from fsc26 and creates DNAbin
#' objects for each simulated population
#'
#' @param file PATH to .arp file
#'
#' @export

read.arp <- function(file) {
  arp <- scan(file, what = "raw", sep = "\n", strip.white = T, quiet = T) # Read in file, file becomes a vector where each line is an entry as a string
  Nsamp <- as.numeric(gsub("NbSamples=", "", arp[grep("NbSamples", arp)])) # Calculate the number of samples by finding the line in the .arp file and getting rid of everything but the number
  #print(Nsamp)
  Nsite <- as.numeric(gsub(".+polymorphic sites: ", "", arp[grep("polymorphic sites", arp)])) # Same as above but to find the number of polymorphic sites
  #print(paste("Nsite = ", Nsite))
  if(Nsamp == 1) {
    linenum <- grep(paste0("SampleName=\"Sample ", 1), arp) # Find the line number to start on for the sample
    N <- as.numeric(gsub("SampleSize=", "", arp[linenum+1])) # Find out how many sequences are in the sample
    out <- matrix(nrow = N, ncol = Nsite)

    for(j in 1:N) { # Loop through the number of sequences in the sample
      index <- grep(paste0(1,"_",j,"\\t"), arp) # Find the line with the desired sequence
      #print(paste("Index = ", index))
      seq <- gsub(paste0(".+\\t"), "", arp[index]) # Extract the sequence as a single character string
      #print(seq)
      seq <- unlist(strsplit(seq, "")) # Split the sequence into a vector where each element is a nucleotide
      #print(length(seq))
      out[j,] <- seq # Add that sequence to the matrix
    }
    out <- ape::as.DNAbin(out)
    return(out)

  } else {

    out.list <- list() # Initialize listfor output

    for(i in 1:Nsamp) { # Loop through the number of samples
      linenum <- grep(paste0("SampleName=\"Sample ", i), arp) # Find the line number to start on for the sample
      N <- as.numeric(gsub("SampleSize=", "", arp[linenum+1])) # Find out how many sequences are in the sample
      out.list[[i]] <- matrix(nrow = N, ncol = Nsite) # Initialize matrix for sample sequences
      #print(out.list)

      for(j in 1:N) { # Loop through the number of sequences in the sample
        index <- grep(paste0(i,"_",j,"\\t"), arp) # Find the line with the desired sequence
        #print(paste("Index = ", index))
        seq <- gsub(paste0(".+\\t"), "", arp[index]) # Extract the sequence as a single character string
        #print(seq)
        seq <- unlist(strsplit(seq, "")) # Split the sequence into a vector where each element is a nucleotide
        #print(length(seq))
        out.list[[i]][j,] <- seq # Add that sequence to the matrix
      }
      rownames(out.list[[i]]) <- paste0(i,"-",1:N) # Add sequence names (following that of the fastsimcoal2 output)
    }
    samples <- paste0(rep("sample",Nsamp),".",1:Nsamp) # Make vector of names for list
    #nam <- c(samples, "full") # Same as above
    #out.list[[Nsamp+1]] <- do.call(rbind, out.list) # Add another matrix including all sequences
    #names(out.list) <- nam # Give the matrix objects names
    #out <- lapply(out.list, as.DNAbin) # Convert the matrices into DNABIN objects
    names(out.list) <- samples
    return(out.list)
    #names(out) <- samples
    #return(out)
  }
}

#' Read in a .param file
#'
#' @param path PATH to .param file
#'
#' @export

read.param <- function(path) {
  input <- readLines(path)
  for(i in 1:length(input)) {
    if(i == 1) {
      header <- unlist(strsplit(input[i], "\t"))
      print(header)
      rmv <- grep("Max", header)
      header <- header[-rmv]
      out <- matrix(nrow = length(input) - 1, ncol = length(header))
    } else {
      inline <- as.numeric(unlist(strsplit(input[i], "\t")))
      out[i - 1, ] <- inline
    }
  }
  colnames(out) <- header
  out <- as.data.frame(out)
  out
}



#' Concatenate sequences into a vector from matrix
#'
#' @param x A matrix of a multiple sequence alignment where each column is a base
#' pair and each row is a sequence
#'
#' @export

concat_seqs <- function(x) {
  y <- apply(x, 1, function(x) {paste(x, collapse = "")})
  return(y)
}

#' Function to calculate the number of shared haplotypes between populations
#'
#' @param x A list of matrices of multiple sequence alignments or DNAbins
#'
#' @export

sharedhaps <- function(x) {
  t <- lapply(x, fscSS::concat_seqs) # Turn matrices into vectors where each element is a single sequence and put all in a list
  combs <- utils::combn(t, 2, simplify = F) # Get all the pairwise combinations between populations
  all.inter <- lapply(combs, function(x) intersect(x[[1]], x[[2]])) # Find the intersection between all pairwise combinations
  out.un <- unlist(all.inter) # Make all shared haplotypes into a single vector
  shared.haps <- intersect(out.un, out.un) # Find which haplotypes are shared
  #print(paste0("Number of Shared Haplotypes = ", length(shared.haps)), quote = F)
  seqs <- unlist(t) # Put all sequences into one vector
  counts <- table(seqs) # Count how many times each haplotype occurs
  shared <- match(shared.haps, names(counts)) # Get counts of only shared haplotypes
  freqs <- as.vector(counts[shared]) # Same as above but convert to normal vector
  freqs <- freqs/length(seqs) # Calculate the percentage of all sequences the shared haplotypes exhibit
  #print(paste0("Percent Sequences Shared = ", sum(freqs)), quote = F)
  return(list(length(shared.haps), sum(freqs)))
}


#' Calculate the average number of nucleotide differences (Tajima 1983 eq. A3)
#'
#' @param x An object of class DNAbin or a matrix where rows are sequences and columns are nucleotide positions (assumed to be already aligned)
#'
#' @export

calc_pi <- function(x) {

  differ <- function(x, y) {
    n <- which(x != y)
    return(length(n))
  }

  l <- split(x, row(x))
  d <- sapply(l, function(m) sapply(l, function(n) differ(m, n)))
  d[upper.tri(d)] <- 0
  if(length(l) == 0) {
    ave.pi <- 0
  } else {
    ave.pi <- sum(d/choose(length(l), 2))
  }
  return(ave.pi)
}


#' Function to calculate the average number of nucleotide differences between two separate populations
#'
#' @param x,y Objects of class DNAbin or matrices where rows are sequences and columns are nucleotide positions from two different populations
#'
#' @export

calc_pi_pops <- function(x, y) {

  differ <- function(x, y) {
    n <- which(x != y)
    return(length(n))
  }

  l.1 <- split(x, row(x))
  l.2 <- split(y, row(y))
  d <- sapply(l.1, function(m) sapply(l.2, function(n) differ(m, n)))
  ave.pi <- sum(d)/length(d)
  return(ave.pi)
}

#' Function to calculate the Watterson estimator for theta
#'
#' @param x Object of class DNAbin or a matrix of multiple sequence alignment where
#' rows are sequences and columns are nucleotide positions
#'
#' @export


theta_watt <- function(x) {
  sites <- apply(x, 2, unique) # Find which sites have multiple nucleotides present
  if(is.matrix(sites)) {
    s <- ncol(sites)
  } else {
    s <- length(which(sapply(sites, length) != 1)) # Calculate the total number of segregating sites
  }
  theta <- s/sum(1/1:(nrow(x)-1))
  theta
}

#' Function to calculate Tajima's D
#'
#' @param x Object of class DNAbin or a matrix of multiple sequence alignment where
#' rows are sequences and columns are nucleotide positions
#'
#' @export

tajimaD <- function(x) {
  n <- nrow(x)
  pi <- fscSS::calc_pi(x)
  sites <- apply(x, 2, unique) # Find which sites have multiple nucleotides present
  if(is.matrix(sites)) {
    s <- ncol(sites)
  } else {
    s <- length(which(sapply(sites, length) != 1)) # Calculate the total number of segregating sites
  }
  tmp <- 1:(n - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (n + 1)/(3*(n - 1))
  b2 <- 2*(n^2 + n + 3)/(9*n*(n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1*n) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  D <- (pi - s/a1)/sqrt(e1*s + e2*s*(s-1))
  D
}

#' Function to calculate the number of haplotypes, number of segregating sites, haplotype diversity, and pi of a matrix of sequences
#'
#' @param x Object of class DNAbin or a matrix of multiple sequence alignment where
#' rows are sequences and columns are nucleotide positions
#'
#' @export

sumstats <- function(x) {
  n <- nrow(x)
  if(length(x) == 0) {
    num.haps <- 1
    hap.freqs <- 1
  } else {
    haps <- pegas::haplotype(x) # Find haplotypes using Pegas's haplotype function
    num.haps <- length(attr(haps, "index")) # Find the total number of haplotypes
    hap.freqs <- lengths(attr(haps, "index"))/n # Make a vector containing the frequencies of each haplotype
  }
  h <- (n*(1-sum(hap.freqs^2)))/(n-1) # Calculate haplotype diversity (Nei 1983, eq. 8.4)
  sites <- apply(x, 2, unique) # Find which sites have multiple nucleotides present
  if(is.matrix(sites)) { # Fixes an error where the output of sites is a matrix if all segregating sites have the same number of mutations
    s <- ncol(sites)
  } else if(is.list(sites)) {
    s <- length(which(sapply(sites, length) != 1)) # Calculate the total number of segregating sites
  } else if(num.haps == 1) {
    s <- 0
  }
  if(length(x[1,]) == 0) { # Add if statement to check if sequences are different so Pi can be calculated
    pi <- 0
  } else {
    pi <- fscSS::calc_pi(x) # Calculate the Average number of nucleotide differences (π)
  }
  theta <- fscSS::theta_watt(x)
  D <- fscSS::tajimaD(x)
  out <- c(num.haps, h, s, pi, theta, D)
  names(out) <- c("K", "H", "S", "Pi", "ThetaW", "TajimaD")
  return(out)
}


#' Function to calculate sumstats on samples with multiple populations
#'
#' @param x List of objects of class DNAbin or matrices of multiple sequence alignment where
#' rows are sequences and columns are nucleotide positions. Each object in list is a separate
#' population.
#'
#' @export


sumstats_pops <- function(x) {
  npops <- length(x)
  y <- lapply(x, fscSS::sumstats) # Calculate sumstats within populations
  y <- unlist(y) # make a vector of sumstats from populations
  ind.K <- grep("K", names(y)) # Find which vector elements correspond to pi
  mean.K <- mean(y[ind.K]) # Calculate the mean pi across populations
  sd.K <- stats::sd(y[ind.K]) # Calculate the standard deviation of pi across populations
  ind.S <- grep("S", names(y)) # Find which vector elements correspond to pi
  mean.S <- mean(y[ind.S]) # Calculate the mean pi across populations
  sd.S <- stats::sd(y[ind.S]) # Calculate the standard deviation of pi across populations
  ind <- grep("Pi", names(y)) # Find which vector elements correspond to pi
  mean.pi <- mean(y[ind]) # Calculate the mean pi across populations
  sd.pi <- stats::sd(y[ind]) # Calculate the standard deviation of pi across populations
  ind.h <- grep("H", names(y))
  mean.H <- mean(y[ind.h])
  sd.H <- stats::sd(y[ind.h])
  ind.theta <- grep("ThetaW", names(y))
  mean.theta <- mean(y[ind.theta])
  sd.theta <- stats::sd(y[ind.theta])
  ind.D <- grep("TajimaD", names(y))
  mean.D <- mean(y[ind.D])
  var.D <- stats::var(y[ind.D])
  all.samps <- do.call(rbind, x) # Make a new matrix pooling the populations together
  n <- nrow(all.samps)
  all.stats <- fscSS::sumstats(all.samps)
  shared.haps <- fscSS::sharedhaps(x)
  tot <- c(all.stats["K"], shared.haps[[1]], shared.haps[[2]], all.stats["S"], all.stats["H"], all.stats["Pi"], all.stats["ThetaW"], all.stats["TajimaD"],mean.K, sd.K, mean.S, sd.S, mean.H, sd.H, mean.pi, sd.pi, mean.theta, sd.theta, mean.D, var.D) # Output pooled summary statistics
  out <- append(y, tot)
  popshead <- vector()
  for(i in 1:npops) {
    popshead <- append(popshead, paste0(c("K", "H", "S", "Pi", "ThetaW", "TajimaD"), i))
  }
  names(out) <- c(popshead, "Tot.K","Num.shared", "Per.shared", "Tot.S", "Tot.H", "Tot.pi", "Tot.theta", "Tot.D","Mean.K", "sd.K", "Mean.S", "sd.S", "Mean.H", "sd.H", "Mean.Pi", "sd.Pi", "Mean.theta", "sd.theta", "Mean.D", "var.D")

  return(out)
}

#' Function to calculate simple transition/transversion ratio and transition rate
#'
#' @param x Object of class DNAbin
#'
#' @export

count_tstv <- function(x) {
  ts <- 0
  tv <- 0

  for(i in 1:ncol(x)) {
    alleles <- unique(as.vector(as.character(x[,i])))
    if(length(alleles) == 1){
      next
    } else if(length(alleles) == 2) {
      if("a" %in% alleles & "g" %in% alleles) {
        ts <- ts + 1
      } else if("c" %in% alleles & "t" %in% alleles) {
        ts <- ts + 1
      } else {
        tv <- tv + 1
      }
    } else if(length(alleles) > 2) {
      warning(paste0("Multiple hits at position ", i, ", position skipped"))
      next
    }
  }
  cat(paste("Transitions:", ts, "\n"))
  cat(paste("Transversions:", tv, "\n"))
  cat(paste("ts/tv:", ts/tv, "\n"))
  cat(paste("ts rate:", ts/sum(ts, tv), "\n"))
}




