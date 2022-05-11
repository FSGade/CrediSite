#' Title
#'
#' @param consensus
#' @param file
#' @param phenotype_transform
#'
#' @return
#' @export
#'
#' @examples
msa_to_model_data <- function (file, consensus, phenotype_transform = function (x) { x }) {
  msa <- read_fasta(file)
  msa$phenotype <- sapply(
    get_values(msa$fasta_header),
    phenotype_transform
  )

  genotype = matrix(NA, nrow = dim(msa)[1], ncol = nchar(msa$sequence[1]))
  consensus_len = nchar(consensus)

  for (i in seq_along(msa$sequence)) {
    genotype[i, ] <- strsplit(msa$sequence[i], "")[[1]]
  }

  genotype <- data.frame(genotype,
                         stringsAsFactors = TRUE,
                         row.names = seq(dim(genotype)[1]))

  colnames(genotype) <- sapply(
    seq(dim(genotype)[2]),
    function (x) {
      paste("P", x, sep="")
    })

  for (i in seq_along(genotype)){
    genotype[,i] <- relevel(genotype[,i], substr(consensus, i, i))
  }
  genotype = genotype[,sapply(genotype, nlevels)>1]

  return(list(genotype = genotype, phenotype = msa$phenotype))
}

#' @rdname msa_to_model_data
#' @export
read_fasta <- function (file) {
  # From SigniSite package
  if (!is.character(file)) {
    stop("'file' has to be a string specifying a file name")
  }
  if (!file.exists(file)) {
    stop(paste("Unable to read file", file))
  }
  lines = readLines(con = file)
  comments = grep(pattern = "^#", x = lines)
  if (length(comments > 0)) {
    lines = lines[-comments]
  }
  empty_lines = grep(pattern = "^$", x = lines)
  if (length(empty_lines) > 0) {
    lines = lines[-empty_lines]
  }
  n_seqs = length(grep("^>", lines))
  headers = rep(NA, n_seqs)
  sequences = rep("", n_seqs)
  entry_no = 1
  for (line in lines) {
    if (grepl(pattern = "^>", x = line)) {
      headers[entry_no] = line
      entry_no = entry_no + 1
    }
    else {
      sequences[entry_no - 1] = paste0(sequences[entry_no -
                                                   1], line)
    }
  }
  sequences = toupper(x = sequences)
  aa_pattern = "[^ARNDCQEGHILKMFPSTWYVX-]"
  if (any(grepl(pattern = aa_pattern, x = sequences))) {
    msg = paste("Non standard character in sequences.", "Only 'ARNDCQEGHILKMFPSTWYVX-' allowed!",
                "Check your alignment!", "Proceeding by replacing with 'X'")
    warning(msg)
    sequences = gsub(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]",
                     replacement = "X", x = sequences)
  }
  #if (length(unique(nchar(ALIGNMENT$sequence))) > 1) {
  #  warning("SigniSite: The read FASTA files contain sequences of different length. To perform a SigniSite analysis, sequences must be pre-aligned")
  #}
  return(data.frame(fasta_header = headers, sequence = sequences,
                    stringsAsFactors = FALSE))
}

#' @rdname msa_to_model_data
#' @export
get_values <- function (fasta_header) {
  # From SigniSite package
  split_lst = strsplit(x = fasta_header, split = "\\s+")
  values = unlist(lapply(split_lst, function(x) {
    return(x[length(x)])
  }))
  return(as.numeric(values))
}
