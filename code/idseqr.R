read_reports <- function(reports_dir,
                         tax_level=1,
                         min_nt_alignment_length=70,
                         min_nr_alignment_length=0) {
  taxlevel <- tax_level #avoid collisions with dataframe column
  list.files(reports_dir) %>%
    lapply(function(x){
      sample_name <- stringr::str_match(x, "^(.*)_\\d+_taxon_report.csv$")[,2]
      fname <- paste(reports_dir, "/", x, sep="")
      if (file.size(fname) == 0) return(NULL)
      cbind(sample_name=sample_name,
            read.csv(fname, stringsAsFactors=F),
            stringsAsFactors=F)
    }) %>%
    do.call(what=dplyr::bind_rows) %>%
    dplyr::filter(tax_id > 0) %>%
    dplyr::filter(tax_level == taxlevel) %>%
    dplyr::filter(nt_alignment_length >= min_nt_alignment_length,
                  nr_alignment_length >= min_nr_alignment_length)
}

filter_background <- function(reports,
                              # sample names (character)
                              controls,
                              # names should correspond to sample names
                              batches=NULL,
                              counts_column="nt_count",
                              regularize_bg=TRUE,
                              # if NULL, normalizes by summing up reads in reports
                              # recommend to use total_ercc_reads when appropriate
                              normalization=NULL) {
  samples <- unique(reports$sample_name)
  taxa <- unique(reports$name)

  counts_mat <- Matrix::sparseMatrix(
    i=as.integer(factor(reports$sample_name, levels=samples)),
    j=as.integer(factor(reports$name, levels=taxa)),
    x=reports[,counts_column],
    dimnames=list(samples, taxa))

  if (is.null(batches)) {
    batches <- rep("batch0", length(samples))
    names(batches) <- samples
  }

  if(is.null(normalization)) {
    normalization <- Matrix::rowSums(counts_mat)
  }

  normalized_counts_mat <- counts_mat / normalization
  sapply(unique(batches[controls]), function(b) {
    controls %>%
      .[batches[.] == b] ->
      batch_controls
    normalized_counts_mat[batch_controls,,drop=FALSE] %>%
      Matrix::colMeans()
  }) %>%
    t() ->
    batch_bg

  global_bg <- colMeans(batch_bg)

  sapply(unique(batches), function(b) {
    controls %>%
      .[batches[.] == b] ->
      batch_controls
    mat <- normalized_counts_mat[batch_controls,,drop=FALSE]
    if (regularize_bg) {
      mat <- rbind(mat, global_bg)
    }
    Matrix::colMeans(mat)
  }) %>%
    t() ->
    batch_bg

  batch_bg %>%
    .[,colSums(.) > 0,drop=FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("batch") %>%
    tidyr::gather(name, background, -batch) ->
    batch_bg_melted

  data.frame(
    sample_name=reports$sample_name,
    name=reports$name,
    count=reports[,counts_column],
    stringsAsFactors=F) ->
    counts_df

  dplyr::inner_join(
    data.frame(sample_name=controls,
               batch=batches[controls],
               norm=normalization[controls],
               stringsAsFactors=F),
    batch_bg_melted
  ) %>%
  dplyr::left_join(
    counts_df %>%
      .[.$sample_name %in% controls,]
  ) %>%
    tidyr::replace_na(list(count=0)) %>%
    dplyr::mutate(expected_bg_count=normalization[.$sample_name] * background) %>%
    dplyr::filter(expected_bg_count > 0) ->
    model_df

  model_df %>%
    {MASS::glm.nb(count ~ offset(expected_bg_count) - 1, data=., link=identity)} ->
    bg_glm

  theta <- MASS::theta.md(model_df$count, fitted(bg_glm), dfr=df.residual(bg_glm))

  counts_df %>%
    dplyr::mutate(background=batch_bg[as.matrix(cbind(batches[sample_name], name))]) %>%
    dplyr::mutate(expected_bg_count=normalization[sample_name] * background) %>%
    dplyr::mutate(p_val=pnbinom(count, size=theta,
                                mu=expected_bg_count, lower.tail=FALSE)) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::ungroup() ->
    counts_df

  stopifnot(counts_df$sample_name == reports$sample_name)
  stopifnot(counts_df$name == reports$name)

  ret <- cbind(reports,
               counts_df[,c("expected_bg_count", "p_val")],
               stringsAsFactors=F)
  attr(ret, "bg_glm") <- bg_glm
  ret
  }

# filter to keep only the top n taxa per sample, for plotting
filter_top_taxa <- function(reports, top_tax_per_sample=5, counts_column="nt_count") {
  cbind(
    reports[,c("sample_name", "name")],
    count=reports[,counts_column],
    stringsAsFactors=F
  ) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::do(head(., n=top_tax_per_sample)) %>%
    dplyr::ungroup() %>%
    .$name %>% unique() ->
    keep

  reports %>%
    dplyr::filter(name %in% keep)
}

# TODO:
# idseq2phyloseq()
