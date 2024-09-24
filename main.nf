#!/usr/bin/env nextflow

// fetch SRA meta data using direct and esearch
// GEO accession number: GSE120011
// nextflow run main.nf --sra_id 'SRP161855' --geo_series 'GSE120011' -with-singularity -resume

params.geo_series = 'GSE120011'
params.sra_id = 'SRP161855'
params.outdir = "$projectDir/results"

log.info """\
  Fetching SRA mata data using esearch containier (docker)
  ================================================
  SRA accession number: ${params.sra_id}
  GEO accession number: ${params.geo_series}
"""

process GEO_SERIES_MATRIX {
  debug true
  container 'ncbi/edirect:latest'
  publishDir params.outdir, mode: 'copy'

  input:
  val geo_series

  output:
  path "${geo_series}_series_matrix.txt"

  script:
  """
  geo_url="https://ftp.ncbi.nlm.nih.gov/geo/series/"
  replace_nnn=\$(echo ${geo_series} | sed 's/...\$/nnn/')
  url=\${geo_url}/\${replace_nnn}/${geo_series}/matrix/${geo_series}_series_matrix.txt.gz
  # wget \$url 
  curl -o ${geo_series}_series_matrix.txt.gz \$url 
  gunzip ${geo_series}_series_matrix.txt.gz 
  """
}

process EDIRECT {
    debug true
    container 'ncbi/edirect:latest'
    publishDir params.outdir, mode: 'copy'


    input:
    val acc_number

    output:
    path "${acc_number}_run_info.csv"

    shell:
    '''
    esearch -db sra -query !{acc_number} | \
      efetch -format runinfo >  \
      !{acc_number}_run_info.csv
    '''

}

process JOIN_SPR_GEO {
  debug true
  tag "join SRA and GEO metadata"
  publishDir params.outdir, mode: 'copy'
  container = 'rocker/tidyverse:latest'

  input:
  path series_mat
  path run_info

  output:
  path 'metadata_sheet.csv'

  script:
  """
  #!/usr/bin/env Rscript

    library(dplyr)
    # GEO seris matrix
    geo_df <- readr::read_delim('${series_mat}', skip=37, 
                                col_names=FALSE, n_max=11, delim='\t') %>%
      dplyr::filter(X1 %in% c('!Sample_title', '!Sample_geo_accession', 
                              '!Sample_type', '!Sample_source_name_ch1',
                              '!Sample_organism_ch1',
                              '!Sample_characteristics_ch1' )) %>%
      dplyr::mutate(X1 = stringr::str_replace(X1, '!', '')) %>%
      dplyr::mutate(X1 = stringr::str_replace(X1, '_ch1', ''))                                                                  

    geo_sample_df <- t(geo_df[, -1]) %>% as.data.frame()  
    names(geo_sample_df)[1:5] <- geo_df\$X1[1:5]

    # SRA run info
    df <- readr::read_delim('${run_info}') %>%
      dplyr::group_by(SampleName) %>%
      dplyr::mutate(replicate = row_number()) %>%
      dplyr::rename(Sample_geo_accession='SampleName') %>%
      dplyr::left_join(geo_sample_df, by='Sample_geo_accession')

    write.csv(df, file='metadata_sheet.csv', row.names = FALSE)

  """
}

workflow {
    /*  fecth geo series matrix txt file */
    geo_ch = channel.of(params.geo_series)
    series_matrix = GEO_SERIES_MATRIX(geo_ch)
    series_matrix.view()

    /* use edirect to get run info */
    sra_ch = channel.of(params.sra_id)
    run_info = EDIRECT(sra_ch)
    run_info.view()
    /*
    results = JOIN_SPR_GEO(series_matrix, run_info)
    results.view()
    */

}

