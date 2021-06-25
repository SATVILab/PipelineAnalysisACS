#' @title Pre-process data
#' TODO: Sum over 040238 samples
preprocess_fn <- function(data_raw, params_dots, dir_proj){

  if('skip_if_html_found' %in% names(params_dots)){
    if(params_dots$skip_if_html_found){
      if(file.exists(file.path(dir_proj, "output.html"))) stop('Stopping run because html already found.')
    }
  }

  print(dir_proj)
  data_raw %<>% tibble::as_tibble()

  # ===========================
  # Filter
  # ===========================

  if(length(params_dots$filter) != 0){
    for(i in seq_along(params_dots$filter)){
      var <- names(params_dots$filter)[i]
      val <- params_dots$filter[[i]]
      data_raw <- data_raw[data_raw[[var]] %in% val,]
    }
  }

  if('filter_non_responders' %in% names(params_dots)){
    if(params_dots$filter_non_responders){

      scores_mat <- readRDS(file.path(usethis::proj_path(), 'analysis', "supporting",
                                      "cytof", paste0(params_dots$stim, "_scores.rds")))
      sampleid_vec <- rownames(scores_mat)
      scores_tbl <- scores_mat %>%
        tibble::as_tibble() %>%
        dplyr::mutate(SampleID = sampleid_vec) %>%
        tidyr::separate(SampleID, c("SubjectID", "VisitType")) %>%
        dplyr::select(SubjectID, VisitType, everything())


      #plot_tbl <- scores_tbl %>%
      #  pivot_longer(-c(SubjectID, VisitType),
      #               names_to = "cyt_combn",
      #               values_to = "score")

      #windows()
      #ggplot(plot_tbl, aes(x = cyt_combn, y = score)) +
      #  geom_boxplot() +
      #  theme(axis.text.x = element_text(angle = 90))

      #params_dots <- list(cyt_combn_min_prob_vec = cyt_combn_min_prob_vec)
      min_prob_ind_list <- map(names(params_dots$cyt_combn_min_prob_vec), function(cyt_combn){
        min_prob <- params_dots$cyt_combn_min_prob_vec[[cyt_combn]]
        which(scores_tbl[[cyt_combn]]  >= min_prob)
      })

      if(params_dots$min_prob_type == "any"){
        min_prob_ind_vec <- Reduce(union, min_prob_ind_list)
      } else if(params_dots$min_prob_type == "all"){
        min_prob_ind_vec <- Reduce(intersect, min_prob_ind_list)
      } else{
        stop("params_dots$min_prob_type is not either 'all' or 'any'")
      }

      scores_filter_vec <- scores_tbl[min_prob_ind_vec,] %>%
        dplyr::mutate(SampleID = paste0(SubjectID, "_", VisitType)) %>%
        magrittr::extract2('SampleID')

      data_raw %<>%
        dplyr::mutate(SampleID = paste0(SubjectID, "_", VisitType)) %>%
        dplyr::filter(SampleID %in% scores_filter_vec) %>%
        dplyr::select(-SampleID)

    }
  }


  # ===========================
  # Transform
  # ===========================

  if('trans' %in% names(params_dots)){
    if(params_dots$trans == 'sqrt_pos'){
      data_raw %<>%
        dplyr::mutate(resp = pmax(0, resp),
               resp = sqrt(resp))
    } else if(params_dots$trans == 'pos'){
      data_raw %<>%
        dplyr::mutate(resp = pmax(0, resp))
    }
  }

  # ===========================
  # Add clinical data info
  # ===========================

  if(params_dots$response_type == 'hladr_med_diff'){
    if('timeToTBFromVisit' %in% colnames(data_raw)){
      data_raw <- data_raw[,-which(colnames(data_raw) == 'timeToTBFromVisit')]
    }
    if('Progressor' %in% colnames(data_raw)){
      data_raw <- data_raw[,-which(colnames(data_raw) == 'Progressor')]
    }
  }

  if(params_dots$assay %in% c('cytof', 'cytof-faust', 'cytof_ag_flowsom')){
    if(params_dots$response_type %in% c('count', 'prob')){
      data_raw %<>%
        dplyr::left_join(clinical_data %>%
                    dplyr::select(SubjectID, Progressor, timeToTBFromVisit, VisitType),
                  by = c("SubjectID", "VisitType"))
    } else if(params_dots$response_type %in% 'compass'){
      data_raw %<>%
        dplyr::left_join(clinical_data %>%
                    dplyr::select(SubjectID, timeToTBFromVisit, VisitType) ,
                  by = c("SubjectID", "VisitType")) %>%
        dplyr::select(-timeToTB)
    } else if(params_dots$response_type == 'hladr_med_diff'){
      data_raw %<>%
        dplyr::left_join(clinical_data %>%
                    dplyr::select(SubjectID, Progressor, timeToTBFromVisit, VisitType) %>%
                    dplyr::group_by(SubjectID, VisitType) %>%
                    dplyr::slice(1) %>%
                    dplyr::ungroup(),
                  by = c("SubjectID", "VisitType"))
    }

  } else if(params_dots$assay == 'rnaseq' || params_dots$assay == 'antibody' || params_dots$response_type == 'hladr_med_diff'){
    data_raw %<>%
      dplyr::left_join(clinical_data %>%
                  dplyr::select(SubjectID, Progressor, timeToTBFromVisit, VisitType) %>%
                  dplyr::group_by(SubjectID, VisitType) %>%
                  dplyr::slice(1) %>%
                  dplyr::ungroup(),
                by = c("SubjectID", "VisitType"))
  }

  # remove samples with time to tb from visit less than 0
  data_raw %<>%
    dplyr::filter(purrr::map_lgl(data_raw$timeToTBFromVisit, function(x){
      if(is.na(x)) return(TRUE)
      x > 0
    }))

  # choose only progressor or non-prog samples
  data_raw %<>% dplyr::filter(Progressor %in% c('yes', 'no'))

  # add other covariates
  data_raw %<>%
    dplyr::left_join(params_dots$clinical_data %>%
                dplyr::group_by(SubjectID, VisitType) %>%
                dplyr::slice(1) %>%
                dplyr::ungroup() %>%
                dplyr::select(SubjectID, VisitType, Sex, Height, Weight, AgeAtLastBirthDay,
                       Ethnicity, DaysSinceEntry, VisitDate, PreviousDiagnosisOfTB),
              by = c("SubjectID", "VisitType")) %>%
    dplyr::mutate(BMI = Weight/Height^2) %>%
    dplyr::mutate(Year = lubridate::year(VisitDate),
           WinterDate = lubridate::ymd(paste0(Year, "-06-21")),
           DaysFromWinter = VisitDate - WinterDate,
           DaysFromWinter = abs(as.numeric(DaysFromWinter)),
           DaysFromWinter = ifelse(DaysFromWinter < 45, 0, DaysFromWinter),
           Winter = ifelse(DaysFromWinter == 0, "yes", "no")) %>%
    dplyr::select(-c(Year, WinterDate)) %>%
    dplyr::select(dataset:VisitType, Sex:DaysSinceEntry, BMI, DaysFromWinter, Winter,
           PreviousDiagnosisOfTB, everything())

  if(any(purrr::map_lgl(params_dots$conf, function(x) stringr::str_detect(x, 'sig6gene_CorScore')))){
    data_raw %<>%
      dplyr::left_join(TuberculomicsCompendium::signature_6gene,
                by = c("SubjectID", "VisitType")) %>%
      dplyr::select(dataset:Winter, sig6gene_CorScore, PreviousDiagnosisOfTB,
             everything())
    data_raw %<>%
      dplyr::filter(!is.na(sig6gene_CorScore))
  }

  if(names(params_dots$var_exp_spline) != 'tfmttb'){
    if(names(params_dots$var_exp_spline) %in% TuberculomicsCompendium::soma_data_tidy$Soma_Target){
      var_tbl_add <- TuberculomicsCompendium::soma_data_tidy %>%
        dplyr::filter(Soma_Target == names(params_dots$var_exp_spline)) %>%
        dplyr::group_by(SubjectID, VisitType) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(SubjectID, VisitType, Soma_TransformedReadout)
      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(params_dots$var_exp_spline))
      data_raw %<>%
        inner_join(var_tbl_add ,
                  by = c("SubjectID", "VisitType"))
    } else if(names(params_dots$var_exp_spline) == 'risk6'){
      var_tbl_add <- TuberculomicsCompendium::signature_6gene %>%
        dplyr::select(SubjectID, VisitType, sig6gene_CorScore)
      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(params_dots$var_exp_spline))
      data_raw %<>%
        inner_join(var_tbl_add,
                   by = c("SubjectID", "VisitType"))
    } else if(names(params_dots$var_exp_spline) == 'risk11'){
      var_tbl_add <- TuberculomicsCompendium::signature_11gene %>%
        dplyr::select(SubjectID, VisitType, sig11gene_CorScore)
      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(params_dots$var_exp_spline))
      data_raw %<>%
        inner_join(var_tbl_add,
                   by = c("SubjectID", "VisitType"))
    } else if(names(params_dots$var_exp_spline) == 'hladr'){
      var_tbl_add <- dataset_list[[params_dots$dataset_name]]$hladr %>%
        dplyr::mutate(SubjectID = stringr::str_sub(fcs, end = 6),
               VisitType = purrr::map_chr(fcs, function(fcs_ind){
                 fcs_ind <- stringr::str_sub(fcs_ind, start = 8)
                 first_us_loc <- stringr::str_locate(fcs_ind, "_")[1,"start"][[1]]
                 stringr::str_sub(fcs_ind, end = first_us_loc - 1)
               })) %>%
        dplyr::select(SubjectID, VisitType, stim, hladr_med_diff)

      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(params_dots$var_exp_spline))

      data_raw %<>%
        inner_join(var_tbl_add,
                   by = c("SubjectID", "VisitType", "stim"))
    }

  }

  # ===========================
  # Make ttb at least a certain amount
  # ===========================

  if('min_ttb' %in% names(params_dots)){
    data_raw %<>%
      dplyr::mutate(timeToTBFromVisit = pmax(params_dots$min_ttb, timeToTBFromVisit))
  }

  # ===========================
  # Calculate time from max time-to-TB
  # ===========================

  if(!'max_ttb' %in% names(params_dots)){
    max_ttb <- rev(sort(data_raw$timeToTBFromVisit[!is.na(data_raw$timeToTBFromVisit)]))[min(10, nrow(data_raw %>%
                                                                                                        dplyr::filter(!is.na(timeToTBFromVisit))))]
  } else max_ttb <- params_dots$max_ttb

  data_raw %<>%
    dplyr::mutate(tfmttb = ifelse(Progressor == 'no', # make zero non-progressor
                                        0,
                                        -9999),
           tfmttb = ifelse(.data$Progressor == 'yes', # make 0 if time to TB is greater than max_ttb; otherwise make
                                        # how many days it is from collapse time to the actual time
                                        pmax(0, .env$max_ttb - .data$timeToTBFromVisit),
                                        .data$tfmttb)) %>%
    dplyr::select(SubjectID:timeToTBFromVisit, tfmttb, everything())


  # ===========================
  # Scale
  # ===========================

  if('scale' %in% names(params_dots)){
    if(params_dots$scale){
      #if("summed_resp" %in% colnames(data_raw)) data_raw %<>%
      #  dplyr::mutate(summed_resp = summed_resp/1e6)
      if("tfmttb" %in% colnames(data_raw)) data_raw %<>%
        dplyr::mutate(origMeantfmttb = mean(tfmttb),
               tfmttb = tfmttb/1e2)
      #if(resp_sd %in% names(params_dots)){
      #  data_raw %<>% dplyr::mutate(resp = resp/params_dots$resp_sd)
      #  }
      #}
      if('DaysSinceEntry' %in% colnames(data_raw)){
        data_raw %<>% dplyr::mutate(DaysSinceEntry = (DaysSinceEntry-mean(DaysSinceEntry))/sd(DaysSinceEntry))
      }
    }
  }

  # ===========================
  # Shift COMPASS values off zero and 1, if needed
  # ===========================

  if(params_dots$assay %in% c('cytof', 'cytof_ag_flowsom') && params_dots$response_type %in% c('compass', 'prob')){
    data_raw %<>% dplyr::mutate(resp = pmax(resp, min(data_raw$resp[data_raw$resp > 0])),
                         resp = pmin(resp, max(data_raw$resp[data_raw$resp < 1])))
  }

  if(params_dots$response_type %in% c('hladr_med_diff')){
    boost <- 0.01 * (diff(range(data_raw$resp)))
  #  boot <- boost + rnorm(nrow(data_raw), sd = boost/3)
    data_raw %<>% dplyr::dplyr::mutate(resp = pmax(resp, boost))
  }

  # ===========================
  # Remove columns with NA entries
  # ===========================

  data_raw %<>% dplyr::select(-c(timeToTBFromVisit))

  if(params_dots$assay %in% c('cytof', 'cytof-faust')){
    data_raw %<>%
      dplyr::mutate(resp = pmax(0, resp))
  }

  # output
  if(params_dots$assay == 'compass' || params_dots$response_type == 'hladr'){
    return(data_raw %>%
             dplyr::select(fcs, SubjectID, VisitType, Progressor, tfmttb, everything()))
  }

  data_raw %>%
    dplyr::select(SubjectID, VisitType, Progressor, origMeantfmttb, tfmttb, everything())

}
