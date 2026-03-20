
using_study_data <- function(dup_factor = 0, trials = 1) {
    # Load ad_hx from the data root directory specified in the .env file
    dotenv::load_dot_env(file.path(getwd(), ".env"))
    data_root <- Sys.getenv("DATA_ROOT")

    print("Loading data...")

    df <- (
        haven::read_sas(file.path(data_root, "ad_hx.sas7bdat"))
        |> dplyr::mutate(trt = dplyr::if_else(RXGRP == 1, "Active", "Placebo"))
        |> dplyr::select(USUBJID, trt, AGE, bmicat := BMICAT_B)
        |> dplyr::mutate(bmicat = tidyr::replace_na(bmicat, 99))

        |> dplyr::right_join(
            haven::read_sas(file.path(data_root, "ad_safety.sas7bdat"))
            #|> dplyr::slice(21:30)
            |> dplyr::filter(SAFEFLG == 1)
            |> dplyr::mutate(censor = TTDEATH)
            |> dplyr::mutate(
                ttdeath = TTDEATH,
                ttmi = dplyr::if_else(NFMI == 1, TTNFMI, censor),
                ttstroke = dplyr::if_else(NFSTROKE == 1, TTNFSTROKE, censor)
            )
            |> dplyr::select(USUBJID, censor, ttdeath, ttmi, ttstroke, DEATH, MI := NFMI, STROKE := NFSTROKE)
        )
    )
    
    # Duplicate all rows in df with new usubjids to test performance of winratio with larger datasets
    if (dup_factor > 0) {
        dups <- lapply(1:dup_factor, function(i) {
            df |> dplyr::mutate(USUBJID = paste0(USUBJID, "_dup", i))
        })
        df <- dplyr::bind_rows(df, dups)
    }

    print("Loaded and prepared starting data:")
    print(df)

    times <- c()
    for (i in 1:trials) {
        print("Computing win ratio...")
        start <- Sys.time()
        res <- WinRatio::winratio(
            id = "USUBJID",
            trt = "trt",
            active = "Active",
            outcomes = list(
                d = c("DEATH", "s", "ttdeath")
                ,mi = c("MI", "s", "ttmi")
                #,bmi = c("bmicat", "c", "<")
                ,stroke = c("STROKE", "s", "ttstroke")
            ),
            fu = "censor",
            data = df
            ,keep.matrix = TRUE
        )
        end <- Sys.time()
        print("Done.")
        print(summary(res, digits = 10))

        print(paste0("Time: ", end - start, " seconds."))
        print(paste0("V: ", res$v))
        print(paste0("Z: ", res$z))

        mat <- res$wr.matrix
        print(sum(apply(mat > 0, 1, sum) - apply(mat < 0, 1, sum) * res$wr))
        print(sum(apply(mat > 0, 2, sum) - apply(mat < 0, 2, sum) * res$wr))
        #print(apply(mat < 0, 2, sum))
        print(sum(res$loss)/(res$n^2))

        times <- c(times, end - start)
    }
    return(times)
}

using_test_survival_data <- function() {
    library(dplyr)
    require("survival")
    # Creation of dataset 'df' with 3 outcomes:
    # Outcome 1: death (survival event)
    # Outcome 2: cancer recurrence (repeated survival event)
    # Outcome 3: size of largest initial tumour (continuous event)
    data1 <- (
        survival::bladder1
        |> mutate(trt = if_else(treatment == "placebo", "Placebo", "Treatment"))
        |> group_by(id)
        |> mutate(
            death = if_else(max(status) %in% c(2, 3), 1, 0),
            t2death = max(stop)
        )
        |> ungroup()
        |> select(id, trt, death, t2death, number, size)
        |> unique()
    )
    print(data1)
}

# Organize results into a dataframe (one row per trial) and write to a CSV file
results <- data.frame()

for (dup_factor in c(5)) {
    print(paste0("Duplication factor: ", dup_factor))
    times <- using_study_data(dup_factor = dup_factor, trials = 1)
    trial_results <- data.frame(
        dup_factor = dup_factor,
        trial = seq_along(times),
        time_seconds = as.numeric(times)
    )
    results <- rbind(results, trial_results)
    print(times)
}

filename <- paste0("benchmark_results_r_dup5_5.csv")
write.csv(results, filename, row.names = FALSE)
print(paste0("Results written to ", filename))