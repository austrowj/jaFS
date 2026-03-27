duerte_ferreira <- function(df) WinRatio::winratio(
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

cui_huang <- function(df) WINS::win.stat(
    data =
        df
        |> dplyr::mutate(id = as.integer(substr(USUBJID, 13, 18)))
        |> dplyr::select(
            id,
            arm = trt,
            Delta_1 = DEATH,
            Delta_2 = MI,
            Delta_3 = STROKE,
            Y_1 = ttdeath,
            Y_2 = ttmi,
            Y_3 = ttstroke
        )
    ,
    ep_type = 'tte',
    arm.name = c("Active", "Placebo"),
    tau = 0,
    np_direction = "larger",
    priority = c(1:3),
    alpha = 0.05,
    digit = 3,
    stratum.weight = "unstratified",
    method = "unadjusted",
    pvalue = "two-sided"
)

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
        res <- cui_huang(df)
        end <- Sys.time()
        print("Done.")
        print(res)
        print(summary(res, digits = 10))

        print(end - start)
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

# Organize results into a dataframe (one row per trial) and write to a CSV file
results <- data.frame()

for (dup_factor in c(0)) {
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

if (FALSE) {
    filename <- paste0("benchmark_results_r_tmp.csv")
    write.csv(results, filename, row.names = FALSE)
    print(paste0("Results written to ", filename))
}