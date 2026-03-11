
using_study_data <- function() {
    # Load ad_hx from the data root directory specified in the .env file
    dotenv::load_dot_env(file.path(getwd(), ".env"))
    data_root <- Sys.getenv("DATA_ROOT")

    df <- (
        haven::read_sas(file.path(data_root, "ad_hx.sas7bdat"))
        |> dplyr::mutate(trt = dplyr::if_else(RXGRP == 1, "Active", "Placebo"))
        |> dplyr::select(USUBJID, trt, AGE, bmicat := BMICAT_B)
        |> dplyr::mutate(bmicat = tidyr::replace_na(bmicat, 99))

        |> dplyr::right_join(
            haven::read_sas(file.path(data_root, "ad_safety.sas7bdat"))
            |> dplyr::filter(SAFEFLG == 1)
            |> dplyr::mutate(
                censor = TTDEATH,
                ttdeath = DEATH * TTDEATH,
                ttmi = NFMI * TTNFMI,
                ttstroke = NFSTROKE * TTNFSTROKE
            )
            |> dplyr::select(USUBJID, censor, ttdeath, ttmi, ttstroke, DEATH, MI := NFMI, STROKE := NFSTROKE)
        )
    )
    
    # Duplicate all rows in df with new usubjids to test performance of winratio with larger datasets
    #df <- dplyr::bind_rows(df, df |> dplyr::mutate(USUBJID = paste0(USUBJID, "_dup")))
    
    print("Loaded and prepared starting data:")
    print(df)

    print("Computing win ratio...")
    res <- WinRatio::winratio(
        id = "USUBJID",
        trt = "trt",
        active = "Active",
        outcomes = list(
            d = c("DEATH", "s", "ttdeath"),
            mi = c("MI", "s", "ttmi"),
            #bmi = c("bmicat", "c", "<"),
            stroke = c("STROKE", "s", "ttstroke")
        ),
        fu = "censor",
        data = df,
    )
    print("Done.")
    print(summary(res))

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

using_study_data()

