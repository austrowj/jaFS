#pyright: basic

import polars as pl
from polars_readstat import ScanReadstat, scan_readstat

import env

def main():
    ad_hx = scan_readstat(env.DATA_ROOT / 'ad_hx.sas7bdat')
    ad_safety = scan_readstat(env.DATA_ROOT / 'ad_safety.sas7bdat')

    df = (
        ad_hx
        .select(
            'USUBJID',
            pl.when(pl.col('RXGRP') == pl.lit(1)).then(pl.lit('Active')).otherwise(pl.lit('Placebo')).alias('trt'),
            pl.col('BMICAT_B').alias('BMICAT')
        )
        .join(on = 'USUBJID', how = 'right', other = 
            ad_safety
            .filter(pl.col('SAFEFLG') == 1)
            .select(
                'USUBJID',
                pl.col('TTDEATH').alias('censor'),
                'DEATH',
                'NFMI',
                'NFSTROKE',
                pl.col('DEATH').mul(pl.col('TTDEATH')).alias('ttdeath'),
                pl.col('NFMI').mul(pl.col('TTNFMI')).alias('ttmi'),
                pl.col('NFSTROKE').mul(pl.col('TTNFSTROKE')).alias('ttstroke'),
            )
        )

        # Convert to a timeline of events.
        .unpivot(
            ['ttdeath', 'ttmi', 'ttstroke', 'censor'],
            index = ['USUBJID', 'trt', 'BMICAT'],
            variable_name = 'event',
            value_name = 'time_to',
        )
        .filter(pl.col('time_to').gt(0))
        .sort('time_to')
    )
    print(df.collect())

if __name__ == '__main__':
    main()
