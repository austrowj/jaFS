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
            pl.when(pl.col('RXGRP') == pl.lit(1)).then(pl.lit(1)).otherwise(pl.lit(0)).alias('active'),
            pl.col('BMICAT_B').fill_null(99).alias('BMICAT')
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
            index = ['USUBJID', 'active',],# 'BMICAT'],
            variable_name = 'event',
            value_name = 'time_to',
        )
        .filter(pl.col('time_to').gt(0))
        .sort('time_to', pl.col('event').eq(pl.lit('censor')))
    )
    df = df.collect()
    print(df)

    print('Computing win ratio...')
    win_ratio(df)

def win_ratio(timeline: pl.DataFrame) -> None:
    """Calculate the win ratio for the given timeline of events."""
    from sortedcontainers import SortedList

    # Start by inserting the initial lexicographical comparison vector for each subject.
    current_score = {} # Map from USUBJID to current score vector.
    scoreboards = [SortedList(), SortedList()] # One for each treatment group.
    wins, losses, ties = [0, 0], [0, 0], [0, 0] # For active group (number 1)

    for row in timeline.filter(pl.col('event') == pl.lit('censor')).iter_rows(named = True):
        usubjid = row['USUBJID']
        active: int = row['active']
        current_score[usubjid] = (0, 0, 0) # No events recorded

        # Insert the new event into the appropriate sorted list.
        scoreboards[active].add(current_score[usubjid])
    
    print(f'Number of subjects: {len(scoreboards[0]) + len(scoreboards[1])}')
    print(f'Number of subjects in \'Active\' group: {len(scoreboards[1])}')
    print(f'Number of subjects in \'Placebo\' group: {len(scoreboards[0])}')
    print(f'Number of paired comparisons: {len(scoreboards[0]) * len(scoreboards[1])}')

    # Process entire timeline of events in order.
    for row in timeline.iter_rows(named = True):
        usubjid = row['USUBJID']
        active: int = row['active']
        event: str = row['event']

        # Remove the current score from the scoreboard.
        scoreboards[active].remove(current_score[usubjid])

        # Update the score vector based on the event type.
        if event == 'ttdeath':
            current_score[usubjid] = (1, current_score[usubjid][1], current_score[usubjid][2])
        elif event == 'ttmi':
            current_score[usubjid] = (current_score[usubjid][0], 1, current_score[usubjid][2])
        elif event == 'ttstroke':
            current_score[usubjid] = (current_score[usubjid][0], current_score[usubjid][1], 1)

        elif event == 'censor':
            # Do not re-add the subject. Instead, find their position in the opposite scoreboard and report that as number of wins.
            r_index = scoreboards[1 - active].bisect_right(current_score[usubjid])
            l_index = scoreboards[1 - active].bisect_left(current_score[usubjid])
            # The number of wins for this subject is the number of subjects in the opposite group
            # that have a lower score vector.
            losses[active] += l_index
            wins[active] += len(scoreboards[1 - active]) - r_index
            ties[active] += r_index - l_index

            continue # Do not re-add the subject to the scoreboard.

        # Re-insert the updated score into the scoreboard.
        scoreboards[active].add(current_score[usubjid])
    
    # Print results.
    print(f'Active: {wins[1]} wins, {losses[1]} losses, {ties[1]} ties')
    print(f'Total: {wins[1] + losses[1] + ties[1]}')

if __name__ == '__main__':
    main()
