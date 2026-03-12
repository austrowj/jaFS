#pyright: basic

import polars as pl
from polars_readstat import ScanReadstat, scan_readstat

import env

def main():
    print('Loading data...')
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
            #.slice(20, 10)
            .select(
                'USUBJID',
                pl.col('TTDEATH').alias('censor'),
                'DEATH',
                'NFMI',
                'NFSTROKE',
                pl.col('DEATH').mul(pl.col('TTDEATH')).alias('death'),
                pl.col('NFMI').mul(pl.col('TTNFMI')).alias('mi'),
                pl.col('NFSTROKE').mul(pl.col('TTNFSTROKE')).alias('stroke'),
            )
        )

        # Convert to a timeline of events.
        .unpivot(
            ['death', 'mi', 'stroke', 'censor'],
            index = ['USUBJID', 'active'],
            variable_name = 'event',
            value_name = 'time_to',
        )
        .filter(pl.col('time_to').gt(0))
        #.filter(pl.col('event').is_in(['censor', 'death']))
        .sort('time_to', pl.col('event').eq(pl.lit('censor')))
    )
    df = df.collect()
    print(df)

    print('Computing win ratio...')
    win_ratio(df)

def win_ratio(timeline: pl.DataFrame) -> None:
    """Calculate the win ratio for the given timeline of events."""
    from sortedcontainers import SortedList
    import timeit

    start = timeit.default_timer()

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
    for i, row in enumerate(timeline.iter_rows(named = True)):
        usubjid = row['USUBJID']
        active: int = row['active']
        event: str = row['event']
        time_to: int = row['time_to']

        # Remove the current score from the scoreboard.
        scoreboards[active].remove(current_score[usubjid])

        # Update the score vector based on the event type.
        if event == 'death':
            current_score[usubjid] = (time_to, current_score[usubjid][1], current_score[usubjid][2])
        elif event == 'mi':
            current_score[usubjid] = (current_score[usubjid][0], time_to, current_score[usubjid][2])
        elif event == 'stroke':
            current_score[usubjid] = (current_score[usubjid][0], current_score[usubjid][1], time_to)

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

        else:
            raise ValueError(f'Unrecognized event type "{event}".')

        # Re-insert the updated score into the scoreboard.
        scoreboards[active].add(current_score[usubjid])
    
    end = timeit.default_timer()
    print(f'Done.')

    # Print results.
    for i in range(2):
        print(f'{['Placebo', ' Active'][i]}: {wins[i]} wins, {losses[i]} losses, {ties[i]} ties, {wins[i] + losses[i] + ties[i]} total.')

    print(f'Total for Active: {wins[1] + losses[0]} wins, {wins[0] + losses[1]} losses, {sum(ties)} ties, {sum(wins) + sum(losses) + sum(ties)} total.')
    print(f'Time: {end - start} seconds.')
        

if __name__ == '__main__':
    main()
