#pyright: basic
from typing import Any, TypeVar, Callable

from sortedcontainers import SortedList

import polars as pl
from polars_readstat import ScanReadstat, scan_readstat

import env
from data_structures import FenwickTreeWithSortedLists

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
    win_ratio(
        df,
        lambda: (0, 0, 0),
        {
            'death':  lambda value, old: (value, old[1], old[2]),
            'mi':     lambda value, old: (old[0], value, old[2]),
            'stroke': lambda value, old: (old[0], old[1], value),
        }
    )

T = TypeVar('T')
def win_ratio[T](
        timeline: pl.DataFrame,
        default_score: Callable[[], T],
        update_score: dict[str, Callable[[Any, T], T]]
) -> None:
    """Calculate the win ratio for the given timeline of events."""
    import timeit

    start = timeit.default_timer()

    # Start by inserting the initial lexicographical comparison vector for each subject.
    current_score = {} # Map from USUBJID to current score vector.
    scoreboards = [SortedList(), SortedList()] # One for each treatment group.
    wins, losses, ties = [0, 0], [0, 0], [0, 0] # For active group (number 1)

    for row in timeline.unique('USUBJID', keep = 'last').iter_rows(named = True):
        usubjid = row['USUBJID']
        active: int = row['active']
        current_score[usubjid] = default_score()

        # Insert the new event into the appropriate sorted list.
        scoreboards[active].add(current_score[usubjid])
    
    subject_wins = {id: 0 for id in current_score.keys()}
    subject_losses = {id: 0 for id in current_score.keys()}
    
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
        if event != 'censor':
            current_score[usubjid] = update_score[event](time_to, current_score[usubjid])
            scoreboards[active].add(current_score[usubjid])

        else:
            # Do not re-add the subject. Instead, find their position in the opposite scoreboard and report that as number of wins.
            r_index = scoreboards[1 - active].bisect_right(current_score[usubjid])
            l_index = scoreboards[1 - active].bisect_left(current_score[usubjid])

            # The number of wins for this subject is the number of subjects in the opposite group
            # that have a lower score vector.
            win_count = len(scoreboards[1 - active]) - r_index
            tie_count = r_index - l_index
            loss_count = l_index

            # Update treatment-level wins and losses (and ties)
            wins[active] += win_count
            ties[active] += tie_count
            losses[active] += loss_count

            # Update subject-level wins
            subject_wins[usubjid] += win_count
            subject_losses[usubjid] += loss_count
    
    end = timeit.default_timer()
    print(f'Done.')

    # Print results.
    for i in range(2):
        print(f'{['Placebo', ' Active'][i]}: {wins[i]} wins, {losses[i]} losses, {ties[i]} ties, {wins[i] + losses[i] + ties[i]} total.')

    print(f'Total for Active: {wins[1] + losses[0]} wins, {wins[0] + losses[1]} losses, {sum(ties)} ties, {sum(wins) + sum(losses) + sum(ties)} total.')
    print(f'Time: {end - start} seconds.')      

if __name__ == '__main__':
    main()
