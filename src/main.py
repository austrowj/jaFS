#pyright: basic
from math import inf, log, sqrt
from typing import Any, TypeVar, Callable

from sortedcontainers import SortedList

import polars as pl
from polars_readstat import ScanReadstat, scan_readstat

import env
from data_structures import FenwickTreeWithSortedLists

def main(dup_factor = 0, trials = 1):
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
    # Duplicate data with new usubjids to test performance with larger datasets
    dups = (df.with_columns(pl.col('USUBJID') + pl.lit(f'_dup{i}')) for i in range(dup_factor))
    df = pl.concat([df] + list(dups)).sort('time_to', pl.col('event').eq(pl.lit('censor')))
    print(df)
    
    times = []
    for i in range(trials):
        print(f'Computing win ratio ({i+1}/{trials})...')
        times.append(win_ratio(
            df,
            lambda: (0, 0, 0),
            {
                'death':  lambda value, old: (value, old[1], old[2]),
                'mi':     lambda value, old: (old[0], value, old[2]),
                'stroke': lambda value, old: (old[0], old[1], value),
            }
        ))
    return times

T = TypeVar('T')
def win_ratio[T](
        timeline: pl.DataFrame,
        default_score: Callable[[], T],
        update_score: dict[str, Callable[[Any, T], T]]
) -> float:
    """Calculate the win ratio for the given timeline of events. Return time taken."""
    import timeit

    start = timeit.default_timer()

    # Start by inserting the initial lexicographical comparison vector for each subject.
    current_score = {} # Map from USUBJID to current score vector.
    scoreboards = [SortedList(), SortedList()] # One for each treatment group.
    favor_active, favor_placebo, favor_neither = 0, 0, 0

    trt: list[list[str]] = [list(), list()]

    for row in timeline.unique('USUBJID', keep = 'last').iter_rows(named = True):
        usubjid = row['USUBJID']
        active: int = row['active']
        current_score[usubjid] = default_score()

        trt[active].append(usubjid)

        # Insert the new event into the appropriate sorted list.
        scoreboards[active].add(current_score[usubjid])

    subject_data = {id: [0, 0, 0] for id in current_score.keys()} # wins, losses, last seen time
    trees = [FenwickTreeWithSortedLists(len(timeline)), FenwickTreeWithSortedLists(len(timeline))]

    # Save basic summary stats
    n_active = len(scoreboards[1])
    n_placebo = len(scoreboards[0])
    n = n_active + n_placebo
    
    print(f'Number of subjects: {n}')
    print(f'Number of subjects in \'Active\' group: {n_active}')
    print(f'Number of subjects in \'Placebo\' group: {n_placebo}')
    print(f'Number of paired comparisons: {n_active * n_placebo}')

    # Process entire timeline of events in order.
    for i, row in enumerate(timeline.iter_rows(named = True)):
        usubjid = row['USUBJID']
        active: int = row['active']
        event: str = row['event']
        time_to: int = row['time_to']

        tree_time = i + 1 # Fenwick trees are 1-indexed.

        # Remove the current score from the scoreboard.
        scoreboards[active].remove(current_score[usubjid])

        # First, check fenwick tree of other treatment group for updates since last seen time.
        update_losses, update_wins = trees[1-active].query_range(current_score[usubjid], subject_data[usubjid][2], tree_time)
        subject_data[usubjid][0] += update_wins
        subject_data[usubjid][1] += update_losses
        subject_data[usubjid][2] = tree_time

        # Update the score vector based on the event type.
        if event != 'censor':
            current_score[usubjid] = update_score[event](10000 - time_to, current_score[usubjid]) # later events are better
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
            favor_neither += tie_count
            if (active):
                favor_active += win_count
                favor_placebo += loss_count
            else:
                favor_active += loss_count
                favor_placebo += win_count

            # Update subject-level Uw
            subject_data[usubjid][0] += win_count
            subject_data[usubjid][1] += loss_count

            # Add to treatment group fenwick tree
            trees[active].add(current_score[usubjid], tree_time)
    
    end = timeit.default_timer()
    print(f'Done.')

    # Compute Z-score and confidence interval.
    ttw = favor_active / (n*n)
    ttl = favor_placebo / (n*n)
    wr = ttw/ttl if ttl > 0 else -inf
    print(f'Win ratio: {wr}')

    record = (
        timeline
        .group_by('USUBJID', 'active').agg(pl.len().alias('num_events'))
        .join(on = 'USUBJID', how = 'inner', other =
            pl.DataFrame({
                'USUBJID': list(subject_data.keys()),
                'wins': [stat[0] for stat in subject_data.values()],
                'losses': [stat[1] for stat in subject_data.values()],
                'net': [stat[0] - stat[1] for stat in subject_data.values()]
            })
        )
        .sort('USUBJID')
    )

    trt_record = record.filter(pl.col('active').eq(pl.lit(1)))
    pla_record = record.filter(pl.col('active').eq(pl.lit(0)))

    # Results for WR, V, and Z agree with Duarte&Ferreira R package to >15 decimal places
    n0 = n_placebo
    n1 = n_active
    v = (
        sum(
            (
                (
                    trt_record['wins'] - trt_record['losses'] * wr
                )
                / n0
            ) ** 2
        )
        * (n0/n)**2 / n
        + sum(
            (
                (
                    pla_record['losses'] - pla_record['wins'] * wr
                )
                / n1
            )
            ** 2
        )
        * (n1/n)**2 / n
    ) / (ttl**2)
    z = sqrt(n)*log(wr)*wr/sqrt(v)
    #print(f'V: {v}')
    #print(f'Z: {z}')

    #print(f'{favor_active} wins, {favor_placebo} losses, {favor_neither} ties.')
    #print(f'Time: {end - start} seconds.')

    return end - start

if __name__ == '__main__':
    import csv
    with open('out/results_py_2.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['dup_factor', 'trial', 'time'])
        for dup_factor in [0, 1, 2, 3, 4]:
            times = main(dup_factor = dup_factor, trials = 10)
            print(f'Duplication factor: {dup_factor}, times: {times}')
            for trial, time in enumerate(times):
                writer.writerow([dup_factor, trial, time])
