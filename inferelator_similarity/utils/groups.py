from .constants import GROUP_COL, SPECIES_COL

def group_stats(
    df,
    group_col=GROUP_COL,
    stat_col=SPECIES_COL,
    stats=('count', 'nunique')
):

        return df.groupby(
            group_col
        )[stat_col].agg(
            stats
        )
