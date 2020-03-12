import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from collections import OrderedDict

def normalize_read_counts(meta, df_freqs, df_num_reads_per_sample, outpath):
    #################### Normalized Read Counts #######################################################
    # Create an additional table of normalized read counts based on average library sequencing depth
    # get the mapping between samples and libraries from meta data
    library_dict = meta.generate_dict('Sequencing library #', tolower=0)
    df_library = pd.DataFrame([library_dict]).transpose().reindex(df_num_reads_per_sample.transpose().index)
    df_library.columns = {'Sequencing library #'}

    average_mapped_reads_per_library = df_num_reads_per_sample.groupby(by=library_dict, axis=1, sort=False).mean()
    overall_average_reads_per_library = average_mapped_reads_per_library.mean(axis=1)
    normalization_factor = average_mapped_reads_per_library.apply(
        lambda x: overall_average_reads_per_library * (1 / x)).transpose()
    # expand normalization factor to a per sample vector
    normalization_factor_per_sample = df_library.merge(normalization_factor, left_on='Sequencing library #',
                                                       right_index=True, how='left')
    normalization_factor_per_sample.drop('Sequencing library #', axis=1, inplace=True)

    # move data from frequencies to num reads:
    df_reads = convert_freqs_to_num_reads(df_freqs, df_num_reads_per_sample)
    df_reads_normalized = pd.DataFrame(
        df_reads.as_matrix() * normalization_factor_per_sample.transpose().as_matrix().astype(float),
        columns=df_reads.columns)
    df_reads_normalized.index = df_reads.index

    # Generate a plot of different numbers of reads per sample - original, mapped, normalized, normalization factor, sample library
    fig, axarr = plt.subplots(4,1)
    axarr[0].plot(df_reads.sum().values, color='blue', label='mapped')
    axarr[1].plot(df_reads_normalized.sum().values, color='green', label='normalized mapped')
    axarr[2].plot(normalization_factor_per_sample.values, color='purple', label='normalizing factor')
    axarr[3].plot(df_library.values, color='black', label='library #')
    axarr[0].legend()
    axarr[1].legend()
    axarr[2].legend()
    axarr[3].legend()
    plt.savefig(os.path.join(outpath, 'Read normalization history'))

    return df_reads, df_reads_normalized

#convert frequency table to num reads table
#multiply each cell by the number of reads of the column
def convert_freqs_to_num_reads(df_freqs, df_num_reads_per_sample):
    df_temp = df_freqs.copy()
    df_temp.loc[:] = df_freqs.as_matrix() * df_num_reads_per_sample.as_matrix().astype(float)
    df_temp.index = df_freqs.index
    return df_temp

def agg_allsame(x):
    if (x == x.iloc[0]).all():
        return x.iloc[0]
    else:
        return 'AMB'

# input can be normalized reads or frequencies
# data is summed across taxonomic levels
def generate_tables_per_taxonomic_level(df_reads, df_taxonomies):
    df_new = pd.DataFrame()

    orig_columns = df_taxonomies.columns.append(df_reads.columns)

    # propogate information up from species level to all higher order taxonomies
    taxon_levels = {0: 'domain', 1: 'phylum', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6: 'species'}
    # taxon_levels = {4: 'family', 5: 'genus'} #for DEBUG, just look at a couple lower tax levels
    offset_levels = {'domain': 80000, 'phylum': 70000, 'class': 60000, 'order': 50000, 'family': 40000, 'genus': 30000}
    # agg_for_taxons = {'domain':agg_allsame, 'phylum':agg_allsame, 'class':agg_allsame, 'order':agg_allsame, 'family':agg_allsame, 'genus':agg_allsame}

    # I originally wanted to groupby taxon level and then aggregate the sample data by sum and the taxon data by agg_allsame
    # but I couldn't get it to work in one shot.
    # d = dict(zip(df_reads.columns, np.repeat('sum', len(df_reads.columns))))
    # agg_for_taxons.update(d)

    # add taxonomic prefix to data and drop species information which we don't need for upper layers
    df_taxonomies_copy = df_taxonomies.drop('species', axis=1)
    df_merged = df_taxonomies_copy.merge(df_reads, left_index=True, right_index=True, how='inner')

    higher_order_taxon_levels = []
    for index, taxon_level in taxon_levels.iteritems():
        if taxon_level == 'species': #skip species
            continue

        # reset agg function list
        agg_for_taxons = OrderedDict(
            [('domain', agg_allsame), ('phylum', agg_allsame), ('class', agg_allsame), ('order', agg_allsame),
             ('family', agg_allsame), ('genus', agg_allsame)])

        # group by desired taxon level with hierarchy
        higher_order_taxon_levels.append(taxon_level)
        df_merged_agg = df_merged.groupby(higher_order_taxon_levels, axis=0).sum()


        df_merged_agg.reset_index(inplace=True)
        #expand taxonomies to full set of 7

        for j in np.arange(index + 1, max(taxon_levels) + 1):
            df_merged_agg = pd.concat([df_merged_agg, pd.Series(np.tile(np.nan, len(df_merged_agg.index)), name=taxon_levels[j], dtype='object')], axis=1)

        # add an offset to the taxonomic level for easy identification (Ravid's idea)
        df_merged_agg.index = df_merged_agg.index + offset_levels[taxon_level]

        # concatenate the tables for all the levels
        df_new = pd.concat([df_new, df_merged_agg], axis=0, sort=True)

    # split again
    # reindex the colums after all the concatenation to be in the 'right' order
    df_new = df_new.reindex_axis(orig_columns, axis=1)
    df_taxonomies_new = df_new.select_dtypes(include={'object'})
    df_reads_new = df_new.select_dtypes(include={'float64'})

    return df_reads_new, df_taxonomies_new


