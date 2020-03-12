from metadata_class_TMB import Metadata
import pandas as pd
import numpy as np
from scipy import stats
from Analysis_pipeline_extra_functions_TMB import normalize_read_counts, generate_tables_per_taxonomic_level
from taxonomy_TMB import Taxonomy
import os
import pickle
import matplotlib.pyplot as plt
from mne.stats import fdr_correction

class dataLoader:

    def __init__(self, results_file_dir):
        self.FILES_PATH = results_file_dir

    def load_data(self, metadata_file, relative_abundance_file, num_reads_file, source_files_dir, results_file_dir):
        # initialize these pieces of data
        self.INFILES_PATH = source_files_dir
        self.FILES_PATH = results_file_dir
        self.metadata_file = metadata_file
        self.relative_abundance_file = relative_abundance_file
        self.num_reads_file = num_reads_file

        self.meta = Metadata(os.path.join(self.INFILES_PATH, self.metadata_file))

        print "Reading frequency data from: " + self.relative_abundance_file
        df = pd.read_table(os.path.join(self.INFILES_PATH, self.relative_abundance_file), sep='\t')

        print "Reading number of reads per sample from: " + self.num_reads_file
        f = open(os.path.join(self.INFILES_PATH, self.num_reads_file))
        num_reads_per_sample = f.readline().strip().split("\t")
        f.close()

        self.df_taxonomies = df.select_dtypes(include={'object'})
        self.df_freqs = df.select_dtypes(include={'float64'})

        self.df_num_reads_per_sample = pd.DataFrame(data=num_reads_per_sample).transpose().astype('float')
        self.df_num_reads_per_sample.columns = self.df_freqs.columns
        self.df_reads_normalized = pd.DataFrame()

        print 'Removing irrelevant samples'
        self.select_samples(self.df_freqs.columns[self.df_freqs.columns.isin(self.meta.df_meta.columns)])


    def get_outfiles_path(self):
        return self.FILES_PATH

    #this might be dangerous because it exposes the metadata object, but there is a class that can be used to manipulate it carefully
    def get_metadata_object(self):
        return self.meta

    #use select_samples with IDs of samples that received a minimum number of reads
    def remove_min_read_samples(self, min_reads):
        print 'Removing samples with less than '+str(min_reads)+' reads'
        df_num_t = self.df_num_reads_per_sample.transpose()
        min_reads_samples = df_num_t[df_num_t[0] >= min_reads].index

        self.select_samples(min_reads_samples)

    #given sample IDs, remove them (generic function)
    def select_samples(self, sample_ids):
        self.meta.reduce(sample_ids)
        self.df_freqs = self.df_freqs[sample_ids]
        self.df_num_reads_per_sample = self.df_num_reads_per_sample[sample_ids]
        if self.df_reads_normalized.empty:
            pass
        else:
            self.df_reads_normalized = self.df_reads_normalized[sample_ids]

        assert (self.df_freqs.columns.isin(self.meta.df_meta.columns).sum() == len(self.df_freqs.columns))
        self.num_samples_orig = len(self.df_num_reads_per_sample.columns)
        print 'Number of samples remaining is %d' % (self.num_samples_orig)


    #use the actual values in the index (loc) and NOT the numeric location in the matrix (iloc)
    def select_rows(self, row_idx):
        # remove the 'all' nan lines from frequency table and global contaminants
        self.df_freqs = self.df_freqs.loc[row_idx]
        # do the same for the correlating taxonomy table
        self.df_taxonomies = self.df_taxonomies.loc[row_idx]
        ##### Should I adjust the total number of reads after removing these lines????
        if self.df_reads_normalized.empty:
            pass
        else:
            self.df_reads_normalized = self.df_reads_normalized.loc[row_idx]

        print 'Removing rows, number of remaining rows is '+str(len(row_idx))

    def floor_abundances(self, abundance_cutoff):
        print 'Flooring abundances at '+str(abundance_cutoff)
        self.df_freqs.replace(0, np.nan, inplace=True)
        # replace all values that are lower than the floor with nan
        floor_locations = (self.df_freqs < abundance_cutoff)

        self.df_freqs[floor_locations] = np.nan
        # mask the normalized reads by the same mask as the df_freqs floor mask - to be consistent
        if self.df_reads_normalized.empty:
            pass
        else:
            self.df_reads_normalized[floor_locations] = 0  # nan's aren't used in this table

        # drop rows that become all nans (zeros originally)
        idx = self.df_freqs.index[~(self.df_freqs.isnull().all(1))]
        self.select_rows(idx)


    def generate_normalized_read_counts(self):
        print 'Normalizing reads according library averages'
        df_reads, self.df_reads_normalized = normalize_read_counts(self.meta, self.df_freqs, self.df_num_reads_per_sample, self.FILES_PATH)


    def save_data(self, filename):
        #save binary object
        f = open(os.path.join(self.FILES_PATH, filename), 'wb')
        pickle.dump(self.__dict__, f)
        f.close()

        #save to text files
        self.df_freqs.to_csv(os.path.join(self.FILES_PATH, 'df_freqs.csv'))
        self.df_reads_normalized.to_csv(os.path.join(self.FILES_PATH, 'df_reads_normalized.csv'))
        self.df_taxonomies.to_csv(os.path.join(self.FILES_PATH, 'df_taxonomies.csv'))
        self.df_num_reads_per_sample.to_csv(os.path.join(self.FILES_PATH, 'df_num_reads_per_sample.csv'))
        df_reads_normalized_merged = self.df_taxonomies.merge(self.df_reads_normalized, left_index=True, right_index=True,
                                                         how='inner')
        df_reads_normalized_merged.to_csv(os.path.join(self.FILES_PATH, 'df_reads_normalized_with_taxonomy.csv'))


    def load_data_from_cache(self, filename):
        f = open(os.path.join(self.FILES_PATH, filename), 'rb')
        tmp = pickle.load(f)
        f.close()

        self.__dict__.update(tmp)


    #change df_freqs, df_reads_normalized, and df_taxonomies to include higher level taxonomies
    #assumption here is that df_reads_normalized exists already because we use it for assertions
    def bubble_up_taxonomic_levels(self):
        df_freqs_hier_taxons, df_taxonomies_hier_taxons = generate_tables_per_taxonomic_level(self.df_freqs, self.df_taxonomies)
        if self. df_reads_normalized.empty:
            print 'Error: df_reads_normalized is not initialized, not bubbling up taxonomies'
            return
        else:
            df_reads_normalized_hier_taxons, dummy = generate_tables_per_taxonomic_level(self.df_reads_normalized, self.df_taxonomies)

        # reset 0 frequencies to nans
        df_freqs_hier_taxons[df_freqs_hier_taxons == 0] = np.nan

        # assert that the columns are in the same order before concatenting instead of using sort=True because we
        # need to preserve the order so that other things downstream don't break
        assert (pd.DataFrame(zip(self.df_reads_normalized.columns, df_reads_normalized_hier_taxons.columns),
                                 columns=['species', 'higher_taxons']).apply(lambda x: x[0] == x[1],
                                                                             axis=1).sum() == len(df_reads_normalized_hier_taxons.columns))

        # reorder index of df_taxonomies to match standard index
        self.df_taxonomies = pd.concat([self.df_taxonomies, df_taxonomies_hier_taxons], axis=0, sort=False)
        self.df_reads_normalized = pd.concat([self.df_reads_normalized, df_reads_normalized_hier_taxons], axis=0,
                                                 sort=False)
        self.df_freqs = pd.concat([self.df_freqs, df_freqs_hier_taxons], axis=0, sort=False)

    #if the tables contain upper level taxonomies, return info for species level only
    def focus_on_species(self):
        species_ids = Taxonomy().filter_ids(self.df_freqs.index, 'species')
        df_freqs_tmp = self.df_freqs.loc[species_ids]
        df_taxonomies_tmp = self.df_taxonomies.loc[species_ids]
        if self.df_reads_normalized.empty:
            pass
        else:
            df_reads_normalized_tmp = self.df_reads_normalized.loc[species_ids]

        return df_taxonomies_tmp, df_freqs_tmp, df_reads_normalized_tmp

    def filter_environmental_contaminants(self, percentage_in_controls=0.1):
        df_contaminants, df_non_contaminants_mask, df_percentages = self.filter_contaminants([('control', 'control'), ('ntc', 'ntc')], 'all_controls', percentage_in_controls)
        return df_contaminants, df_non_contaminants_mask, df_percentages

    def filter_contaminants(self, conditions, desc, percentage_in_controls=0.1):
        meta = self.meta
        if isinstance(conditions, tuple):
            condition_ids = meta.get_condition_ids(conditions)
        else: # a list of condition tuples
            condition_ids = meta.get_condition_ids(conditions[0])
            for i in np.arange(1, len(conditions)):
                condition_ids = condition_ids.append(meta.get_condition_ids(conditions[i]))

        num_controls = len(condition_ids)
        df_controls = self.df_freqs[condition_ids]
        #strong assumption that cells are NaN when there is no bacterial frequency in a sample, this would not work for df_reads_normalized
        df_percentages = df_controls.count(axis=1) / num_controls

        ##########A little graphing#################
        tmp_for_hist = df_percentages.loc[Taxonomy().filter_ids(df_percentages.index, 'species')]
        tmp_for_hist.hist(bins=40)
        plt.xticks(np.arange(-0.025, 1, 0.025), rotation=90)
        plt.yticks(np.arange(-1000, 12000, 1000))
        plt.savefig(os.path.join(self.FILES_PATH, 'Species_distribution_in_'+desc))
        ############################################

        df_contaminants = df_percentages.map(lambda x: 1 if x >= percentage_in_controls else 0)
        #we don't use the upper level taxons for this analysis since the nature of contamination is by a specific bacterial species falling in the tube
        #so erase the marking of contaminations from upper level taxons
        higher_taxon_ids = np.setdiff1d(df_contaminants.index, Taxonomy().filter_ids(df_contaminants.index, 'species'))
        df_contaminants.loc[higher_taxon_ids] = 0

        #by this time, higher level taxons are already zero-ed out so we don't need to filter for species only, the comparison to 1 takes care of it
        self.df_taxonomies.loc[df_contaminants == 1].to_csv(os.path.join(self.FILES_PATH, 'df_contaminants_'+desc+'_'+str(percentage_in_controls)+'.csv'))

        print 'Number of '+desc+' species contaminants that appear in over '+str(100*percentage_in_controls)+'% of controls is: '+str(df_contaminants.sum())
        #flip the mask
        df_non_contaminants_mask = df_contaminants.map(lambda x: 0 if x==1 else 1)

        return df_contaminants, df_non_contaminants_mask, df_percentages

class hitCaller:

    def __init__(self, dataLoader_object):
        self.dataLoader_object = dataLoader_object
        self.meta = self.dataLoader_object.get_metadata_object()
        self.batch_types = ('DNA Extraction Batch', 'PCR - New Batch', 'Sequencing library #')

    def save_hit_caller(self, filename):
        #save binary object
        f = open(os.path.join(self.dataLoader_object.FILES_PATH, filename), 'wb')
        pickle.dump(self.__dict__, f)
        f.close()

    def save_hits(self, hits_table):
        f = open(os.path.join(self.dataLoader_object.FILES_PATH, 'hits'), 'wb')
        pickle.dump(hits_table, f)
        f.close()

    def load_hits(self):
        f = open(os.path.join(self.dataLoader_object.FILES_PATH, 'hits'), 'rb')
        tmp = pickle.load(f)
        f.close()

        return tmp

    def load_hit_caller_from_cache(self, filename):
        f = open(os.path.join(self.dataLoader_object.FILES_PATH, filename), 'rb')
        tmp = pickle.load(f)
        f.close()

        self.__dict__.update(tmp)


    def save_filtered_data(self, filter_result, filename):
        f = open(os.path.join(self.dataLoader_object.FILES_PATH, filename), 'wb')
        pickle.dump(filter_result, f)
        f.close()

    def load_filtered_data(self, filename):
        f = open(os.path.join(self.dataLoader_object.FILES_PATH, filename), 'rb')
        tmp = pickle.load(f)
        f.close()
        return tmp

    def hit_call(self, meta, df_freqs):
        desc1 = 'batch'
        condition_stats_batches = self.batch_filters(meta, df_freqs)
        desc2 = 'centers'
        condition_stats_centers = self.center_filters(meta, df_freqs, 2, pval_cutoff=0.05)
        conditions_stats_combined = self.combine_results(condition_stats_batches, condition_stats_centers,
                                                               'intersect', desc1, desc2)
        desc3 = 'paraffin'
        condition_stats_paraffins = self.paraffin_filters(meta, df_freqs, pval_cutoff=0.05)
        conditions_stats_combined_with_paraffin = self.combine_results(conditions_stats_combined,
                                                                             condition_stats_paraffins,
                                                                             'intersect', desc1 + '_' + desc2, desc3)
        return conditions_stats_combined_with_paraffin

    def hit_call_new(self, meta, df_freqs):
        desc1 = 'new'
        condition_stats_new = self.new_filters(meta, df_freqs, pval_cutoff=0.05,
                                                     fdr_cutoff=0.20, how_many_centers=2)
        desc2 = 'paraffin'
        #trying to test for paraffin hits in every iteration
        #the other option is to set it up once and just intersect with the condition hit calling
        condition_stats_paraffins = self.paraffin_filters(meta, df_freqs, pval_cutoff=0.05)
        conditions_stats_combined_with_paraffin = self.merge_results(condition_stats_new,
                                                                           condition_stats_paraffins, desc1, desc2)
        return conditions_stats_combined_with_paraffin


    def get_hit_ids(self, condition_stats):
        return condition_stats.index[condition_stats.any(axis=1)]

    def new_filters(self, meta, df_freqs, pval_cutoff = 0.05, fdr_cutoff = 0.20, how_many_centers=2):
        # get sliced table per condition
        conditions = meta.get_conditions()
        conditions = filter(lambda x: (x != ('control', 'control') and (x != ('ntc', 'ntc'))), conditions)

        num_controls_total = len(meta.get_control_ids(with_NTCs=1))
        epsilon_p = 1.0 / num_controls_total

        df_filter_tracker_raw = pd.DataFrame()
        df_filter_tracker = pd.DataFrame()

        condition_stats = pd.DataFrame()
        for c in conditions:
            print c
            condition_ids = meta.get_condition_ids(c)
            df_freqs_c = df_freqs[condition_ids]
            df_pvals_c = pd.DataFrame()

            for batch_type in self.batch_types:
                batch_list = list(set(list(meta.df_meta.loc[batch_type, condition_ids])))

                if batch_type == 'DNA Extraction Batch':
                    with_NTCs = 0  # exclude NTC controls from extraction batch test
                else:
                    with_NTCs = 1
                df_controls_c = df_freqs[self.get_relevant_controls_ids(batch_type, batch_list, with_NTCs)]

                #for colon samples that are coming from only one center, we want to apply an FDR here
                if (batch_type == 'Sequencing library #') & ((c == ('colon', 'tumor')) | (c == ('colon', 'normal'))):
                    dummy, pvals_tmp = self.binomial(df_freqs_c, df_controls_c, epsilon_p, apply_fdr=1)
                else:
                    pvals_tmp = self.binomial(df_freqs_c, df_controls_c, epsilon_p, apply_fdr=0)

                df_pvals_c = pd.concat([df_pvals_c, pvals_tmp], axis=1)  ###Check that the concat is OK here

            df_pvals_c.columns = self.batch_types
            #if condition is colon, use the FDR cutoff in the criteria
            if ((c == ('colon', 'tumor')) | (c == ('colon', 'normal'))):
                pass_df = df_pvals_c.apply(
                    lambda x: 1 if ((x[0] <= pval_cutoff) & (x[1] <= pval_cutoff) & (x[2] <= fdr_cutoff)) else 0,
                    axis=1)
            else:
                pass_df = df_pvals_c.apply(lambda x: 1 if ((x[0]<=pval_cutoff) & (x[1]<=pval_cutoff) & (x[2]<=pval_cutoff)) else 0, axis=1)

            #df_pvals_c should contain all taxons / indeces at this point so it is safe to concatenate
            #add on the condition name to the column name
            new_columns = [None] * len(df_pvals_c.columns)
            for i, col in enumerate(df_pvals_c.columns):
                #new_columns[i] = c[0]+'_'+c[1]+'_'+col
                new_columns[i] = (c[0],c[1],col)
            df_pvals_c.columns = new_columns
            df_filter_tracker_raw = pd.concat([df_filter_tracker_raw, df_pvals_c], axis=1) #track the pvalues of each filter type
            df_filter_tracker = pd.concat([df_filter_tracker, df_pvals_c.applymap(lambda x: 1 if x<=pval_cutoff else np.nan)], axis=1) #track the pvalues
            if ((c == ('colon', 'tumor')) | (c == ('colon', 'normal'))):
                df_filter_tracker[(c[0], c[1], 'Sequencing library #')] = df_filter_tracker_raw[(c[0], c[1], 'Sequencing library #')]\
                    .map(lambda x: 1 if x<=fdr_cutoff else np.nan)


            # from now in, dealing with batch hits
            df_freqs_c = df_freqs_c.loc[pass_df==1,:]

            # per center binomial
            center_dict = meta.generate_dict('Center', tolower=1)
            df_freqs_c_center_group = df_freqs_c.groupby(by=center_dict, axis=1, sort=False)
            centers = df_freqs_c_center_group.groups.keys()

            if len(centers) == 1:  # only one center, pass_levels remain at 3
                print 'Condition %s has samples from only one center, not applying center filter' % (str(c))
                #need to merge results because we've downsized the tables
                condition_stats = condition_stats.merge(pass_df.to_frame(), left_index=True, right_index=True, how='outer')  # just add batch filter results
            else:  # condition is represented by more than one center
                center_pvals = pd.DataFrame()
                for cent in centers:
                    # make a dataframe of the relevant center data
                    df_freqs_c_center = df_freqs_c[df_freqs_c_center_group.groups[cent]]
                    df_controls_c_center_ids = list()
                    for batch_type in self.batch_types:
                        batch_list = list(set(list(meta.df_meta.loc[batch_type, condition_ids])))
                        df_controls_c_center_ids.extend(
                            self.get_relevant_controls_ids(batch_type, batch_list, 1))  # with NTCs here
                    df_controls_c_center_ids = list(
                        set(df_controls_c_center_ids))  # remove duplicate sample ids from the list
                    #downsize control list to include only indeces from df_freqs_c_center, or df_freqs_c, should be the same row index
                    df_stats_c_center = self.binomial(df_freqs_c_center, df_freqs.loc[df_freqs_c_center.index, df_controls_c_center_ids], epsilon_p,
                                                      apply_fdr=0)

                    center_pvals = pd.concat([center_pvals, df_stats_c_center], axis=1)

                center_pvals.columns = centers
                statistics = center_pvals.apply(lambda x: stats.combine_pvalues(x), axis=1)
                combined_pvals = statistics.apply(lambda x: pd.Series(x[1]))
                (dummy, combined_pvals_fdr) = fdr_correction(combined_pvals, alpha = 0.05, method='indep')
                pass_ids = center_pvals.apply(lambda x: 1 if (x < pval_cutoff).sum() >= how_many_centers else 0, axis=1)
                pass_ids = pass_ids & (combined_pvals_fdr <= fdr_cutoff)[:,0]
                pass_ids = pass_ids.map(lambda x: 1 if x==True else 0) #just convert boolean back to binary
                condition_stats = condition_stats.merge(pass_ids.to_frame(), left_index=True, right_index=True, how='outer')

                new_columns = [None] * len(center_pvals.columns)
                for i, col in enumerate(center_pvals.columns):
                    #new_columns[i] = c[0] + '_' + c[1] + '_' + col
                    new_columns[i] = (c[0],c[1],col)
                center_pvals.columns = new_columns

                df_filter_tracker_raw = df_filter_tracker_raw.merge(center_pvals, left_index=True, right_index=True, how='left')
                df_filter_tracker_raw = df_filter_tracker_raw.merge(pd.Series(combined_pvals_fdr[:,0], index=combined_pvals.index,
                                                                      name=(c[0],c[1],'combined_pval_fdr')).to_frame(),
                                                            left_index=True, right_index=True, how='left')

                df_filter_tracker = df_filter_tracker.merge(center_pvals.applymap(lambda x: 1 if x <= pval_cutoff else np.nan),
                                                            left_index=True, right_index=True,how='left')
                df_filter_tracker = df_filter_tracker.merge(pd.Series(combined_pvals_fdr[:,0], index=combined_pvals.index,
                                                                      name=(c[0],c[1],'combined_pval_fdr'))
                                                        .map(lambda x: 1 if x <= fdr_cutoff else np.nan).to_frame(),
                                                        left_index=True, right_index=True, how='left')
                #add a columns depicting whether the taxons pass the entire center filter or not
                df_filter_tracker = df_filter_tracker.merge(pd.Series(pass_ids.map(lambda x: 1 if x == 1 else np.nan),
                                                                  name=(c[0], c[1], 'center_filter_summary')).to_frame(),
                                                        left_index=True, right_index=True, how='left')

        # be aware the the same bacteria can get more than one level assignment per condition!!
        # condition_stats.set_index('index', drop=True, inplace=True).drop(('index',''), axis=1, inplace=True)
        condition_stats.columns = conditions
        condition_stats.replace(0, np.nan, inplace=True)
        #this stage is necessary since adding colon adds a bunch of junk
        idx = condition_stats.index[~(condition_stats.isnull().all(1))]
        condition_stats = condition_stats.loc[idx]

        # don't change the condition stats table, just merge for printing
        condition_stats.merge(self.dataLoader_object.df_taxonomies, left_index=True, right_index=True, how='inner').to_csv(
            os.path.join(self.dataLoader_object.FILES_PATH, 'confidence_levels_per_condition.csv'))

        condition_stats.to_pickle(os.path.join(self.dataLoader_object.FILES_PATH, 'confidence_levels_per_condition.dat'))

        return condition_stats, df_filter_tracker_raw, df_filter_tracker


    def batch_filters(self, meta, df_freqs, pval_cutoff = 0.05, fdr_cutoff = 0.25):
        conditions = meta.get_conditions()
        conditions = filter(lambda x: (x != ('control', 'control') and (x != ('ntc', 'ntc'))), conditions)

        num_controls_total = len(meta.get_control_ids(with_NTCs=1))
        epsilon_p = 1.0 / num_controls_total

        condition_stats = pd.DataFrame()
        for c in conditions:
            print c
            condition_ids = meta.get_condition_ids(c)
            df_freqs_c = df_freqs[condition_ids]
            df_pvals_c = pd.DataFrame()

            for batch_type in self.batch_types:
                batch_list = list(set(list(meta.df_meta.loc[batch_type, condition_ids])))

                if batch_type == 'DNA Extraction Batch':
                    with_NTCs = 0 #exclude NTC controls from extraction batch test
                else:
                    with_NTCs = 1
                df_controls_c = df_freqs[self.get_relevant_controls_ids(batch_type, batch_list, with_NTCs)]

                if batch_type == 'Sequencing library #':
                    dummy, pvals_tmp = self.binomial(df_freqs_c, df_controls_c, epsilon_p, apply_fdr=1)
                else:
                    pvals_tmp = self.binomial(df_freqs_c, df_controls_c, epsilon_p, apply_fdr=0)

                df_pvals_c = pd.concat([df_pvals_c, pvals_tmp], axis=1) ###Check that the concat is OK here

            df_pvals_c.columns = self.batch_types
            #pass_df = df_pvals_c.apply(lambda x: 1 if ((x[0]<=pval_cutoff) & (x[1]<=pval_cutoff) & (x[2]<=pval_cutoff)) else 0, axis=1)
            pass_df = df_pvals_c.apply(
                lambda x: 1 if ((x[0] <= pval_cutoff) & (x[1] <= pval_cutoff) & (x[2] <= fdr_cutoff)) else 0, axis=1)
            condition_stats = pd.concat([condition_stats, pass_df], axis=1) #don't need merge becuase the indeces of the different conditions are identical

        condition_stats.columns = conditions #do we need this?
        condition_stats.replace(0, np.nan, inplace=True)

        # don't change the condition stats table, just merge for printing
        condition_stats.merge(self.dataLoader_object.df_taxonomies, left_index=True, right_index=True, how='inner').to_csv(
            os.path.join(self.dataLoader_object.FILES_PATH, 'hit_calling_per_condition_batch_filters.csv'))

        return condition_stats

    def center_filters(self, meta, df_freqs, how_many_centers, pval_cutoff=0.05):
        conditions = meta.get_conditions()
        conditions = filter(lambda x: (x != ('control', 'control') and (x != ('ntc', 'ntc'))), conditions)

        num_controls_total = len(meta.get_control_ids(with_NTCs=1))
        epsilon_p = 1.0 / num_controls_total
        #pval_cutoff = 0.05

        condition_stats = pd.DataFrame()
        for c in conditions:
            #print c
            condition_ids = meta.get_condition_ids(c)
            df_freqs_c = df_freqs[condition_ids]

            center_dict = meta.generate_dict('Center', tolower=1)
            df_freqs_c_center_group = df_freqs_c.groupby(by=center_dict, axis=1, sort=False)
            centers = df_freqs_c_center_group.groups.keys()
            #pass_ids = pd.DataFrame(np.tile(0, len(df_freqs_c.index)))
            pass_ids = pd.DataFrame(np.tile(1, len(df_freqs_c.index))) #Deborah wants the one center bugs to just pass to the final hit calling file
            pass_ids.index = df_freqs_c.index
            if len(centers) == 1:  # only one center, pass_levels remain at 3
                print 'Condition %s has samples from only one center, not applying center filter' % (str(c))
                condition_stats = pd.concat([condition_stats, pass_ids], axis=1)
            else:  # condition is represented by more than one center
                center_pvals = pd.DataFrame()
                for cent in centers:
                    # make a dataframe of the relevant center data
                    df_freqs_c_center = df_freqs_c[df_freqs_c_center_group.groups[cent]]
                    df_controls_c_center_ids = list()
                    for batch_type in self.batch_types:
                        batch_list = list(set(list(meta.df_meta.loc[batch_type, condition_ids])))
                        df_controls_c_center_ids.extend(self.get_relevant_controls_ids(batch_type, batch_list, 1)) #with NTCs here
                    df_controls_c_center_ids = list(set(df_controls_c_center_ids)) #remove duplicate sample ids from the list
                    df_stats_c_center = self.binomial(df_freqs_c_center, df_freqs[df_controls_c_center_ids], epsilon_p, apply_fdr=0)

                    center_pvals = pd.concat([center_pvals, df_stats_c_center], axis=1)
                    #center_pvals.index = df_stats_c_center.index
                # update bacterial confidence level to 2 if the bacteria pass in at least two centers AND have already passed the original binomial test
                center_pvals.columns = centers
                if c == ('paraf control', 'pcontrol'):
                    center_pvals.to_csv(os.path.join(self.dataLoader_object.FILES_PATH, 'hit_calling_paraffin_controls_per_center_pvalues.csv'))
                pass_ids = center_pvals.apply(lambda x: 1 if (x<pval_cutoff).sum() >= how_many_centers else 0, axis=1)
                condition_stats = pd.concat([condition_stats, pass_ids], axis=1)

        condition_stats.replace(0, np.nan, inplace=True)
        condition_stats.columns = conditions
        #condition_stats.index = df_freqs_c.index

        # don't change the condition stats table, just merge for printing
        condition_stats.merge(self.dataLoader_object.df_taxonomies, left_index=True, right_index=True,
                              how='inner').to_csv(
            os.path.join(self.dataLoader_object.FILES_PATH, 'hit_calling_per_condition_two_centers.csv'))

        return condition_stats


    def paraffin_filters(self, meta, df_freqs, pval_cutoff = 0.05):
        conditions = meta.get_conditions()
        conditions = filter(lambda x: (x != ('control', 'control') and (x != ('ntc', 'ntc'))), conditions)
        epsilon_p = 1*10**(-295) #let's try this very small p, does every binomial pass this?

        condition_stats = pd.DataFrame()
        for c in conditions:
            #print c
            condition_ffpe_ids = meta.get_condition_ids_by_rowname(c, 'Material', ['ffpe'], tolower=1) #get only FFPE samples
            #if there are no ffpe samples, then everything passes so make pass_id a series of ones
            if len(condition_ffpe_ids) == 0:
                pass_df = pd.DataFrame(np.tile(1, len(df_freqs.index))) #e.g. pancreas tumor
                pass_df.index = df_freqs.index
            else:
                df_freqs_c_ffpe = df_freqs[condition_ffpe_ids]
                # find the centers where the condition's FFPE samples were taken from
                center_dict = meta.generate_dict('Center', tolower=1)
                df_freqs_c_center_group = df_freqs_c_ffpe.groupby(by=center_dict, axis=1, sort=False)
                centers = df_freqs_c_center_group.groups.keys()

                #now forget about ffpe samples specifically and take all the condition's samples
                condition_ids = meta.get_condition_ids(c)
                df_freqs_c = df_freqs[condition_ids]
                df_pvals_c = pd.DataFrame()

                df_controls_c = df_freqs[self.get_relevant_paraffin_ids(centers)]

                df_pvals_c = self.binomial(df_freqs_c, df_controls_c, epsilon_p, apply_fdr=0)

                pass_df = df_pvals_c.apply(lambda x: 1 if x <=pval_cutoff else 0)

            condition_stats = pd.concat([condition_stats, pass_df], axis=1) #don't need merge becuase the indeces of the different conditions are identical

        condition_stats.columns = conditions #do we need this?
        condition_stats.replace(0, np.nan, inplace=True)

        # don't change the condition stats table, just merge for printing
        condition_stats.merge(self.dataLoader_object.df_taxonomies, left_index=True, right_index=True, how='inner').to_csv(
            os.path.join(self.dataLoader_object.FILES_PATH, 'hit_calling_per_condition_paraffin_filters.csv'))

        return condition_stats


    def binomial(self, df_freqs_test, df_freqs_control, epsilon_p, apply_fdr):

        # expected probability for binomial distribution based on control samples
        num_controls = len(df_freqs_control.columns)
        p = df_freqs_control.count(axis=1) / num_controls
        p[p == 0] = epsilon_p
        p.index = df_freqs_control.index

        N = len(df_freqs_test.columns)
        b_test_df = pd.DataFrame(zip(df_freqs_test.count(axis=1), p))

        b_test_df['pvals'] = b_test_df.apply(lambda x: stats.binom_test(x[0], N, x[1], alternative="greater"), axis=1)
        if apply_fdr:
            (dummy, b_test_df['FDR']) = fdr_correction(b_test_df['pvals'].replace(np.nan, 1), alpha=0.05,
                                                       method='indep')
        b_test_df.index = df_freqs_test.index

        if apply_fdr:
            return b_test_df['pvals'], b_test_df['FDR']
        else:
            return b_test_df['pvals']


    def mask_results(self, condition_stats, turn_off_mask):
        zero_values = np.tile(0, condition_stats.shape[1])
        condition_stats.loc[turn_off_mask==1, :] = zero_values
        return condition_stats

    #current assumption is that the index is the same
    def combine_results(self, conditiona_stats1, condition_stats2, how, desc1, desc2):
        if how == 'intersect':
            p = pd.Panel({'first': conditiona_stats1, 'second': condition_stats2})
            s = p.sum(axis=0)
            s = s.applymap(lambda x: 1 if x > 1 else np.nan)
            s.merge(self.dataLoader_object.df_taxonomies, left_index=True, right_index=True,
                                  how='inner').to_csv(
                os.path.join(self.dataLoader_object.FILES_PATH, 'hit_calling_per_condition_'+desc1+'_'+desc2+'.csv'))
            return s
        else:
            print 'This method is not supported'

    #this function should work even when the indeces are NOT the same
    def merge_results(self, condition_stats1, condition_stats2, desc1, desc2):
        #only consider the indeces in the first variable
        condition_stats2_tmp = condition_stats2.loc[condition_stats1.index, :]
        result = self.combine_results(condition_stats1, condition_stats2_tmp, 'intersect', desc1, desc2)
        return result

    def get_relevant_controls_ids(self, batch_type, batch_list, with_NTCs):
        meta = self.dataLoader_object.get_metadata_object()
        control_ids = meta.get_controls_ids_by_batch(batch_type, batch_list, with_NTCs)
        return control_ids

    #get paraffin controls that came from the same group of centers as the condition itself
    def get_relevant_paraffin_ids(self, centers):
        meta = self.dataLoader_object.get_metadata_object()
        ids = meta.get_condition_ids_by_rowname(('paraf control', 'pcontrol'), 'Center', centers, tolower=1)
        return ids

    def prevalence_analysis(self, add_all_condition):
        df_freqs = self.dataLoader_object.df_freqs
        condition_dict = self.dataLoader_object.meta.generate_condition_dict()
        # get condition sizes
        df_sizes = df_freqs.groupby(by=condition_dict, axis=1, sort=False).size()
        df_sizes[('control', 'control')] = df_sizes[('control', 'control')] + df_sizes[('ntc', 'ntc')]
        if add_all_condition:
            df_sizes['all'] = df_sizes.sum() - (df_sizes[('control', 'control')] + df_sizes[('ntc', 'ntc')] + df_sizes[('paraf control', 'pcontrol')])
            # generate table of prevalence per condition (per bacteria), combine control and ntc numbers
            # use the df_freqs table which has NaNs for bacteria that didn't exist or pass threshold in the condition. That way 'count()' will count the prevalence of a bacteria in a condition when grouping by that condition
        df_counts = df_freqs.groupby(by=condition_dict, axis=1, sort=False).count()
        df_counts[('control', 'control')] = df_counts[('control', 'control')] + df_counts[('ntc', 'ntc')]
        if add_all_condition:
            # df_counts.drop(('ntc','ntc'), axis=0, inplace=True)
            df_counts['all'] = df_counts.sum(axis=1) - (df_counts[('control', 'control')] + df_counts[('ntc', 'ntc')] + df_counts[('paraf control', 'pcontrol')])
        df_prevalences = df_counts.divide(df_sizes)

        self.dataLoader_object.df_taxonomies.merge(df_prevalences, left_index=True, right_index=True, how='inner').to_csv(os.path.join(self.dataLoader_object.FILES_PATH, 'df_prevalences_per_condition.csv'))

    def prevalence_analysis_groupby_criteria(self, groupby_criteria, tolower_list):
        df_freqs = self.dataLoader_object.df_freqs
        groupby_dict = self.dataLoader_object.meta.generate_groupby_dict(groupby_criteria, tolower_list)
        # get group sizes
        df_sizes = df_freqs.groupby(by=groupby_dict, axis=1, sort=False).size()
        df_counts = df_freqs.groupby(by=groupby_dict, axis=1, sort=False).count()
        df_prevalences = df_counts.divide(df_sizes)

        self.dataLoader_object.df_taxonomies.merge(df_prevalences, left_index=True, right_index=True,
                                                   how='inner').to_csv(
            os.path.join(self.dataLoader_object.FILES_PATH, 'df_prevalences_'+'_'.join(groupby_criteria)+'.csv'))

        return df_prevalences

class analysisPipeline:
    def main(self, do_species):
        #calls the data loader

        pipeline_data = dataLoader('./results')

        #calls the data loader functions to tweak the data once loaded
        #the order is as follows:
        if do_species:
            pipeline_data.load_data('metadata.xlsx',
                                    'reconstruction.txt',
                                    'reconstruction_first_line',
                                    './',
                                    './results')

            #normalize read counts via library average
            pipeline_data.generate_normalized_read_counts()

            print 'Number of samples per condition (before min reads filter)'
            condition_dict = pipeline_data.meta.generate_condition_dict()
            print pipeline_data.meta.df_meta.groupby(by=condition_dict, axis=1, sort=False).count().loc['tissue type']

            #remove samples without a minimum of 1000 reads
            pipeline_data.remove_min_read_samples(1000)
            #flooring very low abundances - this happens after normalizing reads, so the normalized reads table is updating accordingly
            pipeline_data.floor_abundances(10**(-4))

            pipeline_data.df_taxonomies.merge(pipeline_data.df_reads_normalized, left_index=True, right_index=True, how='inner').to_csv(
                os.path.join(pipeline_data.FILES_PATH, 'df_reads_normalized_before_global_contaminants_filter.csv'))
            pipeline_data.df_taxonomies.merge(pipeline_data.df_freqs, left_index=True, right_index=True, how='inner').to_csv(
                os.path.join(pipeline_data.FILES_PATH, 'df_freqs_before_global_contaminants_filter.csv'))


            # find environmental contaminations (use default percentage of 7.5% for now)
            df_environmental_contaminants, df_non_environmental_contaminants_mask, df_contols_percentages = \
                pipeline_data.filter_environmental_contaminants(0.075)

            # paraffin controls filter
            df_paraf_contaminants, df_non_paraf_contaminants_mask, df_paraf_percentages = \
                pipeline_data.filter_contaminants(('paraf control', 'pcontrol'), 'paraffin', 0.075)

            idx = pipeline_data.df_freqs.index[(df_non_environmental_contaminants_mask & df_non_paraf_contaminants_mask) ==1]
            pipeline_data.select_rows(idx) # keeps rows that are not contaminants from global or paraffin controls

            # expand the tables to contain the higher level taxonomies - there is a function to focus only on species when relevant
            pipeline_data.bubble_up_taxonomic_levels()
            #save analysis (as a temporary stage)
            pipeline_data.save_data('Analysis_pipeline_16S_results_species')

        else:
           pipeline_data.load_data_from_cache('Analysis_pipeline_16S_results_species')

        #calls the hit caller
        hit_caller = hitCaller(pipeline_data)

        condition_stats_new, df_filter_tracker_raw, df_filter_tracker = hit_caller.new_filters(pipeline_data.meta, pipeline_data.df_freqs, pval_cutoff = 0.05, fdr_cutoff = 0.20, how_many_centers=2)
        hit_caller.save_filtered_data(condition_stats_new, 'new_fdr_0_20')
        desc1 = 'new_fdr_0_20'
        # condition_stats_new = hit_caller.load_filtered_data('new_fdr_0_20') #optional: load cached filter results

        hit_caller.save_filtered_data(df_filter_tracker, 'df_filter_tracker_'+desc1)
        hit_caller.save_filtered_data(df_filter_tracker_raw, 'df_filter_tracker_raw_' + desc1)


        condition_stats_paraffins = hit_caller.paraffin_filters(pipeline_data.meta, pipeline_data.df_freqs, pval_cutoff=0.05)
        hit_caller.save_filtered_data(condition_stats_paraffins, 'paraffin_filtered_data')
        desc3 = 'paraffin'
        #condition_stats_paraffins = hit_caller.load_filtered_data('paraffin_filtered_data') #optional: load cached filter results

        conditions_stats_combined_with_paraffin = hit_caller.merge_results(condition_stats_new,
                                                                             condition_stats_paraffins, desc1, desc3)

        new_columns = [None] * len(condition_stats_paraffins.columns)
        for i, col in enumerate(condition_stats_paraffins.columns):
            # new_columns[i] = c[0] + '_' + c[1] + '_' + col
            new_columns[i] = (col[0], col[1], 'paraffin control')
        condition_stats_paraffins.columns = new_columns

        df_filter_tracker = pd.concat([df_filter_tracker, condition_stats_paraffins], axis=1)


        pipeline_data.df_taxonomies.merge(df_filter_tracker, left_index=True, right_index=True, how='inner').to_csv(os.path.join(pipeline_data.FILES_PATH, 'df_filter_tracker.csv'))



        pipeline_data.df_taxonomies.merge(df_filter_tracker_raw, left_index=True, right_index=True, how='inner').to_csv(os.path.join(pipeline_data.FILES_PATH, 'df_filter_tracker_raw_without_paraffin.csv'))


        conditions_stats_combined_with_paraffin.merge(hit_caller.dataLoader_object.df_taxonomies, left_index=True, right_index=True,
                                       how='inner').to_csv(
                                          os.path.join(hit_caller.dataLoader_object.FILES_PATH,
                                            'hit_calling_per_condition_' + desc1 +'_'+desc3+'.csv'))


        hit_caller.save_hits(conditions_stats_combined_with_paraffin)


if __name__ == "__main__":
    analysisPipeline().main(1)