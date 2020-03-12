import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class Metadata:
    #properties should be:
        #the parsed file (df_meta)
        #the transposed file (df_meta_t)

    # constructor reads in metadata file and parses it to create the datatable with
    # relevant data for tumor microbiome samples
    def __init__(self, metadata_filename):
        #HEADS UP if you change this, you also need to change the 'condition' variable in Analysis_pipeline_16S.py if you are going to run the pipeline
        self.condition_line = 'Tumor_Normal_Control_NTC_Unknown_Pcontrol'
        #self.condition_line = 'Tumor_Pre-Malignant_NAT_Normal_Control_NTC_Unknown_Pcontrol'

        # Get metadata file
        filename_meta = metadata_filename

        print "Reading meta data from: " + filename_meta
        meta = pd.read_excel(filename_meta, "Sheet1", header=1, skiprows=0, index_col=0)
        # fill empty spaces with NaNs so they can easily be dropped
        meta.replace('', np.nan, inplace=True)
        # drop empty rows (0)
        meta.dropna(axis=0, how='all', inplace=True)

        ##### not needed for all metadata '''
        # get rid of columns that are irrelevant to the experiment (controls...)
        # this is defined in the metadata file as 'Not relevant to cancer microbiome'
        ix = meta.index.get_loc("Group # (cancer microbiome project (Deborah))")
        # find the columns for which there is positive data that designates the sample to a different project
        # ix_bool is a Series so we can use its values aen(meta.columnss criteria for the 'meta' dataframe
        ix_bool = meta.iloc[ix, :].notnull()
        relevant_cols = meta.columns[ix_bool.values]
        meta = meta[relevant_cols]

        #change condition names to lower case for accurate grouping
        meta.iloc[meta.index.get_loc(self.condition_line)] = map(lambda x:x.lower(), meta.iloc[meta.index.get_loc(self.condition_line)])
        meta.iloc[meta.index.get_loc('Tumor_Pre-Malignant_NAT_Normal_Control_NTC_Unknown_Pcontrol')] = \
            map(lambda x: x.lower(), meta.iloc[meta.index.get_loc('Tumor_Pre-Malignant_NAT_Normal_Control_NTC_Unknown_Pcontrol')])
        meta.iloc[meta.index.get_loc('tissue type')] = map(lambda x: x.lower(), meta.iloc[
            meta.index.get_loc('tissue type')])

        #fix the condition_line property for breast samples to differentiate between true-normal and nat
        #at this point everything is lower case
        #there is one condition called normal-fa, which we want to consider as normal
        tissue_line_ix = meta.index.get_loc('tissue type')
        condition_line_ix = meta.index.get_loc('Tumor_Normal_Control_NTC_Unknown_Pcontrol')
        condition_line_2_ix = meta.index.get_loc('Tumor_Pre-Malignant_NAT_Normal_Control_NTC_Unknown_Pcontrol')

        breast_tissue_samples = (meta.iloc[tissue_line_ix] == 'breast')
        breast_tissue_samples = breast_tissue_samples.index[breast_tissue_samples == True]
        breast_condition_types_dict = {'normal': 'normal', 'tumor': 'tumor', 'nat': 'nat', 'normal-fa': 'normal',
                                       'pre-malignant': 'pre-malignant'}
        meta.iloc[condition_line_ix][breast_tissue_samples] = meta.iloc[condition_line_2_ix][breast_tissue_samples].map(breast_condition_types_dict)

        self.df_meta = meta
        # transpose to put the data types in the columns and then
        df_meta_t = meta.transpose()
        # reset index to make the sampled IDs a field
        df_meta_t.reset_index(inplace=True, col_fill='SampleID') #this doesn't rename the column to 'Sample ID'

        #df_meta_t['Tumor_Normal_Control_NTC_Unknown_Pcontrol'] = df_meta_t['Tumor_Normal_Control_NTC_Unknown_Pcontrol'].apply(
         #   lambda x: str(x).lower())
        #df_meta_t['Tumor_Pre-Malignant_NAT_Normal_Control_NTC_Unknown_Pcontrol'] = df_meta_t['Tumor_Pre-Malignant_NAT_Normal_Control_NTC_Unknown_Pcontrol'].apply(lambda x: str(x).lower())

        self.df_meta_t =  df_meta_t

    def add_row(self, rowdata):
        self.df_meta = self.df_meta.append(rowdata)
        self.df_meta_t = self.df_meta.transpose()
        self.df_meta_t.reset_index(inplace=True, col_fill='SampleID')

    def generate_condition_dict(self):
        sample_names = self.df_meta.columns

        condition_ix = self.df_meta.index.get_loc(self.condition_line)
        conditions = self.df_meta.iloc[condition_ix, :].str.lower()  # ensure that there aren't differences in the cases

        tissue_str_ix = self.df_meta.index.get_loc("tissue type")
        tissues_str = self.df_meta.iloc[tissue_str_ix, :].str.lower()

        return dict(zip(sample_names, zip(tissues_str, conditions)))

    def generate_groupby_dict(self, rownames, tolower_list):
        sample_names = self.df_meta.columns
        columns = [None] * len(rownames)
        for i, r in enumerate(rownames):
            row_ix = self.df_meta.index.get_loc(r)
            if tolower_list[i] == 1: #string
                columns[i] = self.df_meta.iloc[row_ix, :].str.lower()
                #columns[i] = columns[i].str.strip()
            else: #number
                columns[i] = pd.to_numeric(self.df_meta.iloc[row_ix, :])

        return dict(zip(sample_names, zip(*columns)))

    def generate_dichotomous_meta_in_tissue_dict(self, rowname, tolower):
        sample_names = self.df_meta.columns

        dichotomous_ix = self.df_meta.index.get_loc(rowname)
        if tolower:
            #conditions = self.df_meta.iloc[dichotomous_ix, :].str.lower()  # ensure that there aren't differences in the cases
            conditions = self.df_meta.iloc[dichotomous_ix, :].astype(str).str.lower()
            conditions = conditions.str.strip()
        else:
            conditions = self.df_meta.iloc[dichotomous_ix, :]

        tissue_str_ix = self.df_meta.index.get_loc("tissue type")
        tissues_str = self.df_meta.iloc[tissue_str_ix, :].str.lower()

        return dict(zip(sample_names, zip(tissues_str, conditions)))


    #if row is strings, they should be converted to lower case to prevent redundancies
    def generate_dict(self, row_name, tolower):
        # make list of sample names
        sample_names = self.df_meta.columns
        # make list of conditions
        row_ix = self.df_meta.index.get_loc(row_name)
        if (tolower):
            row_data = self.df_meta.iloc[row_ix, :].astype(str).str.lower() # ensure that there aren't differences in the cases
            #row_data.apply(lambda x: x.replace('nan', np.nan))
        else:
            row_data = pd.to_numeric(self.df_meta.iloc[row_ix, :])  # these are integers so no need to deal with upper/lowercase
        # zip together into dictionary
        condition_dict = dict(zip(sample_names, row_data))

        return condition_dict


    def plot_meta(self):
        #use the meta for queries table
        #plot number of samples per tissue and material type
        df_meta_t_group = self.df_meta_t.groupby(['tissue type'])['Material'].value_counts()
        df = df_meta_t_group.unstack(1) #makes the different materials into different series/colors (legend)
        df = df.loc[:, df.columns.notnull()] #get rid of NaN values for plot
        df.plot(kind='bar', title='Sample material types per tissue')
        plt.tight_layout()

        # plot number of samples per condition (tissue/type) and material type
        df_meta_t_group = self.df_meta_t.groupby(['tissue type', 'Tumor_Normal_Control_NTC_Unknown_Pcontrol'])[
            'Material'].value_counts()
        df = df_meta_t_group.unstack(2)
        df = df.loc[:, df.columns.notnull()]
        df.plot(kind='bar', title='Sample material types per condition')
        plt.tight_layout()

        #more statistics/numbers?
        plt.show()

    #a function to extract the sampleIDs that are relevant to a particular grouping?
    #returns a Series with sampleIDs as values

    def draw_frequencies(self, ax, data, label, color):
        d = ~np.isnan(data.astype('float64'))  # get rid of nan values to make the histogram work
        #ax.hist(np.log10(d), bins=100, normed=True, cumulative=True, histtype='step', label=label, color=color)
        ax.hist(d, bins=10, normed=True, label=label, color=color)

    # metadata is the name of the row of metadata you want to plot by condition
    # type is the the type of metadata,
    #   if it's categorical, the plot will be a bar plot for each category
    #   if it's continuous, the plot will be a boxplot for each condition
    def plot_meta_by_condition(self, metadata, type, desc):
        if type == 'categorical':
            df_meta_t_group = self.df_meta_t.groupby(['tissue type', 'Tumor_Normal_Control_NTC_Unknown_Pcontrol'])[
                metadata].value_counts()
            df = df_meta_t_group.unstack(2)
            df = df.loc[:, df.columns.notnull()]
            df.plot(kind='bar', title=desc+' per condition')

        else: # type=continuous - works for some columns, e.g. material, not clear why not for others....
            local_df_meta_t = self.df_meta_t.reset_index() #get index into first column
            df_meta_t_group = self.df_meta_t.groupby(['tissue type', 'Tumor_Normal_Control_NTC_Unknown_Pcontrol'])
            fig, axarr = plt.subplots(3, 6, sharey=True, sharex=True)
            i = 0
            j = 0

            for g in df_meta_t_group.groups.keys():

                #print 'Group is ' + str(g).replace("u'", "")
                data = local_df_meta_t.loc[df_meta_t_group.groups[g],metadata]
                d = data[~np.isnan(data.astype('float64'))]  # get rid of nan values to make the histogram work

                # we have samples from each type in this condition
                if len(data) > 0:
                    axarr[i,j].hist(d, bins=20, normed=False, label=desc, color='purple') #normed=True

                axarr[i, j].set_title(str(g).replace("u'", ""))

                if j == 0:
                    j = j + 1
                elif j % 5 == 0:
                    j = 0
                    i = i + 1
                else:
                    j = j + 1

            plt.legend(loc='best', shadow=True, bbox_to_anchor=(1, 0.5))
        plt.xlabel('Condition')
        plt.tight_layout()

    def get_sample_ids(self, groupby, group, column, condition):
        df_meta_t_group = self.df_meta_t.groupby(groupby)
        return self.df_meta_t.iloc[df_meta_t_group.groups[group]][self.df_meta_t[column] == condition]['index'] #change index to SampleID if I ever get the column name to change

    #reduce meta data tables to include only relevant columns
    def reduce(self, sample_ids):
        self.df_meta = self.df_meta[sample_ids]
        self.df_meta_t = self.df_meta.transpose().reset_index()

    ##################################
    def get_conditions(self):
        return self.df_meta_t.groupby(['tissue type', self.condition_line]).groups.keys()
        #return [('breast', 'tumor')] # try using only one condition


    def get_condition_ids(self, condition):
       return self.df_meta_t.iloc[self.df_meta_t.groupby(['tissue type', self.condition_line]).groups[condition]]['index']

    # This function could move to the metadata class
    def get_tissue_columns(self, tissue):
        samples = []
        tissue_dict = self.generate_dict('tissue type', tolower=1)
        for key, val in tissue_dict.iteritems():
            if val == tissue:
                samples.append(key)
        return samples

    def get_controls_ids_by_batch(self, batch_type, batch_list, with_NTCs):
        control_ids = self.get_control_ids(with_NTCs)

        controls_meta = self.df_meta[control_ids]
        controls_meta = controls_meta.loc[:, map(lambda x: x in batch_list, controls_meta.loc[batch_type,:])]
        return controls_meta.columns

    def get_condition_ids_by_rowname(self, condition, rowname, row_values, tolower):
        condition_ids = self.get_condition_ids(condition)
        # start with relevant condition
        condition_meta = self.df_meta[condition_ids]
        row_data = condition_meta.loc[rowname, :]
        if tolower:
            row_data = map(lambda x: str(x).lower(), row_data)
        # within the condition, locate the samples that have values in the rowname row that fall inside row_values
        condition_meta = condition_meta.loc[:, map(lambda x: x in row_values, row_data)]

        return condition_meta.columns

    def get_control_ids(self, with_NTCs):
        control_ids = self.get_condition_ids(('control', 'control'))
        # XXXXXXXXXXXXXXXXXX
        if with_NTCs:
            control_ids = control_ids.append(self.get_condition_ids(('ntc', 'ntc')))

        return control_ids

    # given two lists of sample (names)
    # we want to generate a nice printout of their
    # average age +- confidence level
    # %Male/%Female
    # %FFPE/%Snap
    # %- in each center
    # average number of reads in the samples +- confidence interval
    # Tumor size?
    # BMI - if we can get it
    # other Batch-like things
    #   #sequencing library numbers
    #   #DNA extraction batch numbers
    def generate_metadata_statistics(self, group1, desc):
        #subset the metadata to the relevant group of samples
        group1_meta = self.df_meta[group1]
        group1_meta_t = group1_meta.transpose()
        group1_meta_t.reset_index(inplace=True, col_fill='SampleID')

        group1_age_avg = group1_meta_t['Age'].mean()
        group1_age_ci_high = group1_age_avg + group1_meta_t['Age'].sem()*1.96
        group1_age_ci_low =  group1_age_avg - group1_meta_t['Age'].sem()*1.96

        group1_sequencing_libraries = ','.join(str(l) for l in set(group1_meta_t['Sequencing library #']))

        group1_centers = ','.join(str(l) for l in set(group1_meta_t['Center']))

        group1_gender = ((group1_meta_t['Gender (F/M)'] == 'M').sum())*100.0 / len(group1)

        group1_material = ((group1_meta_t['Material'] == 'FFPE').sum()) * 100.0 / len(group1)

        print 'Age: %f [%f, %f]\nSequencing Libraries: (%s)\nGender(M): %f.1%%\nMaterial(FFPE): %f.1%%\n' %(group1_age_avg, group1_age_ci_low, group1_age_ci_high, group1_sequencing_libraries, group1_gender, group1_material)

        s = pd.Series(index = ['Age', 'Age_CI_high', 'Age_CI_low', 'Sequencing Libraries', 'Centers', 'Gender(%M)', 'Material(%FFPE)'], data=(group1_age_avg, group1_age_ci_low, group1_age_ci_high, group1_sequencing_libraries, group1_centers, group1_gender, group1_material), name=desc)

        return s


    def incoporate_new_line(self, new_row_content):
        df_meta_t_copy = self.df_meta_t.copy()  # make sure not to stomp on the bona-fide metadata
        df_meta_t_copy.set_index('index', drop=True, inplace=True)  # go back to sample names as index
        df_meta_t_copy = df_meta_t_copy.merge(new_row_content, left_index=True, right_index=True, how='left')
        self.df_meta = df_meta_t_copy.transpose()
        self.df_meta_t = self.df_meta.transpose()
        self.df_meta_t.reset_index(inplace=True, col_fill='SampleID')

    def update_line_with_dict(self, rowname, my_dict):
        df_meta_t_copy = self.df_meta_t.copy()  # make sure not to stomp on the bona-fide metadata
        df_meta_t_copy[rowname] = df_meta_t_copy[rowname].map(my_dict)
        df_meta_t_copy.set_index('index', drop=True, inplace=True)  # go back to sample names as index
        self.df_meta = df_meta_t_copy.transpose()
        self.df_meta_t = self.df_meta.transpose()
        self.df_meta_t.reset_index(inplace=True, col_fill='SampleID')