
class Taxonomy:

    def __init__(self):
        self.taxon_levels_list = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
        #can convert this to a Dict on zip()
        self.offset_levels = {'domain':80000, 'phylum':70000, 'class':60000, 'order':50000, 'family':40000, 'genus':30000, 'species':0}
        self.taxon_levels = {0:'domain', 1:'phylum', 2:'class', 3:'order', 4:'family', 5:'genus', 6:'species'}
        self.taxon_levels_inv = {'domain':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species':6}

    def filter_ids(self, input_ids, taxon_level):
        if taxon_level not in self.taxon_levels_list:
            print 'Error in taxon level, not filtering'
            return input_ids

        next_taxon_level = self.taxon_levels_inv[taxon_level] - 1
        if next_taxon_level < 0:
            ids = input_ids[input_ids >= self.offset_levels[taxon_level]]
        else:
            ids = input_ids[
                (input_ids >= self.offset_levels[taxon_level]) & (input_ids < self.offset_levels[self.taxon_levels[next_taxon_level]])]
        return ids

    def get_taxon_ids(self, df_taxonomies, taxon_level, taxon_name):
        return df_taxonomies.index[df_taxonomies[taxon_level] == taxon_name]
