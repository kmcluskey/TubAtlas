import pandas as pd
from loguru import logger
from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib import pyplot as plt
from IPython.display import display

import os, sys


class GeneSelector(object):

    """
    This is a class used to select tissue/tissues of interest from FlyMet based on high level abundance
    compared to the surrounding tissues

    """
    def __init__(self, AB_FACTOR=10, FPKM_MIN=10, AVG_MULT=1):

        '''
        :param AB_FACTOR: Abundance factor - how much more abundant you want your tissue to be - default is 10
        '''

        self.DATA_FOLDER = os.path.abspath(os.path.join('.', 'data'))
        self.RESULTS_FOLDER = os.path.abspath(os.path.join('.', 'results'))

        self.enrich_mf, self.enrich_m, self.enrich_f, self.enrich_l, self.f_abundance, \
        self.m_abundance, self.l_abundance, self.mf_abundance = self.set_up_dfs()

        #This file has a few FBgN additions that were not in the sept21 version
        self.all_genes_df = pd.read_csv(os.path.join(self.DATA_FOLDER,'FlyAtlas2_Alltissues_Allgenes_sept21_2.csv'),
                                   encoding='unicode_escape').set_index('FlyBaseID')


        self.df_dict = {("F", "enrich"): self.enrich_f,  ("M", "enrich"): self.enrich_m, ("L", "enrich"): self.enrich_l,
                        ("MF", "enrich"): self.enrich_mf, ("F", "abundance"): self.f_abundance, ("M", "abundance"): self.m_abundance,
                        ("L", "abundance"): self.l_abundance, ("MF", "abundance"): self.mf_abundance,}

        # self.tissue_list = tissue_list
        self.AB_FACTOR = AB_FACTOR
        # If you don't want a primary tissue then just pass zero for this
        self.FPKM_MIN = FPKM_MIN
        # When a tissue/set of tissues are not the only tissue with AB_FACTORxgreater than all other tissues this parameter
        # is used to find out which other tissues have that gene > AVG_MULT*Average
        self.AVG_MULT = AVG_MULT

        self.ab_df_list = [self.f_abundance, self.m_abundance, self.l_abundance]
        self.en_df_list = [self.enrich_f, self.enrich_m, self.enrich_l]

        self.suffix_dict = {'enrich_f': '_en_f', 'enrich_m':'_en_m', 'enrich_l': '_en_l', 'f_abundance':'_ab_f',
                            'm_abundance':'_ab_m', 'l_abundance':'_ab_l' }


    def set_up_dfs(self):

        # Data sets from Andrew Gillen - might be good to grab directly from the FlyAtlas DB
        # Enrichment datasets
        enrich_mf = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Enrichment_AvgAdult.csv'), encoding='unicode_escape').set_index('FPKM')
        enrich_f = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Enrichment_Female.csv'), encoding='unicode_escape').set_index('FPKM')
        enrich_m = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Enrichment_Male.csv'), encoding='unicode_escape').set_index('FPKM')
        enrich_l = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Enrichment_Larval.csv'), encoding='unicode_escape').set_index('FPKM')

        # Abundance datasets
        f_abundance = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Abundance_Female.csv'), encoding='unicode_escape').set_index('FPKM')
        m_abundance = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Abundance_Male.csv'), encoding='unicode_escape').set_index('FPKM')
        l_abundance = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Abundance_Larval.csv'), encoding='unicode_escape').set_index('FPKM')
        mf_abundance = pd.read_csv(os.path.join(self.DATA_FOLDER,'FPKM_Abundance_AvgAdult.csv'), encoding='unicode_escape').set_index('FPKM')


        all_df_list = [enrich_mf, enrich_f, enrich_m, enrich_l, f_abundance, m_abundance, l_abundance, mf_abundance ]

        droplist = ['FBgn0000003', 'FBgn0000022', 'FBgn0000137', 'FBgn0002561', 'FBgn0002563',
                    'FBgn0002564', 'FBgn0002940', 'FBgn0003261', 'FBgn0004034', 'FBgn0004167', 'FBgn0004170',
                    'FBgn0005427', 'FBgn0010019', 'FBgn0010355', 'FBgn0011225', 'FBgn0011230', 'FBgn0011286',
                    'FBgn0011822', 'FBgn0013674', 'FBgn0013686', 'FBgn0013708', 'FBgn0024183', 'FBgn0024989',
                    'FBgn0026197', 'FBgn0026566', 'FBgn0027660', 'FBgn0028583', 'FBgn0029176', 'FBgn0029521',
                    'FBgn0031176', 'FBgn0031981', 'FBgn0032240', 'FBgn0032505', 'FBgn0032940', 'FBgn0034618',
                    'FBgn0035627', 'FBgn0037003', 'FBgn0037004', 'FBgn0037005', 'FBgn0037230', 'FBgn0040370',
                    'FBgn0040371', 'FBgn0041605', 'FBgn0051072', 'FBgn0051130', 'FBgn0051407', 'FBgn0051716',
                    'FBgn0052425', 'FBgn0052626', 'FBgn0052816', 'FBgn0053208', 'FBgn0065083', 'FBgn0065104',
                    'FBgn0065107', 'FBgn0067622', 'FBgn0082990', 'FBgn0083028', 'FBgn0085447', 'FBgn0085771',
                    'FBgn0085813', 'FBgn0085819', 'FBgn0086073', 'FBgn0086906', 'FBgn0250731', 'FBgn0250824',
                    'FBgn0261504', 'FBgn0261836', 'FBgn0261885', 'FBgn0262107', 'FBgn0262810', 'FBgn0262811',
                    'FBgn0263257', 'FBgn0263392', 'FBgn0264556', 'FBgn0265714', 'FBgn0265964', 'FBgn0265965',
                    'FBgn0265966', 'FBgn0266568', 'FBgn0267514', 'FBgn0267515', 'FBgn0267516', 'FBgn0267519',
                    'FBgn0267521', 'FBgn0267522', 'FBgn0267702', 'FBgn0267703', 'FBgn0267742', 'FBgn0283442']

        # Rename columns and index
        for df in all_df_list:
            df.rename(columns={'Thoracicoabdominal ganglion': 'VNC', 'Malpighian Tubules': 'Tubules',
                               'Virgin Spermatheca': 'VSperma', 'Mated Spermatheca': 'MSperma'}, inplace=True)
            df.index.name = "FlyBaseID"
            df.drop(droplist, inplace=True)
            df.drop("Whole body", axis=1, inplace=True)

        # Select the columns that we are interested in from the dataframes/csv files
        mf_columns = ['Head', 'Eye', 'Brain', 'VNC', 'Crop', 'Midgut', 'Hindgut', 'Tubules',
                      'Fat body', 'Salivary gland', 'Heart', 'Carcass', 'Rectal pad',
                      'Testis', 'Accessory glands', 'Ovary', 'VSperma', 'MSperma']
        enrich_mf = enrich_mf[mf_columns].copy()

        male_columns = ['Head', 'Eye', 'Brain', 'VNC', 'Crop', 'Midgut', 'Hindgut', 'Tubules',
                        'Fat body', 'Salivary gland', 'Heart', 'Carcass', 'Rectal pad',
                        'Testis', 'Accessory glands']
        enrich_m = enrich_m[male_columns].copy()

        female_columns = ['Head', 'Eye', 'Brain', 'VNC', 'Crop', 'Midgut', 'Hindgut', 'Tubules',
                          'Fat body', 'Salivary gland', 'Heart', 'Carcass', 'Rectal pad', 'Ovary',
                          'VSperma', 'MSperma']
        enrich_f = enrich_f[female_columns].copy()

        larval_columns = ['CNS', 'Midgut', 'Hindgut', 'Tubules', 'Fat body', 'Salivary gland',
                          'Trachea', 'Carcass', 'Garland cells']
        enrich_l = enrich_l[larval_columns].copy()

        # Name the datasets for ease of use
        enrich_mf.name = 'enrich_mf'
        enrich_f.name = 'enrich_f'
        enrich_m.name = 'enrich_m'
        enrich_l.name = 'enrich_l'

        f_abundance.name = 'f_abundance'
        m_abundance.name = 'm_abundance'
        l_abundance.name = 'l_abundance'
        mf_abundance.name = 'mf_abundance'

        return enrich_mf, enrich_m, enrich_f, enrich_l, f_abundance, m_abundance, l_abundance, mf_abundance

    def write_csv(self, df, filename, tissue):

        try:
            write_dir = os.path.join(self.RESULTS_FOLDER, tissue)

            if os.path.exists(write_dir):
                df.to_csv(os.path.join(write_dir, filename))
                logger.info("Wrote out %s" % os.path.join(write_dir, filename))
            else:
                os.makedirs(write_dir)
                df.to_csv(os.path.join(write_dir, filename))
                logger.info("Wrote out %s" % os.path.join(write_dir, filename))

        except Exception as e:
            logger.warning("Writing failed because of %s " % e)



    def get_abundance_list(self, tissues, pop):

        ab_list = list(self.get_ab_en_df(tissues, pop).index)

        return ab_list

    def analyse_tissues(self, tissues, pop_list):

        """
        This method should eventually write out in Table form all of the data that is useful to a
        set of Tissues and Populations
        :param tissues: List of tissues for analysis
        :param pop: All of the populations to be used in the analysis
        :return: :return: subset_dict - dict for all of the venn diagram subsets and genes in those subsets

        Write and print out all of the information associated with those tissues
        """


        pop_dict = self.get_pop_dict(tissues, pop_list)
        pop_sets = list(pop_dict.values())

        labels = tuple(["Enrich"+x for x in pop_list])
        plot_title = str(self.AB_FACTOR)+"x Abundant genes in "+', '.join(map(str, tissues))

        # Write out the venn diagram and put in the appropriate directory
        tissue = '_'.join(tissues)
        fname = tissue + '_venn_diagram.png'

        if len(pop_list)==3:
            venn3(pop_sets, set_labels=(labels))
            plt.title(plot_title)
            try:
                write_dir = os.path.join(self.RESULTS_FOLDER, tissue)

                if os.path.exists(write_dir):
                    plt.savefig(os.path.join(write_dir, fname))
                    logger.info("Saving figure %s" % os.path.join(write_dir, fname))
                else:
                    os.makedirs(write_dir)
                    plt.savefig(os.path.join(write_dir, fname))
                    logger.info("Saving figure %s" % os.path.join(write_dir, fname))
            except Exception as e:
                logger.warning("Writing failed because of %s " % e)

            plt.show()

        else:
            venn2(pop_sets, set_labels=(labels))
            plt.title(plot_title)
            try:
                write_dir = os.path.join(self.RESULTS_FOLDER, tissue)

                if os.path.exists(write_dir):
                    plt.savefig(os.path.join(write_dir, fname))
                    logger.info("Saving figure %s" % os.path.join(write_dir, fname))
                else:
                    os.makedirs(write_dir)
                    plt.savefig(os.path.join(write_dir, fname))
                    logger.info("Saving figure %s" % os.path.join(write_dir, fname))
            except Exception as e:
                logger.warning("Writing failed because of %s " % e)

            plt.show()

        subset_dict = self.calc_venn_subsets(pop_dict)
        self.write_subset_dfs(tissues, subset_dict)

        return subset_dict

    def calc_venn_subsets(self, pop_dict):
        """
        :param pop_dict: population name (M, F, L) and genes in set
        :return: subset_dict subset_name and genes
        """
        subset_dict = {}
        try:
            set_m = pop_dict['M']
            set_f = pop_dict['F']
            set_l = pop_dict['L']

            subset_dict["m_only"] = list(set_m.difference(set_f.union(set_l)))
            subset_dict["f_only"] = list(set_f.difference(set_m.union(set_l)))
            subset_dict["l_only"] = list(set_l.difference(set_m.union(set_f)))

            # All sets
            subset_dict["mfl"] = list(set.intersection(*map(set, [set_m, set_l, set_f])))

            # M&L but not F
            subset_dict["ml_not_f"] = list(set_l.intersection(set_m).difference(set_f))
            # M&F but not L
            subset_dict["mf_not_l"] = list(set_f.intersection(set_m).difference(set_l))
            # F&L but not M
            subset_dict["fl_not_m"] = list(set_f.intersection(set_l).difference(set_m))

        except KeyError as e:
            logger.error("Currently this method only works when comparing all 3 populations, "
                         "please do that until this is rewritten")


        return subset_dict

    def order_genes(self, genes, pop, ab_en, tissue):
        """

        :param genes: list of genes (str) that need ordering
        :param pop: Population you want to order by: M, F, L, MF (average MF)
        :param tissue: the tissue of interest (that you want to wort on)
        :return: genes ordered by abundance or enrichment (list, of gene strings)
        """

        df = self.df_dict[(pop, ab_en)]

        df_genes = df.loc[genes]
        df_sorted = df_genes.sort_values(tissue, ascending=False)

        return list(df_sorted.index)

    def get_matching_indexes(self, df_list):

        di_list = []
        for d in df_list:
            di = d.index
            di_list.append(di)

        match_list = set.intersection(*map(set, di_list))

        logger.info("There are %d matches in these dfs" % len(match_list))

        return match_list

    def get_pop_dict(self, tissues, pop_list):

        pop_dict = {}

        for pop in pop_list:
            pop_ab = self.get_abundance_list(tissues, pop)
            pop_dict[pop]=set(pop_ab)
        return pop_dict


    def get_ab_en_df(self, tissues, pop):

        """
        :param tissue: Tissue of interest (string)
        :param pop: This is the population of interest (M,F or L)
        :return Abundance/enrichment data frame with highly abundant genes for tissue of interest
        """
        # tissues = ["Tubules", "Testis"]
        ab_df = self.df_dict[(pop, "abundance")]
        en_df = self.df_dict[(pop, "enrich")]

        ab_en_df= pd.DataFrame()

        column_names, tissue_abundance = self.add_tissue_abundance(tissues, ab_df)

        # Grab all of the abundant genes
        ab_genes = tissue_abundance[(tissue_abundance[column_names] == True).all(axis=1)]

        if not ab_genes.empty:

            ### Remove FbGN codes where enrichment is NaN
            # Get all of the FBGN codes from our abundance calc and grab the enrichment values
            matching_enrich_df = en_df.loc[list(ab_genes.index)]

            # Remove the genes where the enrichment is Nan and return the list of genes
            # enric_df[enric_df[tissues].notna().all(axis=1)]
            non_na_indexes = matching_enrich_df[(matching_enrich_df[tissues].notna().all(axis=1))].index

            ab_df_tiss = ab_genes.loc[non_na_indexes][tissues]
            en_df_tiss = en_df.loc[non_na_indexes][tissues]

            ab_en_df = en_df_tiss.merge(ab_df_tiss, on="FlyBaseID", suffixes=['_en', '_ab'])

            #Sort by tubules enrichment of the popultion passed
            ab_en_df = ab_en_df.sort_values(by=ab_en_df.columns[0], ascending=False)


            #Add on the gene names and symbols
            name_df = self.get_gene_code_df(list(ab_en_df.index))
            ab_en_df = pd.merge(name_df, ab_en_df, left_index=True, right_index=True)

        logger.info("We are returning {} genes from {} ordered by {} enrichment ordered by {}".format(ab_en_df.shape[0], pop, tissues, tissues[0]))

        #Save this df when the method is run

        tissue_name = '_'.join(tissues)
        fname = tissue_name + pop + "_all.csv"
        self.write_csv(ab_en_df, fname, tissue_name)


        return ab_en_df


    def split_coding_non_coding_df(self, ab_en_df):

        coding = ab_en_df[ab_en_df.Gene.str.startswith("CG")]
        non_coding = ab_en_df[ab_en_df.Gene.str.startswith("CR")]

        logger.info("Returning {} coding and {} non-coding genes".format(coding.shape[0], non_coding.shape[0]))

        return coding, non_coding

    def add_tissue_abundance(self, tissues, ab_df):

        if all(t in list(ab_df.columns) for t in tissues):

            column_names = []
            tissue_abundance = ab_df.copy()
            primary_tissue = tissues[0]

            # Work only with Tubule tissues where the FPKM value >10
            ab_fpmk = ab_df[ab_df[primary_tissue] > self.FPKM_MIN]

            # Drop both sets of tissues from the DF in order to calculate the max without them

            ab_no_tub = ab_fpmk.drop(tissues, axis=1)
            ab_no_tub['max'] = ab_no_tub.max(axis=1)

            for tissue in tissues:
                # Divide the FPMK tissue value of interest by 10
                tissue_abundance[tissue + '_div10'] = ab_df[tissue] / self.AB_FACTOR
                # Check if tissue_div10 >= the FPKM values in all of the other tissue (as long as it's not zero)
                column_name = tissue + '_ab'
                column_names.append(column_name)

                # Select genes with FPMK >10
                tissue_abundance = tissue_abundance.loc[ab_fpmk.index]
                tissue_abundance[column_name] = (tissue_abundance[tissue + '_div10'] >= ab_no_tub['max']) & (
                        tissue_abundance[tissue + '_div10'] != 0.00)

        else:
            sys.exit("You seem to have chosen a tissue/population combination that is not possible, please try again")


        return column_names, tissue_abundance


    def get_gene_code_dict(self, fbgn_codes):
        """
        :param genes: The fbgn for genes
        :return: A list of gene codes that we can use
        """

        fb_genes = {}
        fb_id, fb_sy, fb_name = None, None, None

        for fbgn in fbgn_codes:
            try:
                row = self.all_genes_df[self.all_genes_df.index == fbgn]
                fb_id = row.Annotation[0]
                fb_sy = row.Symbol[0]
                fb_name = row.Name[0]
            except Exception as e:
                logger.warning("Returning None for {} as {} ".format(fbgn, e))

            fb_genes[fbgn] = fb_id, fb_sy, fb_name

        return fb_genes

    def get_gene_code_df(self, fbgn_codes):
        m_dict = self.get_gene_code_dict(fbgn_codes)
        name_df = pd.DataFrame(m_dict).T
        name_df.columns = ["Gene", "Symbol", "Name"]
        return name_df

    def greater_than_average_tissues(self, pop, ab_en, gene, AVG_MULT=2, GEN_FPMK_MIN = 10):
        """
        :param pop: Population: M, F, L
        :param ab_en: enrichment or abundance
        :param gene: the gene of interest
        :param AVG_MULT: How many times over the average are we interested in
        :return: List of tissues that have greater than average abundance/enrichment for that gene
        """

        # GEN_FPMK_MIN # only return tissues if the FPKM > GEN_FPMK_MIN
        df = self.df_dict[(pop, ab_en)]

        if ab_en=="enrich":
            GEN_FPMK_MIN=0


        #Select gene and drop the main tissue of interest
        gene_info = df.loc[gene]
        gene_info['avg'] = gene_info.mean()

        # Get the results that are AVG_MULT*the average but also with FPKM >10
        gt_10_tissues = gene_info[gene_info > (AVG_MULT * gene_info['avg'])] > GEN_FPMK_MIN
        gt_av_tissues = gene_info[gene_info > (AVG_MULT * gene_info['avg'])]

        gt_av = list(gt_av_tissues.index.values)
        gt_10 = list(gt_10_tissues[gt_10_tissues].index.values)

        # print (gene, "\nAll :", gt_av, "\nFPKM >", GEN_FPMK_MIN,":",gt_10, "\n" )

        return gt_av, gt_10


    def get_tissue_ab_en(self, tissues, genes):
        """
        Given a tissues or list of tissues give the abundance and enrichment values for a list of genes
        :param tissues: Tissue or list of tissues as strings
        :return: a DF of the abundance and enrichment values for the tissues, for the given genes
        """
        en_ab_tissue_df = pd.DataFrame()

        for df_ab, df_en in zip(self.ab_df_list, self.en_df_list):

            ab_name = df_ab.name
            en_name = df_en.name

            for tissue in tissues:
                try:
                    tissue_ab = df_ab[tissue]
                    tissue_en = df_en[tissue]
                    name_ab = tissue_ab.name + self.suffix_dict[ab_name]
                    name_en = tissue_en.name + self.suffix_dict[en_name]

                    tissue_en = tissue_en.rename(name_en)
                    tissue_ab = tissue_ab.rename(name_ab)

                    en_ab_tissue_df = pd.concat([en_ab_tissue_df, tissue_ab, tissue_en], axis=1)

                except KeyError:
                    pass
                    logger.warning("No tissue of this name for these DFs {} {} ".format(ab_name, en_name))

        en_ab_tissue_genes = en_ab_tissue_df.loc[genes]
        gene_codes = self.get_gene_code_df(genes)

        return_df = gene_codes.join(en_ab_tissue_genes)


        return return_df

    def abundant_tissues_df(self, query_genes, ab_en, pop_list, AVG_MULT=1, FPKM_MIN=10):
        """
        Method to return the tissues for which the abundance (or enrichment) is > AVG_MULT* Average of
        that found in the other tissues
        :param query_genes: Genes of interest (FbGn only, strings(
        :param ab_en: enrich or abundance (strings)
        :param pop_list: M, L, and or/F (list of string)
        :param AVG_MULT: How many * the average abundance or enrichment we are interested in
        :param MIN_FPKM: How many * the average abundance or enrichment we are interested in

        :return: DF with higher than average abundant tissues for set of genes in a df gene:population: Tissue_list
        """
        logger.info("Getting tissues greater that {avg_mult}* average with MIN FPMK of {fpmk_min}".format(avg_mult=AVG_MULT,fpmk_min=FPKM_MIN))
        name_dict = {"L": "Larval", "F": "Female", "M": "Male"}
        dict_list = []
        pop_ex_list = []

        for pop in pop_list:
            pop_ex_list.append(name_dict[pop])
            gene_ab_tissue_dict = {}
            for gene in query_genes:
                ab_tissues, ab_10 = self.greater_than_average_tissues(pop, ab_en, gene, AVG_MULT, FPKM_MIN)
                gene_ab_tissue_dict[gene] = ab_10
            dict_list.append(gene_ab_tissue_dict)

        abundant_tissues_df = pd.DataFrame(dict_list, index=pop_ex_list)

        return abundant_tissues_df

    # This orderes the given genes, adds on the ab and en for the selected tissues

    def write_subset_dfs(self, tissues, subset_dict):

        """
        :param tissues: The Tissues of interest
        :param subset_dict: The dictionary with the subset (e.g. m_only) and the genes it contains
        :return: Writes out the subset to a Tissue specfic directory
        """

        #These dictionaries might be better off in some sort of DB structure - this will suffice for mow
        # A disctionary that gives us the other populations to analyse given the subsets.
        other_pops_dict = {'m_only': ["L", "F"], 'f_only': ["L", "M"], 'l_only': ["M", "F"], 'mfl': [],
                           'ml_not_f': ["F"], 'mf_not_l': ["L"], 'fl_not_m': ["M"]}
        # A dictionary to tell us how to order the subsets
        ordered_dict = {'m_only': "M", 'f_only': "F", 'l_only': "L", 'mfl': "MF", 'ml_not_f': "M", 'mf_not_l': "M",
                        'fl_not_m': "F"}

        for k in subset_dict.keys():
            tissue = '_'.join(tissues)
            if subset_dict[k]:
                ordered_mfl = self.order_genes(subset_dict[k], ordered_dict[k], "enrich", tissues[0])
                df = self.get_tissue_ab_en(tissues, ordered_mfl)
                ab_tissue = self.abundant_tissues_df(subset_dict[k], "abundance", other_pops_dict[k], self.AVG_MULT, self.FPKM_MIN).T
                new_df = df.join(ab_tissue)
                file_name = k + '_' + tissue + '.csv'
                logger.info("Writing the DF for {}".format(k))
                display (new_df)
                self.write_csv(new_df, file_name, tissue)
            else:
                logger.info("There are no genes in set %s and hence no output table " % k)








