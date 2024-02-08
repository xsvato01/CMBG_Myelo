import pandas as pd
import argparse


def main(varlist: str, input_table: str, outname: str, total_samples: int):
    samples_vars = pd.read_csv(
        varlist, names=["chr", "pos", "ref", "alt", "sample"], sep='\t', dtype={"chr": str,
                                                                                "pos": int,
                                                                                "ref": str,
                                                                                "alt": str,
                                                                                "sample": str})
    agg_dict = {
        'chr': 'first',
        'pos': 'first',
        'ref': 'first',
        'alt': 'first',
        'sample': 'count',
    }
    grouped_vars = samples_vars.groupby(
        by=["chr", "pos", "ref", "alt"], as_index=False).agg(agg_dict)
    grouped_vars.rename(columns={'sample': 'occurence'}, inplace=True)
    # grouped_vars['occurence'] = grouped_vars['sample_count'].astype(int)
    grouped_vars['occurence_ratio'] = grouped_vars['occurence'] / total_samples

    # if no sex chromosome present in table, chr casts to int64
    table = pd.read_csv(input_table, dtype={'Chromosome': object})
    merged_df = pd.merge(grouped_vars, table,
                         right_on=['Chromosome', 'Start_Position',
                                   'Reference_Allele', 'Tumor_Seq_Allele2'],
                         left_on=["chr", "pos", "ref", "alt"],
                         how='right')
    merged_df.drop(["chr", "pos", "ref", "alt"], axis=1, inplace=True)
    merged_df.to_csv(outname, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Group samples for same variants")
    parser.add_argument(
        '--varlist', help='Input path for .tsv with samples and their vars to group', required=True)
    parser.add_argument(
        '--table', help='Input path for .csv parsed vcf file to concat sample names', required=True)
    parser.add_argument(
        '--n', help='Total number of samples', required=True)
    parser.add_argument(
        '--outname', help='Output name', required=True)
    args = parser.parse_args()
    main(args.varlist, args.table, args.outname, int(args.n))
    # main('/Volumes/lamb/home/450402/000000-My_Documents/MARECKOVA_BTK/launch/TP53_20240126/joined/CXCR4_1_F.joinedvariants.tsv',
    #      '/Volumes/lamb/home/450402/000000-My_Documents/MARECKOVA_BTK/launch/TP53_20240126/create_full_table/CXCR4_1_F.mutect2.filt.norm.vep.full.csv',
    #      'test', int('13'))
