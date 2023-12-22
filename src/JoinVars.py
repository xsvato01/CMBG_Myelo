import pandas as pd
import argparse


def main(varlist, input_table, outname):
    samples_vars = pd.read_csv(
        varlist, names=["chr", "pos", "ref", "alt", "sample"], sep='\t')
    agg_dict = {'chr': 'first', 'pos': 'first',
                'ref': 'first', 'alt': 'first', 'sample': ', '.join, }
    grouped_vars = samples_vars.groupby(
        by=["chr", "pos", "ref", "alt"], as_index=False).agg(agg_dict)
    table = pd.read_csv(input_table)

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
        '--outname', help='Output name', required=True)
    args = parser.parse_args()
    main(args.varlist, args.table, args.outname)
