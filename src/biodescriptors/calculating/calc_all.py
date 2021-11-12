import os
import pandas as pd

from biodescriptors import calculating


def calc_single_file(filedir, filename, clamp_resid, ref):
    """ Forms pandas dataframe with all statistics for single file."""
    cols = ['prot_name']
    data = [filename]
    file_full_name = os.path.join(filedir, filename)
    df = pd.DataFrame(pd.Series(data, index=cols[0:len(data)]), columns=cols)
    df_prothel = calculating.prot_hel_dist_to_pandas(file_full_name, ref, filename)
    df_pairseps = calculating.pairwise_sep_dist_to_pandas(file_full_name, ref, filename)
    df_alphaangle = calculating.COM_Calpha_angles_to_pandas(file_full_name, ref, filename)
    df_sse = calculating.sse_content_to_pandas(file_full_name, filename)
    df_len = calculating.len_of_hel_to_pandas(file_full_name, ref, filename)
    df_cos = calculating.angles_between_hel_to_pandas(file_full_name, ref, filename)
    df_dssp = calculating.dssp_hel_to_pandas(file_full_name, ref, filename)
    df_extra = calculating.dssp_extra_to_pandas(file_full_name, ref, filename)
    df_clamps = calculating.COM_clamp_to_pandas(file_full_name, clamp_resid, filename)
    df_cl_dist = calculating.charge_clamp_dist_to_pandas(file_full_name, clamp_resid, filename)
    df_cl_angles = calculating.charge_clamp_angles_to_pandas(file_full_name, clamp_resid, filename)
    df_acc = calculating.acc_per_hel_to_pandas(file_full_name, ref, filename)
    df_concat = (
        df
        .merge(df_prothel, on='prot_name')
        .merge(df_pairseps,on='prot_name')
        .merge(df_alphaangle,on='prot_name')
        .merge(df_sse, on='prot_name')
        .merge(df_len, on='prot_name')
        .merge(df_cos, on='prot_name')
        .merge(df_dssp, on='prot_name')
        .merge(df_extra, on='prot_name')
        .merge(df_clamps,on='prot_name')
        .merge(df_cl_dist, on='prot_name')
        .merge(df_cl_angles, on='prot_name')
        .merge(df_acc, on='prot_name')
    )
    return df_concat


def calc_all(filedir, output_full_path, clamp_resid, ref):
    """Forms pandas dataframe with all statistics for all files in filedir and saves it to .csv ."""
    filenames = os.listdir(filedir)
    number_files = len(filenames)
    counter = 1
    final_df = pd.DataFrame()
    for filename in filenames:
        print((f'calculating {counter} out of {number_files}, structure - {filename}'))
        final_df = final_df.append(calc_single_file(filedir, filename, clamp_resid, ref))
        counter+=1
    final_df = final_df.reset_index().drop(columns=['index'])
    final_df.to_csv(output_full_path, index=False)
    return final_df
