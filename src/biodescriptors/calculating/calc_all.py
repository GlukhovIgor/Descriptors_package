import os
import pandas as pd

from biodescriptors import calculating


def calc_single_file(filename, 
                    clamp_resid, 
                    ref, 
                    dssp=True,
                    prot_hel_dist=True,
                    pairwise_sep_dist=True,
                    com_calpha_angles=True,
                    len_of_hel=True,
                    angles_between_hel=True,
                    com_clamp=True,
                    charge_clamp_dist=True,
                    charge_clamp_angles=True,
                    acc_per_hel=True,
                    dssp_hel=True,
                    sse_content=True,
                    dssp_extra=True):
    """ Forms pandas dataframe with all statistics for single file.
    
    Parameters:
    ----------
    filename: str
        Name of .pdb file for which descriptors will be calculated. 
    clamp_resid: list of ints
        TODO: describe.
    ref: list of ints
        TODO: describe.
    dssp: bool, default=True
        Determines whether to calculate DSSP-dependent descriptors or not.
        
    Returns:
    -------
    pandas.DataFrame with calculated descriptors.

    """
    cols = ['prot_name']
    base_filename = os.path.basename(filename)
    data = [base_filename]
    descriptors_frames_list = []
    df = pd.DataFrame(pd.Series(data, index=cols[0:len(data)]), columns=cols)
    if prot_hel_dist:
        descriptors_frames_list.append(calculating.prot_hel_dist_to_pandas(filename, ref, base_filename))
    if pairwise_sep_dist:
        descriptors_frames_list.append(calculating.pairwise_sep_dist_to_pandas(filename, ref, base_filename))
    if com_calpha_angles:
        descriptors_frames_list.append(calculating.COM_Calpha_angles_to_pandas(filename, ref, base_filename))
    if len_of_hel:
        descriptors_frames_list.append(calculating.len_of_hel_to_pandas(filename, ref, base_filename))
    if angles_between_hel:
        descriptors_frames_list.append(calculating.angles_between_hel_to_pandas(filename, ref, base_filename))
    if com_clamp:
        descriptors_frames_list.append(calculating.COM_clamp_to_pandas(filename, clamp_resid, base_filename))
    if charge_clamp_dist:
        descriptors_frames_list.append(calculating.charge_clamp_dist_to_pandas(filename, clamp_resid, base_filename))
    if charge_clamp_angles:
        descriptors_frames_list.append(calculating.charge_clamp_angles_to_pandas(filename, clamp_resid, base_filename))
    if dssp:
        if acc_per_hel:
            descriptors_frames_list.append(calculating.acc_per_hel_to_pandas(filename, ref, base_filename))
        if dssp_hel:
            descriptors_frames_list.append(calculating.dssp_hel_to_pandas(filename, ref, base_filename))
        if sse_content:
            descriptors_frames_list.append(calculating.sse_content_to_pandas(filename, base_filename))
        if dssp_extra:
            descriptors_frames_list.append(calculating.dssp_extra_to_pandas(filename, ref, base_filename))

    for descriptors_frame in descriptors_frames_list:
        df = df.merge(descriptors_frame, on='prot_name')
        
    return df


def calc_all(filedir, 
            output_full_path, 
            clamp_resid, 
            ref,
            dssp=True,
            prot_hel_dist=True,
            pairwise_sep_dist=True,
            com_calpha_angles=True,
            len_of_hel=True,
            angles_between_hel=True,
            com_clamp=True,
            charge_clamp_dist=True,
            charge_clamp_angles=True,
            acc_per_hel=True,
            dssp_hel=True,
            sse_content=True,
            dssp_extra=True):
    """Forms pandas dataframe with all statistics for all files in filedir and saves it to .csv ."""
    filenames = os.listdir(filedir)
    number_files = len(filenames)
    counter = 1
    final_df = pd.DataFrame()
    for filename in filenames:
        print((f'calculating {counter} out of {number_files}, structure - {filename}'))
        file_full_name = os.path.join(filedir, filename)
        final_df = final_df.append(calc_single_file(file_full_name, 
                                                    clamp_resid, 
                                                    ref,
                                                    dssp=dssp,
                                                    prot_hel_dist=prot_hel_dist,
                                                    pairwise_sep_dist=pairwise_sep_dist,
                                                    com_calpha_angles=com_calpha_angles,
                                                    len_of_hel=len_of_hel,
                                                    angles_between_hel=angles_between_hel,
                                                    com_clamp=com_clamp,
                                                    charge_clamp_dist=charge_clamp_dist,
                                                    charge_clamp_angles=charge_clamp_angles,
                                                    acc_per_hel=acc_per_hel,
                                                    dssp_hel=dssp_hel,
                                                    sse_content=sse_content,
                                                    dssp_extra=dssp_extra))
        counter+=1
    final_df = final_df.reset_index().drop(columns=['index'])
    final_df.to_csv(output_full_path, index=False)
    return final_df
