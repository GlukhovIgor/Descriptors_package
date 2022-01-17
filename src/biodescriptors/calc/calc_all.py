import os
import pandas as pd

from biodescriptors import calc


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
    """
    Forms pandas dataframe with all statistics for single file.

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
    prot_hel_dist: bool, default=True
        If true, then the descriptor calculated by :code:'calc_prot_hel_dist' function will be added.
    pairwise_sep_dist: bool, default=True
        If true, then the descriptor calculated by :code:'calc_pairwise_sep_dist' function will be added.
    com_calpha_angles: bool, default=True
        If true, then the descriptor calculated by :code:'calc_COM_Calpha_angles' function will be added.
    len_of_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_len_of_hel' function will be added.
    angles_between_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_angles_between_hel' function will be added.
    com_clamp: bool, default=True
        If true, then the descriptor calculated by :code:'calc_COM_clamp' function will be added.
    charge_clamp_dist: bool, default=True
        If true, then the descriptor calculated by :code:'calc_charge_clamp_dist' function will be added.
    charge_clamp_angles: bool, default=True
        If true, then the descriptor calculated by :code:'calc_charge_clamp_angles' function will be added.
    acc_per_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_acc_per_hel' function will be added.
        Only calculated if 'dssp' parameter is set to True. 
    dssp_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_dssp_hel' function will be added.
        Only calculated if 'dssp' parameter is set to True.
    sse_content: bool, default=True
        If true, then the descriptor calculated by :code:'calc_sse_content' function will be added.
        Only calculated if 'dssp' parameter is set to True.
    dssp_extra: bool, default=True
        If true, then the descriptor calculated by :code:'calc_dssp_hel' "extra" function will be added.
        Only calculated if 'dssp' parameter is set to True.

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
        descriptors_frames_list.append(calc.prot_hel_dist_to_pandas(filename, ref, base_filename))
    if pairwise_sep_dist:
        descriptors_frames_list.append(calc.pairwise_sep_dist_to_pandas(filename, ref, base_filename))
    if com_calpha_angles:
        descriptors_frames_list.append(calc.COM_Calpha_angles_to_pandas(filename, ref, base_filename))
    if len_of_hel:
        descriptors_frames_list.append(calc.len_of_hel_to_pandas(filename, ref, base_filename))
    if angles_between_hel:
        descriptors_frames_list.append(calc.angles_between_hel_to_pandas(filename, ref, base_filename))
    if com_clamp:
        descriptors_frames_list.append(calc.COM_clamp_to_pandas(filename, clamp_resid, base_filename))
    if charge_clamp_dist:
        descriptors_frames_list.append(calc.charge_clamp_dist_to_pandas(filename, clamp_resid, base_filename))
    if charge_clamp_angles:
        descriptors_frames_list.append(calc.charge_clamp_angles_to_pandas(filename, clamp_resid, base_filename))
    if dssp:
        if acc_per_hel:
            descriptors_frames_list.append(calc.acc_per_hel_to_pandas(filename, ref, base_filename))
        if dssp_hel:
            descriptors_frames_list.append(calc.dssp_hel_to_pandas(filename, ref, base_filename))
        if sse_content:
            descriptors_frames_list.append(calc.sse_content_to_pandas(filename, base_filename))
        if dssp_extra:
            descriptors_frames_list.append(calc.dssp_extra_to_pandas(filename, ref, base_filename))

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
    """
    Forms pandas dataframe with all statistics for all files in filedir and saves it to .csv.

    Parameters:
    ----------
    filedir: str
        Path to folder with all .PDB files for which frame will be constructed.
    output_full_path: str
        Path where resulting frame with descriptors will be saved.
    clamp_resid: list of ints
        TODO: describe.
    ref: list of ints
        TODO: describe.
    dssp: bool, default=True
        Determines whether to calculate DSSP-dependent descriptors or not.
    prot_hel_dist: bool, default=True
        If true, then the descriptor calculated by :code:'calc_prot_hel_dist' function will be added.
    pairwise_sep_dist: bool, default=True
        If true, then the descriptor calculated by :code:'calc_pairwise_sep_dist' function will be added.
    com_calpha_angles: bool, default=True
        If true, then the descriptor calculated by :code:'calc_COM_Calpha_angles' function will be added.
    len_of_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_len_of_hel' function will be added.
    angles_between_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_angles_between_hel' function will be added.
    com_clamp: bool, default=True
        If true, then the descriptor calculated by :code:'calc_COM_clamp' function will be added.
    charge_clamp_dist: bool, default=True
        If true, then the descriptor calculated by :code:'calc_charge_clamp_dist' function will be added.
    charge_clamp_angles: bool, default=True
        If true, then the descriptor calculated by :code:'calc_charge_clamp_angles' function will be added.
    acc_per_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_acc_per_hel' function will be added.
        Only calculated if 'dssp' parameter is set to True.
    dssp_hel: bool, default=True
        If true, then the descriptor calculated by :code:'calc_dssp_hel' function will be added.
        Only calculated if 'dssp' parameter is set to True.
    sse_content: bool, default=True
        If true, then the descriptor calculated by :code:'calc_sse_content' function will be added.
        Only calculated if 'dssp' parameter is set to True.
    dssp_extra: bool, default=True
        If true, then the descriptor calculated by :code:'calc_dssp_hel' "extra" function will be added.
        Only calculated if 'dssp' parameter is set to True.

    Returns:
    -------
    pandas.DataFrame with calculated descriptors.

    """
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
        counter += 1
    final_df = final_df.reset_index().drop(columns=['index'])
    final_df.to_csv(output_full_path, index=False)
    return final_df
