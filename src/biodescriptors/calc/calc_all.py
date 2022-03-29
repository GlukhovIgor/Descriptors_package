import os
import pandas as pd
import multiprocessing as mp
from tqdm import tqdm

from biodescriptors import calc

class DescCalculator:

    def __init__(
        self,
        ref=None,
        clamp_resid=None,
        config=calc.constraints.DESCRIPTOR_CONFIG_PATH
    ):
        """
        Initialize calculator.
        
        Parameters:
        ----------
        clamp_resid: list of ints, default=None
            Charge clamp residues list.
        ref: list of ints, default=None
            List of amino acid numbers pairs (start, end) for each helix.
        """
        self.ref = ref
        self.clamp_resid = clamp_resid
        self.config = calc.utils.load_config(config) if isinstance(config, str) else config

    def calc_single_file(
        self,
        filename,
    ):
        """
        Forms pandas dataframe with all statistics for single file.

        Parameters:
        ----------
        filename: str
            Name of .pdb file for which descriptors will be calculated.
        Returns:
        -------
        pandas.DataFrame with calculated descriptors.

        """
        cols = ['prot_name']
        base_filename = os.path.basename(filename)
        data = [base_filename]
        descriptors_frames_list = []
        df = pd.DataFrame(pd.Series(data, index=cols[0:len(data)]), columns=cols)

        descriptor_name_list = self.config['descriptors']

        for descriptor in descriptor_name_list:
            descriptors_frames_list.append(
                NAME_TO_FUNC_MAPPING[descriptor](
                    pdb_file=filename,
                    ref=self.ref,
                    clamp_resid=self.clamp_resid,
                    protein_name=base_filename
                )
            )

        for descriptors_frame in descriptors_frames_list:
            df = df.merge(descriptors_frame, on='prot_name')
            
        return df


    def calc_all(
        self,
        filedir,
        save_to_csv=True,
        output_full_path=None,
    ):
        """
        Forms pandas dataframe with all statistics for all files in filedir and saves it to .csv.

        Parameters:
        ----------
        filedir: str
            Path to folder with all .PDB files for which frame will be constructed.
        save_to_csv: bool, default=True
            If true, then resulting dataframe is saved to .csv file based on output_full_path parameter.
        output_full_path: str, default
            Path where resulting frame with descriptors will be saved.
        Returns:
        -------
        pandas.DataFrame with calculated descriptors.

        """
        filenames = os.listdir(filedir)
        full_filenames = [os.path.join(filedir, filename) for filename in filenames if '.pdb' in filename]
        number_files = len(full_filenames)

        results_dfs = []
        with mp.Pool() as pool:
            for result in tqdm(pool.imap(self.calc_single_file, full_filenames), total=number_files):
                results_dfs.append(result)
        
        final_df = pd.concat(results_dfs, ignore_index=True)
        
        if save_to_csv:
            final_df.to_csv(output_full_path, index=False)
        return final_df


#____________NAME_TO_DESCRIPTOR_MAPPING________________#


NAME_TO_FUNC_MAPPING = {
    'prot_hel_dist': calc.prot_hel_dist_to_pandas,
    'pairwise_sep_dist': calc.pairwise_sep_dist_to_pandas,
    'com_calpha_angles': calc.COM_Calpha_angles_to_pandas,
    'len_of_hel': calc.len_of_hel_to_pandas,
    'angles_between_hel': calc.angles_between_hel_to_pandas,
    'com_clamp': calc.COM_clamp_to_pandas,
    'charge_clamp_dist': calc.charge_clamp_dist_to_pandas,
    'charge_clamp_angles': calc.charge_clamp_angles_to_pandas,
    'acc_per_hel': calc.acc_per_hel_to_pandas,
    'dssp_hel': calc.dssp_hel_to_pandas,
    'sse_content': calc.sse_content_to_pandas,
    'dssp_extra': calc.dssp_extra_to_pandas,
}