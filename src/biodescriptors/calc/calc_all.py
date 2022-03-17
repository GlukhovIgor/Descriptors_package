import os
import pandas as pd

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
                calc.utils.NAME_TO_FUNC_MAPPING[descriptor](
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
        number_files = len(filenames)
        counter = 1
        final_df = pd.DataFrame()
        for filename in filenames:
            print((f'calculating {counter} out of {number_files}, structure - {filename}'))
            file_full_name = os.path.join(filedir, filename)
            final_df = final_df.append(self.calc_single_file(file_full_name))
            counter += 1
        final_df = final_df.reset_index().drop(columns=['index'])
        if save_to_csv:
            final_df.to_csv(output_full_path, index=False)
        return final_df
