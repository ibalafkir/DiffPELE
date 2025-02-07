import os 
from src.utils.pdb_utils import PdbDf, PdbHandler
from src.utils.logger_factory import LoggerFactory
from src.utils.os_utils import display_full_dataframe
from src.pipeline.interface import InterfaceAnalyzer
from biopandas.pdb import PandasPdb
import textwrap


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


class RFdiffContigs:
    
    def __init__(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float
    ):
        
        """
        Initialize the RFdiffContigs class.
        This class pretends to obtain contigs from a system
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        receptor_chains (str): Chain ID for the receptor protein chains.
                               Accepts one chain: 'A'
                               or 2 comma-sepparated: 'A,B'
        ligand_chain (str): Chain ID for the ligand protein chain.
                            Accepts one chain only
        distance_cutOff (float): Distance cutoff for the interface analysis
        """
        
        # Validate inputs
        self._validate_inputs(
            pdb_path,
            receptor_chains,
            ligand_chain,
            distance_cutOff
        )
        
        # Set attributes
        self.pdb_path = pdb_path
        self.distance_cutOff = distance_cutOff
        self.ligand_chain = ligand_chain
        self.receptor_chain_1 = None # will be None if only one chain is provided
        self.receptor_chain_2 = None # will be None if only one chain is provided
        
        if len(receptor_chains) == 1:
            self.receptor_chains = receptor_chains
        if len(receptor_chains) == 3:
            self.receptor_chains = receptor_chains.split(',')
            self.receptor_chain_1 = self.receptor_chains[0]
            self.receptor_chain_2 = self.receptor_chains[1]
    
    
    def _validate_inputs(self,
                         pdb_path, 
                         receptor_chains, 
                         ligand_chain, 
                         distance_cutOff):
        
       # Validate string inputs
        for var_name, var_value in [
            ("pdb_path", pdb_path),
            ("receptor_chains", receptor_chains),
            ("ligand_chain", ligand_chain)
            ]:
            if not isinstance(var_value, str):
                raise TypeError(f"{var_name} should be a string")
        
        # Validate float inputs
        if not isinstance(distance_cutOff, float):
            raise TypeError(f"{distance_cutOff} should be a float")
        
        # Validate PDB
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"File not found: {pdb_path}"
                                    )       
        # Validate distance cut-off
        if distance_cutOff <= 0:
            raise ValueError(f"{distance_cutOff} should be positive")
        
        # Validate chains in PDB
        atom_df = PdbDf(pdb_path)
        atom_df.get_atoms()
        atom_df_chains = atom_df.get_chains_id()
        if receptor_chains[0] not in atom_df_chains or receptor_chains[::-1][0] not in atom_df_chains:
            raise ValueError(f"{receptor_chains} is not a valid chain ID in {pdb_path}")
        if ligand_chain not in atom_df_chains:
            raise ValueError(f"{ligand_chain} is not a valid chain ID in {pdb_path}")
        if receptor_chains == ligand_chain or ligand_chain in receptor_chains:
            raise ValueError(f"Chain IDs should be different")
        if receptor_chains == ' ' or ligand_chain == ' ':
            raise ValueError(f"Chain IDs should not be empty")
        if len(receptor_chains) == 2 or len(receptor_chains) > 3 or len(ligand_chain) > 1:
            raise ValueError(f"Chain IDs should be 1 character for the ligand chain and 1 (letter) \
                            or 3 (letter,letter) characters for the receptor chain")       
                  
        
    def _get_chain_res(self, chain_id):
        
        """
        Get the residue numbers of a chain.
        
        Parameters:
        chain_id (str): Chain ID.
                        e.g. 'A'
        Returns:
        atom_df_ca_chain_lst (list): List of residue numbers of the chain.
                                     e.g. [1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        
        pdb = PdbDf(self.pdb_path)
        atom_df = pdb.get_atoms()
        atom_df_ca = pdb.get_atoms_ca()
        atom_df_ca_chain = PdbDf.get_atoms_chain(atom_df_ca, chain_id)
        atom_df_ca_chain_lst = atom_df_ca_chain['residue_number'].tolist()
                
        return atom_df_ca_chain_lst
        
        
    def _get_selected_res(self):
        
        """
        Get the selected residues from the interface analysis.
        
        Returns:
        interaction_rec_lst (list): List of interacting residues in the receptor
                                    e.g. ['A1', 'A2', 'A3'] in case of one chain ID
                                         ['A1', 'A2', ..., 'B1', 'B2', 'B3'] in case of two chain IDs
        interaction_lig_lst (list): List of interacting residues in the ligand
                                    e.g. ['B1', 'B2', 'B3']
        """

        ########## One receptor chain - One ligand chain
        if len(self.receptor_chains) == 1:
            
            analyzer = InterfaceAnalyzer(
                pdb_path=self.pdb_path,
                receptor_chain=self.receptor_chains,
                ligand_chain=self.ligand_chain,
                distance_cutOff=self.distance_cutOff
            )
        
            _ = analyzer.get_interaction_matrix()
            interaction_rec_df , interaction_lig_df = analyzer.expand_interface(neighborhood=3)
            self.interaction_rec_lst, self.interaction_lig_lst = analyzer.get_interacting_residues(mode='e')
                    
        ######### Two receptor chains - One ligand chain
        if len(self.receptor_chains) == 2:

            # Analyze both interfaces per sepparate
            analyzer_1 = InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=str(self.receptor_chain_1),
                    ligand_chain=self.ligand_chain,
                    distance_cutOff=self.distance_cutOff
                )

            analyzer_2= InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=str(self.receptor_chain_2),
                    ligand_chain=self.ligand_chain,
                    distance_cutOff=self.distance_cutOff
                )

            _ = analyzer_1.get_interaction_matrix()
            interaction_rec_1_df , interaction_lig_df = analyzer_1.expand_interface(neighborhood=3)
            interaction_rec_1_lst, interaction_lig_lst_1 = analyzer_1.get_interacting_residues(mode='e')
            
            _ = analyzer_2.get_interaction_matrix()
            interaction_rec_2_df , interaction_lig_df = analyzer_2.expand_interface(neighborhood=3)
            interaction_rec_2_lst, interaction_lig_lst_2 = analyzer_2.get_interacting_residues(mode='e')
            
            self.interaction_rec_lst = interaction_rec_1_lst + interaction_rec_2_lst # sum the lists of interacting residues in receptor
                
            # Combine the lists of interacting residues in ligand
            interaction_lig_lst = interaction_lig_lst_1 + interaction_lig_lst_2
            interaction_lig_lst_uniq = list(set(interaction_lig_lst))
            self.interaction_lig_lst =  sorted(interaction_lig_lst_uniq, key=lambda x: int(x[1:]))
        
        return self.interaction_rec_lst, self.interaction_lig_lst
    
    
    def _get_and_filter_chunks(self, res_lst):
        
        """
        Get the chunks of residues from the selected residues.
        
        Parameters:
        res_lst (list): List of selected residues.
                        e.g. ['A1', 'A2', 'A3', 'A15', 'A16', 'A17', 'A20', 'A21', 'A22']
                        
        Returns:
        res_lst_chunk (list): List of chunks of residues.
                              e.g. [[1, 2, 3], [5, 6, 7], [10, 11, 12]]
        """
        
        # Get the residue numbers without the chain ID
        res_lst_nochainID = [int(i[1:]) for i in res_lst]
        
        # Get the chunks of residues
        result = []
        subgroup = []

        for res_num in res_lst_nochainID:
            if not subgroup or res_num == subgroup[-1] + 1:
                subgroup.append(res_num)
            else:
                result.append(subgroup)
                subgroup = [res_num]
                
        if subgroup:
            result.append(subgroup)
            
        # Deletes isolate res_numbers e.g. 5 in [1,2,3], [5], [8, 9, 10]
        # since it makes no sense to diffuse one residue
        for i in result:
            if len(i)==1:
                result.remove(i)
        
        # Filter chunks shorter than 3
        result_filter = [i for i in result if len(i) >= 3]
        res_lst_chunk = result_filter
        
        return res_lst_chunk
    

    def _convert_to_contig_language(chain_id, chain_allres_lst, chain_selected_lst):
        """
        Given the residues of a chain and the residues of that chain that are selected
        as interacting, this function converts the information into a contig language.
        
        Parameters:
        chain_id (str): Chain ID of the chain.
                        e.g. 'A'
        chain_allres_lst (list): List of all residue NUMBERS of the chain.
                                 e.g. [1, 2, 3, 4, 5, 6, 7, 8, 9]
        chain_selected_lst (list): List of lists of selected residue NUMBERS grouped of the chain.
                                      e.g. [[1, 2, 3], [5, 6, 7], [10, 11, 12]]
        """
    
        chain_start = chain_allres_lst[0]
        chain_end = chain_allres_lst[-1]
        
        # Flat the list
        chain_selected_lst_extended = []
        for sublst in chain_selected_lst:
            chain_selected_lst_extended.extend(sublst)
        
        # IF THERE ARE SELECTED RESIDUES
        if chain_selected_lst:
            
            # All residues are selected
            if len(chain_selected_lst_extended) == len(chain_allres_lst):
                contig = f"{len(chain_allres_lst)}-{len(chain_allres_lst)}/"

            # First residue is selected but the last isn't
            elif chain_start in chain_selected_lst[0] and chain_end not in chain_selected_lst[-1]:
                contig = f"{len(chain_selected_lst[0])}-{len(chain_selected_lst[0])}/"
                contig += f"{chain_id}{chain_selected_lst[0][-1]+1}-"
                for curr_list in chain_selected_lst[1:]:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    contig += f"{chain_id}{curr_list_end+1}-"
                contig += f"{chain_end}"+"/"

            # First and last residue are selected
            elif chain_start in chain_selected_lst[0] and chain_end in chain_selected_lst[-1]:
                contig = f"{len(chain_selected_lst[0])}-{len(chain_selected_lst[0])}/"
                contig += f"{chain_id}{chain_selected_lst[0][-1]+1}-"
                for curr_list in chain_selected_lst[1:]:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    if curr_list == chain_selected_lst[-1]:
                        pass
                    else:
                        contig += f"{chain_id}{curr_list_end+1}-"

            # First residue isn't selected and last is selected
            elif chain_start not in chain_selected_lst[0] and chain_end in chain_selected_lst[-1]:
                contig = f"{chain_id}{chain_start}-"
                for curr_list in chain_selected_lst:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    if curr_list == chain_selected_lst[-1]:
                        pass
                    else:
                        contig += f"{chain_id}{curr_list_end+1}-"
                        
            # Neither the first nor the last residues are selected
            else:
                contig = f"{chain_id}{chain_start}-"
                for curr_list in chain_selected_lst:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    contig += f"{chain_id}{curr_list_end+1}-"
                contig += f"{chain_end}"+"/"

        # IF THERE ARE NO SELECTED RESIDUES
        else:
            contig = f"{chain_id}{chain_start}-{chain_end}"+"/"

        return contig
    
    
    def get_contigs(self):
        
        """
        Get the contigs of the system.
        
        Returns:
        contigs (str): Contigs of the system.
                       e.g. '[A1-3/A5-7/A10-12]'
        """
        
        # 2 chains ergo interface 
        # ------------------------
        if len(self.receptor_chains) == 1 and self.receptor_chain_1 == None and \
            self.receptor_chain_2 == None: # See init method for explanation
        
            # Get chainid, list of all residues and list of chunks of selected residues
            contigs = RFdiffContigs(
                pdb_path=self.pdb_path,
                receptor_chains=self.receptor_chains,
                ligand_chain=self.ligand_chain,
                distance_cutOff=self.distance_cutOff
                )
            
            # Get selected residues in each chain
            rec_selected, lig_selected = contigs._get_selected_res()
            
            # Get needed inputs for contig conversion
            rec_chainid = rec_selected[0][0]
            rec_allres = contigs._get_chain_res(rec_chainid)      
            rec_chunkres = contigs._get_and_filter_chunks(rec_selected)
            lig_chainid = lig_selected[0][0]    
            lig_allres = contigs._get_chain_res(lig_chainid)
            lig_chunkres = contigs._get_and_filter_chunks(lig_selected)
                    
            # Figure out order of chains in PDB. The order of the final contig need to coincide
            # with the order of the chains in the PDB file
            chains_input = [rec_chainid, lig_chainid]
            pdb = PdbDf(self.pdb_path)
            _ = pdb.get_atoms()
            chains_ordered = pdb.get_chains_id()

            if chains_input == chains_ordered: # if parsed and pdb order coincides, receptor goes 1st
                id1 = rec_chainid
                id2 = lig_chainid
                chain1_allres = rec_allres
                chain2_allres = lig_allres
                chain1_chunkres = rec_chunkres
                chain2_chunkres = lig_chunkres
            
            else: # if parsed and pdb order coincides, ligand goes first
                id1 = lig_chainid
                id2 = rec_chainid
                chain1_allres = lig_allres
                chain2_allres = rec_allres
                chain1_chunkres = lig_chunkres
                chain2_chunkres = rec_chunkres
                
            # Convert to contig language
            contigs_1 = RFdiffContigs._convert_to_contig_language(id1, chain1_allres, chain1_chunkres)
            contigs_2 = RFdiffContigs._convert_to_contig_language(id2, chain2_allres, chain2_chunkres)
                
            # Build contig
            contigs = f"[{contigs_1}0 {contigs_2[:-1]}]"
        
        # 3 chains ergo 2 interfaces 
        # ------------------------
        if len(self.receptor_chains) == 2 and self.receptor_chain_1 is not None and \
            self.receptor_chain_2 is not None: # See init method for explanation

            # Get chainid, list of all residues and list of chunks of selected residues
            contigs = RFdiffContigs(
                pdb_path=self.pdb_path,
                receptor_chains=f'{self.receptor_chain_1},{self.receptor_chain_2}', # this is done to pass e.g. 'A,B' rather than ['A', 'B']
                ligand_chain=self.ligand_chain,
                distance_cutOff=self.distance_cutOff
                )

            # Get selected residues in each chain, including splitting receptor selected list
            rec_selected, lig_selected = contigs._get_selected_res()
            rec_selected_1 = []
            rec_selected_2 = []
            for i in rec_selected:
                if i[0] == self.receptor_chain_1[0]:
                    rec_selected_1.append(i)
                else:
                    rec_selected_2.append(i)
            
            # Get needed inputs for contig conversion
            rec_chainid_1 = rec_selected_1[0][0]
            rec_allres_1 = contigs._get_chain_res(rec_chainid_1)      
            rec_chunkres_1 = contigs._get_and_filter_chunks(rec_selected_1)            
            
            rec_chainid_2 = rec_selected_2[0][0]
            rec_allres_2 = contigs._get_chain_res(rec_chainid_2)      
            rec_chunkres_2 = contigs._get_and_filter_chunks(rec_selected_2)                  
            
            lig_chainid = lig_selected[0][0]    
            lig_allres = contigs._get_chain_res(lig_chainid)
            lig_chunkres = contigs._get_and_filter_chunks(lig_selected)           
            
            # Convert to contig language
            contigs_1 = RFdiffContigs._convert_to_contig_language(rec_chainid_1, rec_allres_1, rec_chunkres_1)
            contigs_2 = RFdiffContigs._convert_to_contig_language(rec_chainid_2, rec_allres_2, rec_chunkres_2)
            contigs_3 = RFdiffContigs._convert_to_contig_language(lig_chainid, lig_allres, lig_chunkres)
            
            # Figure out order of chains in PDB. The order of the final contig need to coincide
            # with the order of the chains in the PDB file
            chains_input = [rec_chainid_1,rec_chainid_2, lig_chainid]
            pdb = PdbDf(self.pdb_path)
            _ = pdb.get_atoms()
            chains_ordered = pdb.get_chains_id()           
        
            if chains_ordered[0] in contigs_1:
                first = contigs_1
            if chains_ordered[0] in contigs_2:
                first = contigs_2
            if chains_ordered[0] in contigs_3:
                first = contigs_3
            if chains_ordered[1] in contigs_1:
                second = contigs_1
            if chains_ordered[1] in contigs_2:
                second = contigs_2
            if chains_ordered[1] in contigs_3:
                second = contigs_3
            if chains_ordered[2] in contigs_1:
                third = contigs_1
            if chains_ordered[2] in contigs_2:
                third = contigs_2
            if chains_ordered[2] in contigs_3:
                third = contigs_3
            
            # Build            
            contigs = f"[{first}0 {second}0 {third[:-1]}]"     

        return contigs


class RFdiffSetUp:

    def __init__(self,
                 pdb_path: str,
                 run_inference_path: str = "/gpfs/projects/bsc72/Repos/RFdiffusion/scripts/run_inference.py"
                 ):
        
        """
        Initialize the RFdiffSetUp class.
        This class pretends to generate the runner for RFdiffusion (partial diffusion)
        This runner is adapted to be used in our cluster MareNostrum V
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        run_inference_path (str): Path to the run_inference.py script from RFdiffusion
        """
        
        # Validate inputs
        self._validate_inputs(pdb_path, 
                              run_inference_path)
        
        # Set attributes
        self.pdb_path = pdb_path
        self.run_inference_path = run_inference_path

    def _validate_inputs(
        self,
        pdb_path: str,
        run_inference_path: str
    ):
        
        # Validate string inputs
        for var_name, var_value in [
            ("pdb_path", pdb_path),
            ("run_inference_path", run_inference_path)
            ]:
            if not isinstance(var_value, str):
                raise TypeError(f"{var_name} should be a string")
        
        # Validate PDB
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"File not found: {pdb_path}")
        #if not os.path.exists(run_inference_path):
         #   raise FileNotFoundError(f"File not found: {run_inference_path}")
    
    def make_runner(self):
        
        """
        Build RFdiffusion partial diffusion runner for MNV
        
        Returns:
        runner_path (str): Path to the runner script
        """
        
        content = r"""
        #!/bin/bash
        #SBATCH --job-name __DJRN__
        #SBATCH -D .
        #SBATCH --output=logs/%x_%j.out
        #SBATCH --error=logs/%x_%j.err
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task=20
        #SBATCH --gres=gpu:1
        #SBATCH --time=01:00:00
        #SBATCH --qos=acc_bscls
        #SBATCH --nodes=1
        #SBATCH --account=bsc72

        # Changes to find ACC's modules
        module purge
        export MODULEPATH=/apps/ACC/modulefiles/applications:/apps/ACC/modulefiles/compilers:/apps/ACC/modulefiles/environment:/apps/ACC/modulefiles/lib

        # Environment stuff
        module load anaconda/2024.02
        eval "$(conda shell.bash hook)"
        conda activate RFDiffusion

        # Parameters
        complex_pdb=__PDB__ # PDB hbond-opt, relaxed, fixed (removed insertions)
        run_inference=__RINF__ # Run inference path

        # Organize dir
        sys_basename=$(basename "$complex_pdb" .pdb)
        w_dir=$(dirname "$complex_pdb")
        contigs_file="${w_dir}/${sys_basename}_contigs.txt"
        out_prefix_noise0="${w_dir}/${sys_basename}_diff/N0/${sys_basename}_diff_N0"
        out_prefix_noise05="${w_dir}/${sys_basename}_diff/N05/${sys_basename}_diff_N05"
        out_prefix_noise075="${w_dir}/${sys_basename}_diff/N075/${sys_basename}_diff_N075"
        out_prefix_noise1="${w_dir}/${sys_basename}_diff/N1/${sys_basename}_diff_N1"

        # RFdiffusion stuff
        n_structs=10
        n_diff_iter=20
        fix_contigs=$(cat ${contigs_file})
        seq_cov="[0-1000]" # Selects all aminoacids to keep seq (it doesn't complain if this selection is longer)
        echo Partial Diffusion for $sys_basename
        echo Reading contigs file $contigs_file

        # No Noise
        ${run_inference} inference.input_pdb=$complex_pdb \
                         inference.output_prefix=$out_prefix_noise0 \
                         "contigmap.contigs=$fix_contigs" \
                         "contigmap.provide_seq=$seq_cov" \
                         inference.num_designs=$n_structs \
                         denoiser.noise_scale_ca=0 \
                         denoiser.noise_scale_frame=0 \
                         diffuser.partial_T=$n_diff_iter

        #  Noisy
        ${run_inference} inference.input_pdb=$complex_pdb \
                         inference.output_prefix=$out_prefix_noise05 \
                         "contigmap.contigs=$fix_contigs" \
                         "contigmap.provide_seq=$seq_cov" \
                         inference.num_designs=$n_structs \
                         denoiser.noise_scale_ca=0.5 \
                         denoiser.noise_scale_frame=0.5 \
                         diffuser.partial_T=$n_diff_iter

        #  Noisy
        ${run_inference} inference.input_pdb=$complex_pdb \
                         inference.output_prefix=$out_prefix_noise075 \
                         "contigmap.contigs=$fix_contigs" \
                         "contigmap.provide_seq=$seq_cov" \
                         inference.num_designs=$n_structs \
                         denoiser.noise_scale_ca=0.75 \
                         denoiser.noise_scale_frame=0.75 \
                         diffuser.partial_T=$n_diff_iter

        #  Noisy
        ${run_inference} inference.input_pdb=$complex_pdb \
                         inference.output_prefix=$out_prefix_noise1 \
                         "contigmap.contigs=$fix_contigs" \
                         "contigmap.provide_seq=$seq_cov" \
                         inference.num_designs=$n_structs \
                         denoiser.noise_scale_ca=1.0 \
                         denoiser.noise_scale_frame=1.0 \
                         diffuser.partial_T=$n_diff_iter

        # Copy all diffusion models
        diff_dir="${w_dir}/${sys_basename}_diffModels"
        mkdir -p $diff_dir

        cp ${w_dir}/${sys_basename}_diff/N0/*pdb $diff_dir
        cp ${w_dir}/${sys_basename}_diff/N05/*pdb $diff_dir
        cp ${w_dir}/${sys_basename}_diff/N075/*pdb $diff_dir
        cp ${w_dir}/${sys_basename}_diff/N1/*pdb $diff_dir
        """
        content = textwrap.dedent(content)
        content = content.lstrip("\n")

        content = content.replace("__PDB__", os.path.abspath(self.pdb_path))
        content = content.replace("__RINF__", os.path.abspath(self.run_inference_path))
        content = content.replace("__DJRN__", os.path.basename(self.pdb_path)[:-4]+'_diff')
        
        return content        


class RFdiffFix:
    
    def __init__(
                self, 
                original_pdb_path: str, 
                diffusion_models_input: str
                ):
        """
        Correct PDB format of diffusion models according
        to the original input. It re-adds chain IDs and 
        residue numbers.
        
        Parameters:
        original_pdb_path (str): Path to the original PDB file.
        diffusion_models_input (str): Path to the diffusion models (path/to/diffusion_models/)
                                      OR
                                      a single diffusion model (path/to/diffusion_model.pdb)
                                      This is automatically handled by the class.
        """
        # Validate inputs
        self._validate_inputs(original_pdb_path, diffusion_models_input)
        
        # Set attributes
        self.original_pdb_path = original_pdb_path
        self.diffusion_models_input = diffusion_models_input
    
    
    def _validate_inputs(
        self,
        original_pdb_path: str,
        diffusion_models_input: str
    ):
            
        # Validate string inputs
        for var_name, var_value in [
            ("original_pdb_path", original_pdb_path),
            ("diffusion_models_input", diffusion_models_input)
            ]:
            if not isinstance(var_value, str):
                raise TypeError(f"{var_name} should be a string")
            
        # Validate PDB
        if not os.path.exists(original_pdb_path):
            raise FileNotFoundError(f"File not found: {original_pdb_path}")
        if diffusion_models_input.endswith('.pdb'):
            if not os.path.exists(diffusion_models_input):
                raise FileNotFoundError(f"File not found: {diffusion_models_input}")
        if diffusion_models_input.endswith('.pdb') == False:
            if not os.path.exists(diffusion_models_input):
                raise FileNotFoundError(f"Directory not found: {diffusion_models_input}")


    def correct_diff_models(self):

        
        if self.diffusion_models_input.endswith('.pdb'): # the input is a single pdb file
            
            logger.info("Correcting a single diffusion model")
            
            # Reduce original PDB to backbone atoms
            PdbHandler.get_backbone(self.original_pdb_path)
            PdbHandler.get_atom(self.original_pdb_path[:-4]+'_backbone.pdb')
            
            # Read as PDB dataframe original and diffusion PDB
            original_pdb_backbone_atoms_df = PdbDf(self.original_pdb_path[:-4]+'_backbone_atom.pdb').get_atoms()
            diffusion_model_atoms_df = PdbDf(self.diffusion_models_input).get_atoms()
            
            if len(original_pdb_backbone_atoms_df) != len(diffusion_model_atoms_df):
                raise ValueError(f"Original PDB and diffusion model have different number of atoms: \
                                {len(original_pdb_backbone_atoms_df)} vs {len(diffusion_model_atoms_df)}")
            
            # Now the diffusion model and original system have the same atom lines in the same order, thus
            # chain ids and residue numbers can be added to the diffusion model(s)
            diffusion_model_atoms_df['chain_id'] = original_pdb_backbone_atoms_df['chain_id']
            diffusion_model_atoms_df['residue_number'] = original_pdb_backbone_atoms_df['residue_number']
            
            # Write to PDB
            diffusion_model_fix = PandasPdb().read_pdb(self.diffusion_models_input)
            diffusion_model_fix.df['ATOM'] = diffusion_model_atoms_df
            diffusion_model_fix.to_pdb(path= self.diffusion_models_input[:-4]+'_tmpfix.pdb', records = ['ATOM'], gz=False)  
            
            # The resulting PDB needs TER and END lines
            PdbHandler.tidy(self.diffusion_models_input[:-4]+'_tmpfix.pdb')
            
            # Remove temporal files
            for i in [
                self.original_pdb_path[:-4]+'_backbone.pdb',
                self.original_pdb_path[:-4]+'_backbone_atom.pdb',
                self.diffusion_models_input[:-4]+'_tmpfix.pdb'
                ]:
                os.remove(i)
            
            # Rename the fixed diffusion model
            os.rename(
                self.diffusion_models_input[:-4]+'_tmpfix_tidied.pdb',
                self.diffusion_models_input[:-4]+'_fix.pdb'
            )
    
        else: # the input is a path to a directory with diffusion models
            
            logger.info("Correcting multiple diffusion models")
            
            items = os.listdir(self.diffusion_models_input)
            diffmodels_to_correct = [os.path.join(self.diffusion_models_input,pdb) for pdb in items if pdb.startswith(os.path.basename(self.original_pdb_path)[:-4]+'_diff')]
            
            if len(diffmodels_to_correct) == 0:
                raise ValueError('No diffusion models detected, diffusion models names must start be original_pdb_name_diff')
            
            # Reduce original PDB to backbone atoms and read as dataframe
            PdbHandler.get_backbone(self.original_pdb_path)
            PdbHandler.get_atom(self.original_pdb_path[:-4]+'_backbone.pdb')
            original_pdb_backbone_atoms_df = PdbDf(self.original_pdb_path[:-4]+'_backbone_atom.pdb').get_atoms()
            
            # Iterate over the diffusion models and correct
            for diffmodel in diffmodels_to_correct:
                
                diffmodel_df = PdbDf(diffmodel).get_atoms()

                if len(original_pdb_backbone_atoms_df) != len(diffmodel_df):
                    raise ValueError(f"Original PDB and diffusion model have different number of atoms: \
                                    {len(original_pdb_backbone_atoms_df)} vs {len(diffmodel_df)}")
                    
                diffmodel_df['chain_id'] = original_pdb_backbone_atoms_df['chain_id']
                diffmodel_df['residue_number'] = original_pdb_backbone_atoms_df['residue_number']
                
                diffmodel_fix = PandasPdb().read_pdb(diffmodel)
                diffmodel_fix.df['ATOM'] = diffmodel_df
                diffmodel_fix.to_pdb(path= diffmodel[:-4]+'_tmpfix.pdb', records = ['ATOM'], gz=False)  
                
                PdbHandler.tidy(diffmodel[:-4]+'_tmpfix.pdb')
                
                os.remove(diffmodel[:-4]+'_tmpfix.pdb')
                os.rename(diffmodel[:-4]+'_tmpfix_tidied.pdb', 
                        diffmodel[:-4]+'_fix.pdb')
                
            for i in [
                self.original_pdb_path[:-4]+'_backbone.pdb',
                self.original_pdb_path[:-4]+'_backbone_atom.pdb'
                ]:
                os.remove(i)

