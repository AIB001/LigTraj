# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# """
# GAFFMaker Python - A tool for building protein-ligand systems for GROMACS MD simulations

# This script streamlines the process of modeling protein-ligand molecular dynamics simulations
# with GROMACS by automating the preparation of force field parameters using AmberTools
# and setting up the simulation system.

# Usage:
#     python GAFF_Maker.py <protein_path> <ligand_path> --output <output_dir>

# Author: Adapted from GAFFMakerV1.1 bash scripts
# """

# import os
# import sys
# import subprocess
# import shutil
# import re
# import tempfile
# import argparse
# import glob
# import time
# from pathlib import Path


# class GAFFMaker:
#     def __init__(self, protein_path, ligand_path, output_dir):
#         """
#         Initialize the GAFFMaker with protein and ligand paths
        
#         Parameters:
#         -----------
#         protein_path : str
#             Path to the protein PDB file
#         ligand_path : str
#             Path to the ligand MOL2 file
#         output_dir : str
#             Directory where output files will be stored
#         """
#         self.protein_path = os.path.abspath(protein_path)
#         self.ligand_path = os.path.abspath(ligand_path)
#         self.output_dir = os.path.abspath(output_dir)
        
#         # Extract filenames without extensions
#         self.protein_name = os.path.basename(protein_path).split('.')[0]
#         self.ligand_name = os.path.basename(ligand_path).split('.')[0]
#         self.ligand_id = self._extract_ligand_id()
#         # 存储所有可能的分子名称变体，用于替换
#         self.possible_molecule_names = [
#             self.ligand_id,           # 从MOL2提取的ID (例如 1unl_ligand)
#             self.ligand_name,         # 文件名 (例如 ligand)
#             "RRC",                    # acpype可能生成的名称
#             self.ligand_id.upper(),   
#             self.ligand_name.upper(),
#             "LIG"                     # 标准化目标名称
#         ]

#         # Create output directory if it doesn't exist
#         os.makedirs(self.output_dir, exist_ok=True)
        
#         # Directory for MDP files
#         self.mdps_dir = os.path.join(self.output_dir, "mdps")
#         os.makedirs(self.mdps_dir, exist_ok=True)
        
#         print(f"Initialized GAFFMaker for protein: {self.protein_path}")
#         print(f"                       and ligand: {self.ligand_path}")
#         print(f"                   Output directory: {self.output_dir}")
        
#         # Debug level
#         self.debug = True

#     def _extract_ligand_id(self):
#         """Extract ligand ID from MOL2 file"""
#         try:
#             with open(self.ligand_path, 'r') as f:
#                 for line in f:
#                     if line.startswith("@<TRIPOS>MOLECULE"):
#                         # Read the next line which contains the molecule name
#                         ligand_id = next(f).strip()
#                         print(f"Extracted ligand ID from mol2 file: {ligand_id}")
#                         return ligand_id
#         except Exception as e:
#             print(f"Failed to extract ligand ID from MOL2 file: {e}")
        
#         # Fallback to filename
#         return self.ligand_name

#     def run_command(self, cmd, cwd=None, env=None, check=True, input_text=None, shell=False):
#         """
#         Run a shell command and return its output
        
#         Parameters:
#         -----------
#         cmd : str or list
#             Command to run
#         cwd : str, optional
#             Working directory
#         env : dict, optional
#             Environment variables
#         check : bool
#             Whether to check for errors
#         input_text : str, optional
#             Text to pass to stdin
#         shell : bool
#             Whether to use shell to execute command
            
#         Returns:
#         --------
#         str
#             Command output
#         """
#         print(f"Running command: {cmd}")
        
#         if isinstance(cmd, str):
#             cmd_str = cmd
#             if not shell:
#                 cmd_list = cmd.split()
#         else:
#             cmd_list = cmd
#             cmd_str = " ".join(cmd)
            
#         try:
#             # 对于需要通过管道提供输入的命令使用Popen
#             if input_text is not None and ("genion" in cmd_str or "pdb2gmx" in cmd_str):
#                 print(f"Using Popen with input: {input_text}")
#                 process = subprocess.Popen(
#                     cmd_str if shell else cmd_list,
#                     cwd=cwd,
#                     env=env,
#                     shell=shell,
#                     stdin=subprocess.PIPE,
#                     stdout=subprocess.PIPE,
#                     stderr=subprocess.PIPE,
#                     text=True
#                 )
#                 stdout, stderr = process.communicate(input=input_text)
                
#                 if process.returncode != 0 and check:
#                     print(f"Error executing command: {cmd}")
#                     print(f"Return code: {process.returncode}")
#                     print(f"Error output: {stderr}")
#                     raise subprocess.CalledProcessError(process.returncode, cmd, stdout, stderr)
                
#                 print(f"Command completed with return code: {process.returncode}")
#                 if stdout.strip():
#                     print(f"Output snippet: {stdout.strip()[:100]}...")
#                 return stdout
#             else:
#                 # 对于其他命令使用run
#                 kwargs = {
#                     'cwd': cwd,
#                     'env': env,
#                     'check': check,
#                     'text': True,
#                     'stdout': subprocess.PIPE,
#                     'stderr': subprocess.PIPE,
#                     'shell': shell
#                 }
                
#                 # 只有在非特殊命令时添加input
#                 if input_text is not None:
#                     kwargs['input'] = input_text
                    
#                 result = subprocess.run(
#                     cmd_str if shell else cmd_list,
#                     **kwargs
#                 )
                
#                 print(f"Command completed with return code: {result.returncode}")
#                 if result.stdout.strip():
#                     print(f"Output snippet: {result.stdout.strip()[:100]}...")
#                 return result.stdout
                
#         except subprocess.CalledProcessError as e:
#             print(f"Error executing command: {cmd}")
#             print(f"Return code: {e.returncode}")
#             print(f"Error output: {e.stderr}")
#             if check:
#                 raise
#             return e.stdout

#     def clean_protein(self):
#         """
#         Clean the protein PDB file by removing HETATM lines
        
#         Returns:
#         --------
#         str
#             Path to the cleaned protein PDB file
#         """
#         print("Cleaning protein PDB file...")
        
#         # Create path for cleaned PDB
#         cleaned_pdb = os.path.join(self.output_dir, f"{self.protein_name}_clean.pdb")
        
#         # Read PDB file and remove HETATM lines
#         with open(self.protein_path, 'r') as f_in, open(cleaned_pdb, 'w') as f_out:
#             for line in f_in:
#                 if not line.startswith('HETATM'):
#                     f_out.write(line)
        
#         # Fix HN1/HN2/HN3 atom names
#         with open(cleaned_pdb, 'r') as f:
#             content = f.read()
        
#         content = content.replace('HN1', 'H1 ').replace('HN2', 'H2 ').replace('HN3', 'H3 ')
        
#         with open(cleaned_pdb, 'w') as f:
#             f.write(content)
        
#         print(f"Protein cleaned and saved to {cleaned_pdb}")
#         return cleaned_pdb

#     def generate_ligand_forcefield(self):
#         """
#         Generate force field parameters for the ligand
        
#         Returns:
#         --------
#         dict
#             Paths to the generated files
#         """
#         print("Generating ligand force field parameters...")
        
#         # Create forcefield directory
#         ff_dir = os.path.join(self.output_dir, "forcefield")
#         os.makedirs(ff_dir, exist_ok=True)
        
#         # Paths for intermediate files
#         amber_mol2 = os.path.join(ff_dir, f"{self.ligand_name}.mol2")
#         prep_file = os.path.join(ff_dir, f"{self.ligand_name}.prep")
#         frcmod_file = os.path.join(ff_dir, f"{self.ligand_name}.frcmod")
#         prmtop_file = os.path.join(ff_dir, f"{self.ligand_name}.prmtop")
#         rst7_file = os.path.join(ff_dir, f"{self.ligand_name}.rst7")
        
#         # Step 1: Generate Amber-compatible mol2
#         self.run_command([
#             "antechamber", 
#             "-i", self.ligand_path, 
#             "-fi", "mol2", 
#             "-o", amber_mol2, 
#             "-fo", "mol2", 
#             "-c", "bcc", 
#             "-s", "2"
#         ])
        
#         # Step 2: Generate prep file
#         self.run_command([
#             "antechamber", 
#             "-i", amber_mol2, 
#             "-fi", "mol2", 
#             "-o", prep_file, 
#             "-fo", "prepi", 
#             "-c", "bcc", 
#             "-s", "2"
#         ])
        
#         # Step 3: Generate frcmod
#         self.run_command([
#             "parmchk2", 
#             "-i", prep_file, 
#             "-f", "prepi", 
#             "-o", frcmod_file
#         ])
        
#         # Step 4: Prepare tleap input
#         tleap_input = os.path.join(ff_dir, f"tleap_{self.ligand_name}.in")
#         with open(tleap_input, 'w') as f:
#             f.write(f"""source leaprc.protein.ff14SB
# source leaprc.gaff

# LIG = loadmol2 {amber_mol2}
# loadamberparams {frcmod_file}

# saveamberparm LIG {prmtop_file} {rst7_file}

# quit
# """)
        
#         # Step 5: Run tleap
#         self.run_command([
#             "tleap", 
#             "-f", tleap_input
#         ])
        
#         # Step 6: Convert to GROMACS format with acpype
#         current_dir = os.getcwd()
#         os.chdir(ff_dir)
        
#         # 开启调试模式
#         if self.debug:
#             print(f"Before acpype, current directory: {os.getcwd()}")
#             print(f"Directory contents: {os.listdir('.')}")
        
#         # 运行acpype
#         self.run_command([
#             "acpype", 
#             "-p", os.path.basename(prmtop_file), 
#             "-x", os.path.basename(rst7_file), 
#             "-d"
#         ])
        
#         # 识别acpype的输出目录
#         if self.debug:
#             print(f"After acpype, current directory: {os.getcwd()}")
#             print(f"Directory contents: {os.listdir('.')}")
        
#         # 尝试多种可能的acpype输出目录名
#         possible_acpype_dirs = [
#             os.path.join(ff_dir, f"{self.ligand_name}.acpype"),  # ligand.acpype
#             os.path.join(ff_dir, f"{self.ligand_id}.acpype"),    # 1unl_ligand.acpype
#             os.path.join(ff_dir, "RRC.acpype"),                  # RRC.acpype 
#             os.path.join(ff_dir, "RRC.amb2gmx"),                 # RRC.amb2gmx
#             os.path.join(ff_dir, f"{self.ligand_name.upper()}.acpype"), 
#             os.path.join(ff_dir, f"{self.ligand_id.upper()}.acpype"),
#             os.path.join(ff_dir, f"{self.ligand_name.upper()}.amb2gmx"),
#             os.path.join(ff_dir, f"{self.ligand_id.upper()}.amb2gmx")
#         ]
        
#         acpype_dir = None
#         for dir_path in possible_acpype_dirs:
#             if os.path.exists(dir_path):
#                 acpype_dir = dir_path
#                 print(f"Found acpype output directory: {acpype_dir}")
#                 if self.debug:
#                     print(f"Contents of acpype directory: {os.listdir(acpype_dir)}")
#                 break
                
#         if not acpype_dir:
#             print(f"Warning: Could not find acpype output directory in expected locations")
#             # 搜索可能的acpype/amb2gmx目录
#             possible_dirs = [d for d in os.listdir('.') if (d.endswith('.acpype') or d.endswith('.amb2gmx')) and os.path.isdir(d)]
#             if possible_dirs:
#                 acpype_dir = os.path.join(ff_dir, possible_dirs[0])
#                 print(f"Using alternative acpype directory: {acpype_dir}")
#                 if self.debug:
#                     print(f"Contents of alternative acpype directory: {os.listdir(acpype_dir)}")
        
#         # 直接在acpype目录中查找文件
#         itp_file = None
#         gro_file = None
#         top_file = None
#         posre_file = None
        
#         if acpype_dir and os.path.exists(acpype_dir):
#             # 在acpype目录中查找文件
#             for filename in os.listdir(acpype_dir):
#                 filepath = os.path.join(acpype_dir, filename)
#                 if filename.endswith('_GMX.itp'):
#                     itp_file = filepath
#                     print(f"Found ITP file: {filepath}")
#                 elif filename.endswith('_GMX.gro'):
#                     gro_file = filepath
#                     print(f"Found GRO file: {filepath}")
#                 elif filename.endswith('_GMX.top'):
#                     top_file = filepath
#                     print(f"Found TOP file: {filepath}")
#                 elif filename.startswith('posre_') and filename.endswith('.itp'):
#                     posre_file = filepath
#                     print(f"Found POSRE file: {filepath}")
        
#         # 如果没有找到文件，进行递归搜索
#         if not all([gro_file, top_file]):
#             print("Warning: Some required files not found. Performing recursive search...")
#             for root, dirs, files in os.walk(ff_dir):
#                 for file in files:
#                     if file.endswith('_GMX.itp') and not itp_file:
#                         itp_file = os.path.join(root, file)
#                         print(f"Found ITP file in recursive search: {itp_file}")
#                     elif file.endswith('_GMX.gro') and not gro_file:
#                         gro_file = os.path.join(root, file)
#                         print(f"Found GRO file in recursive search: {gro_file}")
#                     elif file.endswith('_GMX.top') and not top_file:
#                         top_file = os.path.join(root, file)
#                         print(f"Found TOP file in recursive search: {top_file}")
#                     elif file.startswith('posre_') and file.endswith('.itp') and not posre_file:
#                         posre_file = os.path.join(root, file)
#                         print(f"Found POSRE file in recursive search: {posre_file}")
        
#         ff_files = {
#             "itp_file": itp_file,
#             "gro_file": gro_file,
#             "top_file": top_file,
#             "posre_file": posre_file
#         }
        
#         os.chdir(current_dir)
        
#         print(f"Generated ligand force field files: {ff_files}")
#         return ff_files

#     def standardize_files(self, ff_files):
#         """
#         Standardize file names and directory structure
        
#         Parameters:
#         -----------
#         ff_files : dict
#             Dictionary with paths to force field files
        
#         Returns:
#         --------
#         dict
#             Paths to standardized files
#         """
#         print("Standardizing files...")
        
#         # Create LIG.amb2gmx directory
#         lig_dir = os.path.join(self.output_dir, "LIG.amb2gmx")
#         os.makedirs(lig_dir, exist_ok=True)
        
#         # Copy and rename files
#         std_files = {}
        
#         # Process gro file
#         if ff_files["gro_file"] and os.path.exists(ff_files["gro_file"]):
#             std_files["gro_file"] = os.path.join(lig_dir, "LIG.gro")
#             shutil.copy(ff_files["gro_file"], std_files["gro_file"])
#             print(f"Copied {ff_files['gro_file']} to {std_files['gro_file']}")
            
#             # 替换所有可能的分子名称为LIG
#             for mol_name in self.possible_molecule_names:
#                 if mol_name != "LIG":  # 不替换已经是LIG的名称
#                     self._replace_in_file(std_files["gro_file"], mol_name, "LIG")
#             print(f"Replaced molecule names with 'LIG' in {std_files['gro_file']}")
#         else:
#             print(f"Warning: GRO file not found, skipping. Path was: {ff_files['gro_file']}")
        
#         # Process top file
#         if ff_files["top_file"] and os.path.exists(ff_files["top_file"]):
#             std_files["top_file"] = os.path.join(lig_dir, "LIG.top")
#             shutil.copy(ff_files["top_file"], std_files["top_file"])
#             print(f"Copied {ff_files['top_file']} to {std_files['top_file']}")
            
#             # 替换所有可能的分子名称为LIG
#             for mol_name in self.possible_molecule_names:
#                 if mol_name != "LIG":
#                     self._replace_in_file(std_files["top_file"], mol_name, "LIG")
#             print(f"Replaced molecule names with 'LIG' in {std_files['top_file']}")
#         else:
#             print(f"Warning: TOP file not found, skipping. Path was: {ff_files['top_file']}")
        
#         # 创建或处理ITP文件
#         # 如果存在ITP文件，直接复制并修改；如果不存在，从TOP文件创建
#         if ff_files["itp_file"] and os.path.exists(ff_files["itp_file"]):
#             # 正常处理已存在的ITP文件
#             std_files["itp_file"] = os.path.join(lig_dir, "LIG.itp")
#             self._create_itp_from_existing(ff_files["itp_file"], std_files["itp_file"])
#             print(f"Created {std_files['itp_file']} from existing {ff_files['itp_file']}")
#         elif ff_files["top_file"] and os.path.exists(ff_files["top_file"]):
#             # 从TOP文件创建ITP文件
#             std_files["itp_file"] = os.path.join(lig_dir, "LIG.itp")
#             self._create_itp_from_top(ff_files["top_file"], std_files["itp_file"])
#             print(f"Created {std_files['itp_file']} from TOP file {ff_files['top_file']}")
#         else:
#             print(f"Warning: Cannot create ITP file - no source file available")
        
#         # Process posre file
#         if ff_files["posre_file"] and os.path.exists(ff_files["posre_file"]):
#             std_files["posre_file"] = os.path.join(lig_dir, "posre_LIG.itp")
#             shutil.copy(ff_files["posre_file"], std_files["posre_file"])
            
#             # 替换所有可能的分子名称为LIG
#             for mol_name in self.possible_molecule_names:
#                 if mol_name != "LIG":
#                     self._replace_in_file(std_files["posre_file"], mol_name, "LIG")
#             print(f"Created {std_files['posre_file']} from {ff_files['posre_file']} and replaced molecule names with 'LIG'")
#         else:
#             print(f"Warning: POSRE file not found, skipping. Path was: {ff_files['posre_file']}")
        
#         # Extract atomtypes block from the top file
#         if ff_files["top_file"] and os.path.exists(ff_files["top_file"]):
#             with open(ff_files["top_file"], 'r') as f:
#                 content = f.read()
            
#             # 使用更宽松的正则表达式查找atomtypes部分
#             atomtypes_match = re.search(r'\[ *atomtypes *\](.*?)(\[ *[a-z]|$)', content, re.DOTALL)
#             if atomtypes_match:
#                 atomtypes_block = "[ atomtypes ]\n" + atomtypes_match.group(1).strip()
#                 std_files["atomtypes_file"] = os.path.join(lig_dir, "atomtypes_block.itp")
#                 with open(std_files["atomtypes_file"], 'w') as f:
#                     f.write(atomtypes_block)
#                 print(f"Created atomtypes file {std_files['atomtypes_file']}")
#             else:
#                 print(f"Warning: Could not extract atomtypes block from TOP file: {ff_files['top_file']}")
#                 # 尝试创建一个空的atomtypes文件，以便后续流程能继续
#                 std_files["atomtypes_file"] = os.path.join(lig_dir, "atomtypes_block.itp")
#                 with open(std_files["atomtypes_file"], 'w') as f:
#                     f.write("[ atomtypes ]\n; Empty atomtypes file created because extraction failed\n")
#                 print(f"Created empty atomtypes file as fallback")
        
#         # 检查所有必需的文件是否都生成了
#         required_files = ["gro_file", "itp_file", "atomtypes_file"]
#         missing_files = [f for f in required_files if f not in std_files]
        
#         if missing_files:
#             print(f"WARNING: The following required files are missing: {missing_files}")
#             print("This may cause problems in model building!")
            
#         print(f"Standardized files created in {lig_dir}")
#         return std_files
    
#     def _create_itp_from_existing(self, source_itp, target_itp):
#         """从现有的ITP文件创建标准化的ITP文件"""
#         with open(source_itp, 'r') as f_in, open(target_itp, 'w') as f_out:
#             content = f_in.read()
            
#             # 替换所有可能的分子名称为LIG
#             for mol_name in self.possible_molecule_names:
#                 if mol_name != "LIG":
#                     content = re.sub(r'\b{}\b'.format(re.escape(mol_name)), "LIG", content)
            
#             # 确保文件末尾有posre包含
#             if "#ifdef POSRES" not in content:
#                 content += "\n#ifdef POSRES\n"
#                 content += '#include "posre_LIG.itp"\n'
#                 content += "#endif\n"
            
#             f_out.write(content)
    
#     def _create_itp_from_top(self, top_file, itp_file):
#         """从TOP文件中提取分子部分以创建ITP文件"""
#         try:
#             with open(top_file, 'r') as f:
#                 content = f.read()
            
#             # 查找moleculetype部分直到system部分（或文件结束）
#             mol_match = re.search(r'\[ *moleculetype *\](.*?)(\[ *system|$)', content, re.DOTALL)
            
#             if mol_match:
#                 # 提取分子部分
#                 mol_content = "[ moleculetype ]\n" + mol_match.group(1).strip()
                
#                 # 替换所有可能的分子名称为LIG
#                 for mol_name in self.possible_molecule_names:
#                     if mol_name != "LIG":
#                         mol_content = re.sub(r'\b{}\b'.format(re.escape(mol_name)), "LIG", mol_content)
                
#                 # 添加posre包含
#                 mol_content += "\n\n#ifdef POSRES\n"
#                 mol_content += '#include "posre_LIG.itp"\n'
#                 mol_content += "#endif\n"
                
#                 with open(itp_file, 'w') as f:
#                     f.write(mol_content)
                
#                 return True
#             else:
#                 # 如果找不到moleculetype部分，创建一个基本的
#                 print(f"Warning: Could not find moleculetype section in {top_file}")
#                 basic_itp = """[ moleculetype ]
# ; Name            nrexcl
# LIG              3

# [ atoms ]
# ; Note: This is a placeholder. The actual ITP content could not be extracted.
# ; Please check your ligand structure and force field parameters.

# #ifdef POSRES
# #include "posre_LIG.itp"
# #endif
# """
#                 with open(itp_file, 'w') as f:
#                     f.write(basic_itp)
                
#                 return False
                
#         except Exception as e:
#             print(f"Error creating ITP from TOP file: {e}")
#             return False
    
#     def _replace_in_file(self, file_path, old_str, new_str):
#         """替换文件中的文本"""
#         try:
#             with open(file_path, 'r') as f:
#                 content = f.read()
            
#             # 使用正则表达式替换，确保只替换完整的单词
#             updated_content = re.sub(r'\b{}\b'.format(re.escape(old_str)), new_str, content)
            
#             # 如果内容有变化，才写回文件
#             if updated_content != content:
#                 with open(file_path, 'w') as f:
#                     f.write(updated_content)
#                 return True
#             return False
#         except Exception as e:
#             print(f"Error replacing text in {file_path}: {e}")
#             return False

#     def build_model(self, cleaned_protein, std_files):
#         """
#         Build the GROMACS model
        
#         Parameters:
#         -----------
#         cleaned_protein : str
#             Path to the cleaned protein PDB file
#         std_files : dict
#             Dictionary with paths to standardized files
        
#         Returns:
#         --------
#         str
#             Path to the model directory
#         """
#         print("Building GROMACS model...")
        
#         # 检查必需文件
#         required_files = {
#             "gro_file": "Ligand GRO file",
#             "itp_file": "Ligand ITP file", 
#             "atomtypes_file": "Atomtypes file"
#         }
        
#         missing_files = []
#         for file_key, file_desc in required_files.items():
#             if file_key not in std_files or not std_files[file_key] or not os.path.exists(std_files[file_key]):
#                 missing_files.append(f"{file_desc} ({file_key})")
                
#         if missing_files:
#             error_msg = "Missing required files: " + ", ".join(missing_files)
#             print(f"ERROR: {error_msg}")
#             return None
            
#         # Create model directory
#         model_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
#         os.makedirs(model_dir, exist_ok=True)
        
#         # Step 1: Run pdbfixer to clean the protein further
#         fixed_pdb = os.path.join(model_dir, "fixed_clean_protein.pdb")
#         self.run_command([
#             "pdbfixer", 
#             cleaned_protein, 
#             "--output", fixed_pdb, 
#             "--add-atoms", "heavy", 
#             "--keep-heterogens", "none"
#         ])
        
#         # Step 2: Generate topology using gmx pdb2gmx
#         current_dir = os.getcwd()
#         os.chdir(model_dir)
        
#         # 直接使用Popen运行pdb2gmx
#         self.run_command(
#             ["gmx", "pdb2gmx", "-f", "fixed_clean_protein.pdb", "-o", "pro_lig.gro", "-ignh"],
#             input_text="6\n1\n"  # 6 for AMBER99SB, 1 for TIP3P water
#         )
        
#         # Step 3: Merge ligand GRO coordinates
#         lig_gro = std_files.get("gro_file")
#         if lig_gro and os.path.exists(lig_gro):
#             print(f"Merging ligand coordinates from {lig_gro}")
            
#             # Get ligand atom count and coordinates
#             with open(lig_gro, 'r') as f:
#                 lines = f.readlines()
#                 lig_atom_count = int(lines[1].strip())
#                 lig_coords = ''.join(lines[2:-1])  # All lines between header+count and box
            
#             # Get protein atom count and coordinates
#             with open(os.path.join(model_dir, "pro_lig.gro"), 'r') as f:
#                 lines = f.readlines()
#                 protein_atom_count = int(lines[1].strip())
#                 gro_header = lines[0]
#                 gro_box = lines[-1]
#                 protein_coords = ''.join(lines[2:-1])
            
#             # Calculate total atom count
#             total_atom_count = protein_atom_count + lig_atom_count
            
#             # Write merged GRO file
#             with open(os.path.join(model_dir, "pro_lig_merged.gro"), 'w') as f:
#                 f.write(gro_header)
#                 f.write(f"{total_atom_count}\n")
#                 f.write(protein_coords)
#                 f.write(lig_coords)
#                 f.write(gro_box)
            
#             # Backup original and replace with merged
#             shutil.move(
#                 os.path.join(model_dir, "pro_lig.gro"),
#                 os.path.join(model_dir, "pro_lig_backup.gro")
#             )
#             shutil.move(
#                 os.path.join(model_dir, "pro_lig_merged.gro"),
#                 os.path.join(model_dir, "pro_lig.gro")
#             )
            
#             print("Ligand coordinates successfully merged into pro_lig.gro")
#         else:
#             print(f"Error: Ligand GRO file not found at {lig_gro}")
#             return None
        
#         # Step 4: Modify topology to include ligand parameters
#         topol_path = os.path.join(model_dir, "topol.top")
#         if os.path.exists(topol_path) and "atomtypes_file" in std_files:
#             print("Updating topology file to include ligand parameters")
            
#             # Read the atomtypes block
#             with open(std_files["atomtypes_file"], 'r') as f:
#                 atomtypes_block = f.read()
            
#             # Read the topology file
#             with open(topol_path, 'r') as f:
#                 topol_content = f.read()
            
#             # Insert atomtypes and ligand include after forcefield include
#             force_field_include_pattern = r'(#include "amber99sb.ff/forcefield.itp")'
#             modified_content = re.sub(
#                 force_field_include_pattern,
#                 f'\\1\n\n{atomtypes_block}\n\n#include "../LIG.amb2gmx/LIG.itp"',
#                 topol_content
#             )
            
#             # Add LIG to molecules section if not already there
#             if "[ molecules ]" in modified_content and "LIG" not in modified_content.split("[ molecules ]")[1]:
#                 molecules_section_pattern = r'(\[ *molecules *\].*?)$'
#                 modified_content = re.sub(
#                     molecules_section_pattern,
#                     f'\\1\nLIG\t1',
#                     modified_content,
#                     flags=re.DOTALL
#                 )
            
#             # Write modified topology
#             with open(topol_path, 'w') as f:
#                 f.write(modified_content)
            
#             print("Topology updated successfully")
#         else:
#             if not os.path.exists(topol_path):
#                 print(f"Error: Topology file not found at {topol_path}")
#             if "atomtypes_file" not in std_files:
#                 print(f"Error: Atomtypes file missing from std_files")
#             return None
        
#         # Step 5: Set up the system
#         print("Setting up the system (box, solvent, ions)...")
        
#         # Create box
#         self.run_command([
#             "gmx", "editconf",
#             "-f", "pro_lig.gro",
#             "-o", "pro_lig_newbox.gro",
#             "-c",
#             "-box", "10", "10", "10"
#         ])
        
#         # Solvate
#         self.run_command([
#             "gmx", "solvate",
#             "-cp", "pro_lig_newbox.gro",
#             "-p", "topol.top",
#             "-o", "solv.gro"
#         ])
        
#         # 创建ions.mdp如果不存在
#         ions_mdp = os.path.join(model_dir, "ions.mdp")
#         if not os.path.exists(ions_mdp):
#             with open(ions_mdp, 'w') as f:
#                 f.write("""title                   = Ions
# ; Run parameters
# integrator              = steep
# nsteps                  = 5000
# ; energy minimization
# emtol                   = 1000.0
# emstep                  = 0.01
# nstcomm                 = 100
# ; Output control
# nstxout                 = 500
# nstvout                 = 500
# nstenergy               = 500
# nstlog                  = 500
# ; Neighbor searching
# cutoff-scheme           = Verlet
# ns_type                 = grid
# nstlist                 = 10
# rcoulomb                = 1.0
# rvdw                    = 1.0
# ; Electrostatics
# coulombtype             = PME
# pme_order               = 4
# fourierspacing          = 0.16
# ; Constraints
# constraints             = h-bonds
# constraint-algorithm    = LINCS
# """)
        
#         # Generate ions.tpr
#         self.run_command([
#             "gmx", "grompp",
#             "-f", ions_mdp,
#             "-c", "solv.gro",
#             "-p", "topol.top",
#             "-o", "ions.tpr",
#             "-maxwarn", "10"
#         ])
        
#         # 使用Popen直接向genion提供输入
#         try:
#             print("Running genion with direct input...")
#             # 直接使用Popen运行genion并提供输入
#             self.run_command(
#                 ["gmx", "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top", 
#                 "-pname", "NA", "-nname", "CL", "-neutral"],
#                 input_text="15"  # SOL group
#             )
#         except Exception as e:
#             print(f"Error running genion: {e}")
#             print("Trying alternative method with echo...")
            
#             # 尝试使用echo和管道
#             genion_cmd = "echo 15 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral"
#             self.run_command(genion_cmd, shell=True)
        
#         os.chdir(current_dir)
#         print(f"Model build completed. System files are in {model_dir}")
#         return model_dir

#     def cleanup(self):
#         """Clean up temporary files"""
#         print("Cleaning up temporary files...")
        
#         # Remove Amber temporary files
#         temp_files = [
#             "sqm*", "ANTECHAMBER*", "PREP*", "ATOMTYPE*", "NEWPDB*", "leap*"
#         ]
        
#         for pattern in temp_files:
#             for file_path in Path(self.output_dir).glob(pattern):
#                 try:
#                     os.remove(file_path)
#                     print(f"Removed {file_path}")
#                 except OSError as e:
#                     print(f"Failed to remove {file_path}: {e}")
        
#         # Clean forcefield directory
#         ff_dir = os.path.join(self.output_dir, "forcefield")
#         if os.path.exists(ff_dir):
#             for ext in ["*frcmod", "*prep", "*prmtop", "*rst7", "*log", "*in"]:
#                 for file_path in Path(ff_dir).glob(ext):
#                     try:
#                         os.remove(file_path)
#                         print(f"Removed {file_path}")
#                     except OSError as e:
#                         print(f"Failed to remove {file_path}: {e}")
                        
#         print("Cleanup completed")

#     def run(self):
#         """
#         Run the entire GAFFMaker workflow
        
#         Returns:
#         --------
#         str
#             Path to the output directory
#         """
#         print(f"Starting GAFFMaker workflow...")
        
#         try:
#             # Step 1: Clean protein
#             cleaned_protein = self.clean_protein()
            
#             # Step 2: Generate ligand forcefield
#             ff_files = self.generate_ligand_forcefield()
            
#             # Step 3: Standardize files
#             std_files = self.standardize_files(ff_files)
            
#             # Step 4: Build model
#             model_dir = self.build_model(cleaned_protein, std_files)
#             if not model_dir:
#                 raise RuntimeError("Failed to build model. See previous errors.")
            
#             # Step 5: Cleanup
#             self.cleanup()
            
#             print(f"\nGAFFMaker workflow completed successfully!")
#             print(f"Output files are in {self.output_dir}")
#             print(f"You can now run MD simulations using the files in {os.path.join(self.output_dir, 'GMX_PROLIG_MD')}")
            
#             return self.output_dir
            
#         except Exception as e:
#             print(f"Error during GAFFMaker workflow: {e}")
#             import traceback
#             traceback.print_exc()
#             raise


# def main():
#     """Main function to parse arguments and run GAFFMaker"""
#     parser = argparse.ArgumentParser(description="GAFFMaker: Build protein-ligand systems for GROMACS")
    
#     parser.add_argument("protein", help="Path to protein PDB file")
#     parser.add_argument("ligand", help="Path to ligand MOL2 file")
#     parser.add_argument("--output", "-o", default="output", help="Output directory")
    
#     args = parser.parse_args()
    
#     # Create and run GAFFMaker
#     gaffmaker = GAFFMaker(args.protein, args.ligand, args.output)
#     gaffmaker.run()


# if __name__ == "__main__":
#     main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GAFFMaker - A tool for building protein-ligand systems for GROMACS MD simulations

This script streamlines the process of modeling protein-ligand molecular dynamics simulations
with GROMACS by automating the preparation of force field parameters using AmberTools
and setting up the simulation system.

Usage:
    python GAFF_Maker.py <protein_path> <ligand_path> --output <output_dir>

Author: Adapted from GAFFMakerV1.1 bash scripts
"""

import os
import sys
import subprocess
import shutil
import re
import tempfile
import argparse
import glob
from pathlib import Path


class GAFFMaker:
    def __init__(self, protein_path, ligand_path, output_dir, overwrite=False):
        """
        Initialize the GAFFMaker with protein and ligand paths
        
        Parameters:
        -----------
        protein_path : str
            Path to the protein PDB file
        ligand_path : str
            Path to the ligand MOL2 file
        output_dir : str
            Directory where output files will be stored
        overwrite : bool
            Whether to overwrite existing files
        """
        self.protein_path = os.path.abspath(protein_path)
        self.ligand_path = os.path.abspath(ligand_path)
        self.output_dir = os.path.abspath(output_dir)
        self.overwrite = overwrite
        
        # Extract filenames without extensions
        self.protein_name = os.path.basename(protein_path).split('.')[0]
        self.ligand_name = os.path.basename(ligand_path).split('.')[0]
        self.ligand_id = self._extract_ligand_id()
        
        # Initialize list of possible molecule names (will be expanded later)
        self.possible_molecule_names = [
            self.ligand_id,           # ID from MOL2 file (e.g., 1unl_ligand)
            self.ligand_name,         # Filename (e.g., ligand)
            self.ligand_id.upper(),   
            self.ligand_name.upper(),
            "LIG"                     # Standard name
        ]

        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Directory for MDP files
        self.mdps_dir = os.path.join(self.output_dir, "mdps")
        os.makedirs(self.mdps_dir, exist_ok=True)
        
        print(f"Initialized GAFFMaker for protein: {self.protein_path}")
        print(f"                       and ligand: {self.ligand_path}")
        print(f"                   Output directory: {self.output_dir}")
        
        # Debug level
        self.debug = True

    def _extract_ligand_id(self):
        """Extract ligand ID from MOL2 file"""
        try:
            with open(self.ligand_path, 'r') as f:
                for line in f:
                    if line.startswith("@<TRIPOS>MOLECULE"):
                        # Read the next line which contains the molecule name
                        ligand_id = next(f).strip()
                        print(f"Extracted ligand ID from mol2 file: {ligand_id}")
                        return ligand_id
        except Exception as e:
            print(f"Failed to extract ligand ID from MOL2 file: {e}")
        
        # Fallback to filename
        return self.ligand_name

    def run_command(self, cmd, cwd=None, env=None, check=True, input_text=None, shell=False):
        """Run a shell command and return its output"""
        print(f"Running command: {cmd}")
        
        if isinstance(cmd, str):
            cmd_str = cmd
            if not shell:
                cmd_list = cmd.split()
        else:
            cmd_list = cmd
            cmd_str = " ".join(cmd)
            
        try:
            # Use Popen for commands that need input via pipe
            if input_text is not None and ("genion" in cmd_str or "pdb2gmx" in cmd_str):
                print(f"Using Popen with input: {input_text}")
                process = subprocess.Popen(
                    cmd_str if shell else cmd_list,
                    cwd=cwd,
                    env=env,
                    shell=shell,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                stdout, stderr = process.communicate(input=input_text)
                
                if process.returncode != 0 and check:
                    print(f"Error executing command: {cmd}")
                    print(f"Return code: {process.returncode}")
                    print(f"Error output: {stderr}")
                    raise subprocess.CalledProcessError(process.returncode, cmd, stdout, stderr)
                
                print(f"Command completed with return code: {process.returncode}")
                if stdout.strip():
                    print(f"Output snippet: {stdout.strip()[:100]}...")
                return stdout
            else:
                # Use run for other commands
                kwargs = {
                    'cwd': cwd,
                    'env': env,
                    'check': check,
                    'text': True,
                    'stdout': subprocess.PIPE,
                    'stderr': subprocess.PIPE,
                    'shell': shell
                }
                
                # Add input only for non-special commands
                if input_text is not None:
                    kwargs['input'] = input_text
                    
                result = subprocess.run(
                    cmd_str if shell else cmd_list,
                    **kwargs
                )
                
                print(f"Command completed with return code: {result.returncode}")
                if result.stdout.strip():
                    print(f"Output snippet: {result.stdout.strip()[:100]}...")
                return result.stdout
                
        except subprocess.CalledProcessError as e:
            print(f"Error executing command: {cmd}")
            print(f"Return code: {e.returncode}")
            print(f"Error output: {e.stderr}")
            if check:
                raise
            return e.stdout

    def clean_protein(self):
        """Clean the protein PDB file by removing HETATM lines"""
        print("Cleaning protein PDB file...")
        
        # Create path for cleaned PDB
        cleaned_pdb = os.path.join(self.output_dir, f"{self.protein_name}_clean.pdb")
        
        # Check if cleaned file already exists
        if os.path.exists(cleaned_pdb) and not self.overwrite:
            print(f"Cleaned protein file already exists at {cleaned_pdb}, skipping cleaning step")
            return cleaned_pdb
        
        # Read PDB file and remove HETATM lines
        with open(self.protein_path, 'r') as f_in, open(cleaned_pdb, 'w') as f_out:
            for line in f_in:
                if not line.startswith('HETATM'):
                    f_out.write(line)
        
        # Fix HN1/HN2/HN3 atom names
        with open(cleaned_pdb, 'r') as f:
            content = f.read()
        
        content = content.replace('HN1', 'H1 ').replace('HN2', 'H2 ').replace('HN3', 'H3 ')
        
        with open(cleaned_pdb, 'w') as f:
            f.write(content)
        
        print(f"Protein cleaned and saved to {cleaned_pdb}")
        return cleaned_pdb

    def generate_ligand_forcefield(self):
        """Generate force field parameters for the ligand"""
        print("Generating ligand force field parameters...")
        
        # Create forcefield directory
        ff_dir = os.path.join(self.output_dir, "forcefield")
        os.makedirs(ff_dir, exist_ok=True)
        
        # Check if already processed
        acpype_check = [d for d in os.listdir(ff_dir) 
                      if (d.endswith('.acpype') or d.endswith('.amb2gmx')) and os.path.isdir(os.path.join(ff_dir, d))]
        
        if acpype_check and not self.overwrite:
            print(f"Found existing acpype output in {ff_dir}, skipping force field generation")
            # Still need to find the files
            return self._find_forcefield_files(ff_dir)
        
        # Paths for intermediate files
        amber_mol2 = os.path.join(ff_dir, f"{self.ligand_name}.mol2")
        prep_file = os.path.join(ff_dir, f"{self.ligand_name}.prep")
        frcmod_file = os.path.join(ff_dir, f"{self.ligand_name}.frcmod")
        prmtop_file = os.path.join(ff_dir, f"{self.ligand_name}.prmtop")
        rst7_file = os.path.join(ff_dir, f"{self.ligand_name}.rst7")
        
        # Step 1: Generate Amber-compatible mol2
        self.run_command([
            "antechamber", 
            "-i", self.ligand_path, 
            "-fi", "mol2", 
            "-o", amber_mol2, 
            "-fo", "mol2", 
            "-c", "bcc", 
            "-s", "2"
        ])
        
        # Step 2: Generate prep file
        self.run_command([
            "antechamber", 
            "-i", amber_mol2, 
            "-fi", "mol2", 
            "-o", prep_file, 
            "-fo", "prepi", 
            "-c", "bcc", 
            "-s", "2"
        ])
        
        # Step 3: Generate frcmod
        self.run_command([
            "parmchk2", 
            "-i", prep_file, 
            "-f", "prepi", 
            "-o", frcmod_file
        ])
        
        # Step 4: Prepare tleap input
        tleap_input = os.path.join(ff_dir, f"tleap_{self.ligand_name}.in")
        with open(tleap_input, 'w') as f:
            f.write(f"""source leaprc.protein.ff14SB
source leaprc.gaff

LIG = loadmol2 {amber_mol2}
loadamberparams {frcmod_file}

saveamberparm LIG {prmtop_file} {rst7_file}

quit
""")
        
        # Step 5: Run tleap
        self.run_command([
            "tleap", 
            "-f", tleap_input
        ])
        
        # Step 6: Convert to GROMACS format with acpype
        current_dir = os.getcwd()
        os.chdir(ff_dir)
        
        if self.debug:
            print(f"Before acpype, current directory: {os.getcwd()}")
            print(f"Directory contents: {os.listdir('.')}")
        
        self.run_command([
            "acpype", 
            "-p", os.path.basename(prmtop_file), 
            "-x", os.path.basename(rst7_file), 
            "-d"
        ])
        
        if self.debug:
            print(f"After acpype, current directory: {os.getcwd()}")
            print(f"Directory contents: {os.listdir('.')}")
        
        os.chdir(current_dir)
        
        # Find the generated files
        return self._find_forcefield_files(ff_dir)
    
    def _find_forcefield_files(self, ff_dir):
        """Find force field files generated by acpype"""
        # Try multiple possible acpype output directories based on ligand name/ID
        possible_acpype_dirs = [
            os.path.join(ff_dir, f"{self.ligand_name}.acpype"),
            os.path.join(ff_dir, f"{self.ligand_id}.acpype"),
            os.path.join(ff_dir, f"{self.ligand_name.upper()}.acpype"),
            os.path.join(ff_dir, f"{self.ligand_id.upper()}.acpype"),
            os.path.join(ff_dir, f"{self.ligand_name.upper()}.amb2gmx"),
            os.path.join(ff_dir, f"{self.ligand_id.upper()}.amb2gmx")
        ]
        
        # Find any additional *.acpype or *.amb2gmx directories
        for item in os.listdir(ff_dir):
            item_path = os.path.join(ff_dir, item)
            if os.path.isdir(item_path) and (item.endswith('.acpype') or item.endswith('.amb2gmx')):
                if item_path not in possible_acpype_dirs:
                    possible_acpype_dirs.append(item_path)
                    # Extract potential molecule name for future replacements
                    mol_name = item.split('.')[0]  # Remove .acpype or .amb2gmx
                    if mol_name and mol_name not in self.possible_molecule_names:
                        self.possible_molecule_names.append(mol_name)
                        print(f"Added potential molecule name from directory: {mol_name}")
        
        acpype_dir = None
        for dir_path in possible_acpype_dirs:
            if os.path.exists(dir_path):
                acpype_dir = dir_path
                print(f"Found acpype output directory: {acpype_dir}")
                
                # Extract molecule name for substitution
                mol_name = os.path.basename(dir_path).split('.')[0]
                if mol_name and mol_name not in self.possible_molecule_names:
                    self.possible_molecule_names.append(mol_name)
                    print(f"Added molecule name from acpype directory: {mol_name}")
                
                if self.debug:
                    print(f"Contents of acpype directory: {os.listdir(acpype_dir)}")
                break
                
        if not acpype_dir:
            print("Warning: Could not find acpype output directory in expected locations")
            # Search for any .acpype/.amb2gmx directory
            possible_dirs = [d for d in os.listdir(ff_dir) if (d.endswith('.acpype') or d.endswith('.amb2gmx')) and os.path.isdir(os.path.join(ff_dir, d))]
            if possible_dirs:
                acpype_dir = os.path.join(ff_dir, possible_dirs[0])
                print(f"Using alternative acpype directory: {acpype_dir}")
                
                # Extract molecule name for substitution
                mol_name = possible_dirs[0].split('.')[0]
                if mol_name and mol_name not in self.possible_molecule_names:
                    self.possible_molecule_names.append(mol_name)
                    print(f"Added molecule name from alternative directory: {mol_name}")
                
                if self.debug:
                    print(f"Contents of alternative acpype directory: {os.listdir(acpype_dir)}")
        
        # Find force field files
        itp_file = None
        gro_file = None
        top_file = None
        posre_file = None
        
        if acpype_dir and os.path.exists(acpype_dir):
            # Look for files in acpype directory
            for filename in os.listdir(acpype_dir):
                filepath = os.path.join(acpype_dir, filename)
                if filename.endswith('_GMX.itp'):
                    itp_file = filepath
                    # Extract residue name from ITP file
                    try:
                        mol_name = os.path.basename(filename).split('_GMX.itp')[0]
                        if mol_name and mol_name not in self.possible_molecule_names:
                            self.possible_molecule_names.append(mol_name)
                            print(f"Added molecule name from ITP file: {mol_name}")
                    except:
                        pass
                    print(f"Found ITP file: {filepath}")
                elif filename.endswith('_GMX.gro'):
                    gro_file = filepath
                    print(f"Found GRO file: {filepath}")
                elif filename.endswith('_GMX.top'):
                    top_file = filepath
                    print(f"Found TOP file: {filepath}")
                elif filename.startswith('posre_') and filename.endswith('.itp'):
                    posre_file = filepath
                    # Extract residue name from posre file
                    try:
                        mol_name = os.path.basename(filename).replace('posre_', '').replace('.itp', '')
                        if mol_name and mol_name not in self.possible_molecule_names:
                            self.possible_molecule_names.append(mol_name)
                            print(f"Added molecule name from posre file: {mol_name}")
                    except:
                        pass
                    print(f"Found POSRE file: {filepath}")
        
        # If files not found, perform recursive search
        if not all([gro_file, top_file]):
            print("Warning: Some required files not found. Performing recursive search...")
            for root, dirs, files in os.walk(ff_dir):
                for file in files:
                    if file.endswith('_GMX.itp') and not itp_file:
                        itp_file = os.path.join(root, file)
                        # Extract residue name
                        try:
                            mol_name = file.split('_GMX.itp')[0]
                            if mol_name and mol_name not in self.possible_molecule_names:
                                self.possible_molecule_names.append(mol_name)
                                print(f"Added molecule name from recursive search ITP: {mol_name}")
                        except:
                            pass
                        print(f"Found ITP file in recursive search: {itp_file}")
                    elif file.endswith('_GMX.gro') and not gro_file:
                        gro_file = os.path.join(root, file)
                        print(f"Found GRO file in recursive search: {gro_file}")
                    elif file.endswith('_GMX.top') and not top_file:
                        top_file = os.path.join(root, file)
                        print(f"Found TOP file in recursive search: {top_file}")
                    elif file.startswith('posre_') and file.endswith('.itp') and not posre_file:
                        posre_file = os.path.join(root, file)
                        # Extract residue name
                        try:
                            mol_name = file.replace('posre_', '').replace('.itp', '')
                            if mol_name and mol_name not in self.possible_molecule_names:
                                self.possible_molecule_names.append(mol_name)
                                print(f"Added molecule name from recursive search posre: {mol_name}")
                        except:
                            pass
                        print(f"Found POSRE file in recursive search: {posre_file}")
        
        print(f"Final molecule name list for replacements: {self.possible_molecule_names}")
        
        ff_files = {
            "itp_file": itp_file,
            "gro_file": gro_file,
            "top_file": top_file,
            "posre_file": posre_file
        }
        
        print(f"Force field files found: {ff_files}")
        return ff_files

    def standardize_files(self, ff_files):
        """Standardize file names and directory structure"""
        print("Standardizing files...")
        
        # Create LIG.amb2gmx directory
        lig_dir = os.path.join(self.output_dir, "LIG.amb2gmx")
        
        # Check if standardization already done
        if os.path.exists(lig_dir) and os.path.exists(os.path.join(lig_dir, "LIG.itp")) and not self.overwrite:
            print(f"Standardized files already exist in {lig_dir}, skipping standardization")
            # Return the paths to the standardized files
            return {
                "gro_file": os.path.join(lig_dir, "LIG.gro") if os.path.exists(os.path.join(lig_dir, "LIG.gro")) else None,
                "itp_file": os.path.join(lig_dir, "LIG.itp") if os.path.exists(os.path.join(lig_dir, "LIG.itp")) else None,
                "top_file": os.path.join(lig_dir, "LIG.top") if os.path.exists(os.path.join(lig_dir, "LIG.top")) else None,
                "posre_file": os.path.join(lig_dir, "posre_LIG.itp") if os.path.exists(os.path.join(lig_dir, "posre_LIG.itp")) else None,
                "atomtypes_file": os.path.join(lig_dir, "atomtypes_block.itp") if os.path.exists(os.path.join(lig_dir, "atomtypes_block.itp")) else None
            }
        
        os.makedirs(lig_dir, exist_ok=True)
        std_files = {}
        
        # Process GRO file
        if ff_files["gro_file"] and os.path.exists(ff_files["gro_file"]):
            std_files["gro_file"] = os.path.join(lig_dir, "LIG.gro")
            shutil.copy(ff_files["gro_file"], std_files["gro_file"])
            print(f"Copied {ff_files['gro_file']} to {std_files['gro_file']}")
            
            # Replace all possible molecule names with LIG
            for mol_name in self.possible_molecule_names:
                if mol_name != "LIG":
                    self._replace_in_file(std_files["gro_file"], mol_name, "LIG")
            print(f"Replaced molecule names with 'LIG' in {std_files['gro_file']}")
        else:
            print(f"Warning: GRO file not found, skipping. Path was: {ff_files['gro_file']}")
        
        # Process TOP file
        if ff_files["top_file"] and os.path.exists(ff_files["top_file"]):
            std_files["top_file"] = os.path.join(lig_dir, "LIG.top")
            shutil.copy(ff_files["top_file"], std_files["top_file"])
            print(f"Copied {ff_files['top_file']} to {std_files['top_file']}")
            
            # Replace all possible molecule names with LIG
            for mol_name in self.possible_molecule_names:
                if mol_name != "LIG":
                    self._replace_in_file(std_files["top_file"], mol_name, "LIG")
            print(f"Replaced molecule names with 'LIG' in {std_files['top_file']}")
        else:
            print(f"Warning: TOP file not found, skipping. Path was: {ff_files['top_file']}")
        
        # Create or process ITP file
        if ff_files["itp_file"] and os.path.exists(ff_files["itp_file"]):
            # Process existing ITP file
            std_files["itp_file"] = os.path.join(lig_dir, "LIG.itp")
            self._create_itp_from_existing(ff_files["itp_file"], std_files["itp_file"])
            print(f"Created {std_files['itp_file']} from existing {ff_files['itp_file']}")
        elif ff_files["top_file"] and os.path.exists(ff_files["top_file"]):
            # Create ITP from TOP file
            std_files["itp_file"] = os.path.join(lig_dir, "LIG.itp")
            self._create_itp_from_top(ff_files["top_file"], std_files["itp_file"])
            print(f"Created {std_files['itp_file']} from TOP file {ff_files['top_file']}")
        else:
            print(f"Warning: Cannot create ITP file - no source file available")
        
        # Process POSRE file
        if ff_files["posre_file"] and os.path.exists(ff_files["posre_file"]):
            std_files["posre_file"] = os.path.join(lig_dir, "posre_LIG.itp")
            shutil.copy(ff_files["posre_file"], std_files["posre_file"])
            
            # Replace all possible molecule names with LIG
            for mol_name in self.possible_molecule_names:
                if mol_name != "LIG":
                    self._replace_in_file(std_files["posre_file"], mol_name, "LIG")
            print(f"Created {std_files['posre_file']} from {ff_files['posre_file']} and replaced molecule names with 'LIG'")
        else:
            print(f"Warning: POSRE file not found, skipping. Path was: {ff_files['posre_file']}")
        
        # Extract atomtypes block from the TOP file
        if ff_files["top_file"] and os.path.exists(ff_files["top_file"]):
            with open(ff_files["top_file"], 'r') as f:
                content = f.read()
            
            # Use more flexible regex to find atomtypes section
            atomtypes_match = re.search(r'\[ *atomtypes *\](.*?)(\[ *[a-z]|$)', content, re.DOTALL)
            if atomtypes_match:
                atomtypes_block = "[ atomtypes ]\n" + atomtypes_match.group(1).strip()
                std_files["atomtypes_file"] = os.path.join(lig_dir, "atomtypes_block.itp")
                with open(std_files["atomtypes_file"], 'w') as f:
                    f.write(atomtypes_block)
                print(f"Created atomtypes file {std_files['atomtypes_file']}")
            else:
                print(f"Warning: Could not extract atomtypes block from TOP file: {ff_files['top_file']}")
                # Create an empty atomtypes file as fallback
                std_files["atomtypes_file"] = os.path.join(lig_dir, "atomtypes_block.itp")
                with open(std_files["atomtypes_file"], 'w') as f:
                    f.write("[ atomtypes ]\n; Empty atomtypes file created because extraction failed\n")
                print(f"Created empty atomtypes file as fallback")
        
        # Check for missing required files
        required_files = ["gro_file", "itp_file", "atomtypes_file"]
        missing_files = [f for f in required_files if f not in std_files]
        
        if missing_files:
            print(f"WARNING: The following required files are missing: {missing_files}")
            print("This may cause problems in model building!")
            
        print(f"Standardized files created in {lig_dir}")
        return std_files
    
    def _create_itp_from_existing(self, source_itp, target_itp):
        """Create standardized ITP file from existing ITP file"""
        with open(source_itp, 'r') as f_in, open(target_itp, 'w') as f_out:
            content = f_in.read()
            
            # Replace all possible molecule names with LIG
            for mol_name in self.possible_molecule_names:
                if mol_name != "LIG":
                    content = re.sub(r'\b{}\b'.format(re.escape(mol_name)), "LIG", content)
            
            # Ensure posre include is present at the end
            if "#ifdef POSRES" not in content:
                content += "\n#ifdef POSRES\n"
                content += '#include "posre_LIG.itp"\n'
                content += "#endif\n"
            
            f_out.write(content)
    
    def _create_itp_from_top(self, top_file, itp_file):
        """Create ITP file by extracting moleculetype section from TOP file"""
        try:
            with open(top_file, 'r') as f:
                content = f.read()
            
            # Find moleculetype section up to system section (or end of file)
            mol_match = re.search(r'\[ *moleculetype *\](.*?)(\[ *system|$)', content, re.DOTALL)
            
            if mol_match:
                # Extract moleculetype section
                mol_content = "[ moleculetype ]\n" + mol_match.group(1).strip()
                
                # Replace all possible molecule names with LIG
                for mol_name in self.possible_molecule_names:
                    if mol_name != "LIG":
                        mol_content = re.sub(r'\b{}\b'.format(re.escape(mol_name)), "LIG", mol_content)
                
                # Add posre include
                mol_content += "\n\n#ifdef POSRES\n"
                mol_content += '#include "posre_LIG.itp"\n'
                mol_content += "#endif\n"
                
                with open(itp_file, 'w') as f:
                    f.write(mol_content)
                
                return True
            else:
                # If moleculetype section not found, create a basic placeholder
                print(f"Warning: Could not find moleculetype section in {top_file}")
                basic_itp = """[ moleculetype ]
; Name            nrexcl
LIG              3

[ atoms ]
; Note: This is a placeholder. The actual ITP content could not be extracted.
; Please check your ligand structure and force field parameters.

#ifdef POSRES
#include "posre_LIG.itp"
#endif
"""
                with open(itp_file, 'w') as f:
                    f.write(basic_itp)
                
                return False
                
        except Exception as e:
            print(f"Error creating ITP from TOP file: {e}")
            return False
    
    def _replace_in_file(self, file_path, old_str, new_str):
        """Replace text in a file"""
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Use regex to replace complete words only
            updated_content = re.sub(r'\b{}\b'.format(re.escape(old_str)), new_str, content)
            
            # Only write back if content changed
            if updated_content != content:
                with open(file_path, 'w') as f:
                    f.write(updated_content)
                return True
            return False
        except Exception as e:
            print(f"Error replacing text in {file_path}: {e}")
            return False

    def build_model(self, cleaned_protein, std_files):
        """Build the GROMACS model"""
        print("Building GROMACS model...")
        
        # Create model directory
        model_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
        
        # Check if model building already completed
        if os.path.exists(os.path.join(model_dir, "solv_ions.gro")) and not self.overwrite:
            print(f"Final model file solv_ions.gro already exists in {model_dir}, skipping model building")
            return model_dir
            
        os.makedirs(model_dir, exist_ok=True)
        
        # Check for required files
        required_files = {
            "gro_file": "Ligand GRO file",
            "itp_file": "Ligand ITP file", 
            "atomtypes_file": "Atomtypes file"
        }
        
        missing_files = []
        for file_key, file_desc in required_files.items():
            if file_key not in std_files or not std_files[file_key] or not os.path.exists(std_files[file_key]):
                missing_files.append(f"{file_desc} ({file_key})")
                
        if missing_files:
            error_msg = "Missing required files: " + ", ".join(missing_files)
            print(f"ERROR: {error_msg}")
            return None
        
        # Step 1: Run pdbfixer to clean the protein
        fixed_pdb = os.path.join(model_dir, "fixed_clean_protein.pdb")
        if not os.path.exists(fixed_pdb) or self.overwrite:
            self.run_command([
                "pdbfixer", 
                cleaned_protein, 
                "--output", fixed_pdb, 
                "--add-atoms", "heavy", 
                "--keep-heterogens", "none"
            ])
        else:
            print(f"Fixed protein PDB already exists at {fixed_pdb}, skipping pdbfixer")
        
        # Step 2: Generate topology using gmx pdb2gmx
        current_dir = os.getcwd()
        os.chdir(model_dir)
        
        pro_lig_gro = os.path.join(model_dir, "pro_lig.gro")
        if not os.path.exists(pro_lig_gro) or self.overwrite:
            # Run pdb2gmx with input for AMBER99SB and TIP3P water
            self.run_command(
                ["gmx", "pdb2gmx", "-f", "fixed_clean_protein.pdb", "-o", "pro_lig.gro", "-ignh"],
                input_text="6\n1\n"
            )
        else:
            print(f"Protein topology already generated at {pro_lig_gro}, skipping pdb2gmx")
        
        # Step 3: Merge ligand GRO coordinates
        merged_gro = os.path.join(model_dir, "pro_lig_merged.gro")
        if (not os.path.exists(merged_gro) and not os.path.exists(os.path.join(model_dir, "pro_lig_backup.gro"))) or self.overwrite:
            lig_gro = std_files.get("gro_file")
            if lig_gro and os.path.exists(lig_gro):
                print(f"Merging ligand coordinates from {lig_gro}")
                
                # Get ligand atom count and coordinates
                with open(lig_gro, 'r') as f:
                    lines = f.readlines()
                    lig_atom_count = int(lines[1].strip())
                    lig_coords = ''.join(lines[2:-1])
                
                # Get protein atom count and coordinates
                with open(pro_lig_gro, 'r') as f:
                    lines = f.readlines()
                    protein_atom_count = int(lines[1].strip())
                    gro_header = lines[0]
                    gro_box = lines[-1]
                    protein_coords = ''.join(lines[2:-1])
                
                # Calculate total atom count
                total_atom_count = protein_atom_count + lig_atom_count
                
                # Write merged GRO file
                with open(merged_gro, 'w') as f:
                    f.write(gro_header)
                    f.write(f"{total_atom_count}\n")
                    f.write(protein_coords)
                    f.write(lig_coords)
                    f.write(gro_box)
                
                # Backup original and replace with merged
                shutil.move(pro_lig_gro, os.path.join(model_dir, "pro_lig_backup.gro"))
                shutil.move(merged_gro, pro_lig_gro)
                
                print("Ligand coordinates successfully merged into pro_lig.gro")
            else:
                print(f"Error: Ligand GRO file not found at {lig_gro}")
                return None
        else:
            print(f"Ligand coordinates already merged, skipping merging step")
        
        # Step 4: Modify topology to include ligand parameters
        topol_path = os.path.join(model_dir, "topol.top")
        if os.path.exists(topol_path) and "atomtypes_file" in std_files:
            if "LIG.amb2gmx/LIG.itp" not in open(topol_path).read() or self.overwrite:
                print("Updating topology file to include ligand parameters")
                
                # Read the atomtypes block
                with open(std_files["atomtypes_file"], 'r') as f:
                    atomtypes_block = f.read()
                
                # Read the topology file
                with open(topol_path, 'r') as f:
                    topol_content = f.read()
                
                # Insert atomtypes and ligand include after forcefield include
                force_field_include_pattern = r'(#include "amber99sb.ff/forcefield.itp")'
                modified_content = re.sub(
                    force_field_include_pattern,
                    f'\\1\n\n{atomtypes_block}\n\n#include "../LIG.amb2gmx/LIG.itp"',
                    topol_content
                )
                
                # Add LIG to molecules section if not already there
                if "[ molecules ]" in modified_content and "LIG" not in modified_content.split("[ molecules ]")[1]:
                    molecules_section_pattern = r'(\[ *molecules *\].*?)$'
                    modified_content = re.sub(
                        molecules_section_pattern,
                        f'\\1\nLIG\t1',
                        modified_content,
                        flags=re.DOTALL
                    )
                
                # Write modified topology
                with open(topol_path, 'w') as f:
                    f.write(modified_content)
                
                print("Topology updated successfully")
            else:
                print(f"Topology already includes ligand parameters, skipping update")
        else:
            if not os.path.exists(topol_path):
                print(f"Error: Topology file not found at {topol_path}")
            if "atomtypes_file" not in std_files:
                print(f"Error: Atomtypes file missing from std_files")
            return None
        
        # Step 5: Set up the system (box, solvent, ions)
        print("Setting up the system (box, solvent, ions)...")
        
        # Create box
        box_gro = os.path.join(model_dir, "pro_lig_newbox.gro")
        if not os.path.exists(box_gro) or self.overwrite:
            self.run_command([
                "gmx", "editconf",
                "-f", "pro_lig.gro",
                "-o", "pro_lig_newbox.gro",
                "-c",
                "-box", "10", "10", "10"
            ])
        else:
            print(f"Box already created at {box_gro}, skipping")
        
        # Solvate
        solv_gro = os.path.join(model_dir, "solv.gro")
        if not os.path.exists(solv_gro) or self.overwrite:
            self.run_command([
                "gmx", "solvate",
                "-cp", "pro_lig_newbox.gro",
                "-p", "topol.top",
                "-o", "solv.gro"
            ])
        else:
            print(f"System already solvated at {solv_gro}, skipping")
        
        # Create ions.mdp if it doesn't exist
        ions_mdp = os.path.join(model_dir, "ions.mdp")
        if not os.path.exists(ions_mdp) or self.overwrite:
            with open(ions_mdp, 'w') as f:
                f.write("""title                   = Ions
; Run parameters
integrator              = steep
nsteps                  = 5000
; energy minimization
emtol                   = 1000.0
emstep                  = 0.01
nstcomm                 = 100
; Output control
nstxout                 = 500
nstvout                 = 500
nstenergy               = 500
nstlog                  = 500
; Neighbor searching
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
; Electrostatics
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
; Constraints
constraints             = h-bonds
constraint-algorithm    = LINCS
""")
        
        # Generate ions.tpr
        ions_tpr = os.path.join(model_dir, "ions.tpr")
        if not os.path.exists(ions_tpr) or self.overwrite:
            self.run_command([
                "gmx", "grompp",
                "-f", ions_mdp,
                "-c", "solv.gro",
                "-p", "topol.top",
                "-o", "ions.tpr",
                "-maxwarn", "10"
            ])
        else:
            print(f"Ions TPR already generated at {ions_tpr}, skipping")
        
        # Add ions
        solv_ions_gro = os.path.join(model_dir, "solv_ions.gro")
        if not os.path.exists(solv_ions_gro) or self.overwrite:
            try:
                print("Running genion with direct input...")
                # Use Popen to provide input directly to genion
                self.run_command(
                    ["gmx", "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top", 
                    "-pname", "NA", "-nname", "CL", "-neutral"],
                    input_text="15"  # SOL group
                )
            except Exception as e:
                print(f"Error running genion: {e}")
                print("Trying alternative method with echo...")
                
                # Fallback to shell method with pipe
                genion_cmd = "echo 15 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral"
                self.run_command(genion_cmd, shell=True)
        else:
            print(f"System already ionized at {solv_ions_gro}, skipping")
        
        os.chdir(current_dir)
        print(f"Model build completed. System files are in {model_dir}")
        return model_dir

    def cleanup(self):
        """Clean up temporary files"""
        print("Cleaning up temporary files...")
        
        # Remove Amber temporary files
        temp_files = [
            "sqm*", "ANTECHAMBER*", "PREP*", "ATOMTYPE*", "NEWPDB*", "leap*"
        ]
        
        for pattern in temp_files:
            for file_path in Path(self.output_dir).glob(pattern):
                try:
                    os.remove(file_path)
                    print(f"Removed {file_path}")
                except OSError as e:
                    print(f"Failed to remove {file_path}: {e}")
        
        # Clean forcefield directory
        ff_dir = os.path.join(self.output_dir, "forcefield")
        if os.path.exists(ff_dir):
            for ext in ["*frcmod", "*prep", "*prmtop", "*rst7", "*log", "*in"]:
                for file_path in Path(ff_dir).glob(ext):
                    try:
                        os.remove(file_path)
                        print(f"Removed {file_path}")
                    except OSError as e:
                        print(f"Failed to remove {file_path}: {e}")
                        
        print("Cleanup completed")

    def run(self):
        """Run the entire GAFFMaker workflow"""
        print(f"Starting GAFFMaker workflow...")
        
        try:
            # Step 1: Clean protein
            cleaned_protein = self.clean_protein()
            
            # Step 2: Generate ligand forcefield
            ff_files = self.generate_ligand_forcefield()
            
            # Step 3: Standardize files
            std_files = self.standardize_files(ff_files)
            
            # Step 4: Build model
            model_dir = self.build_model(cleaned_protein, std_files)
            if not model_dir:
                raise RuntimeError("Failed to build model. See previous errors.")
            
            # Step 5: Cleanup
            self.cleanup()
            
            print(f"\nGAFFMaker workflow completed successfully!")
            print(f"Output files are in {self.output_dir}")
            print(f"You can now run MD simulations using the files in {os.path.join(self.output_dir, 'GMX_PROLIG_MD')}")
            
            return self.output_dir
            
        except Exception as e:
            print(f"Error during GAFFMaker workflow: {e}")
            import traceback
            traceback.print_exc()
            raise


def main():
    """Main function to parse arguments and run GAFFMaker"""
    parser = argparse.ArgumentParser(description="GAFFMaker: Build protein-ligand systems for GROMACS")
    
    parser.add_argument("protein", help="Path to protein PDB file")
    parser.add_argument("ligand", help="Path to ligand MOL2 file")
    parser.add_argument("--output", "-o", default="output", help="Output directory")
    parser.add_argument("--overwrite", "-f", action="store_true", help="Overwrite existing files")
    
    args = parser.parse_args()
    
    # Create and run GAFFMaker
    gaffmaker = GAFFMaker(args.protein, args.ligand, args.output, args.overwrite)
    gaffmaker.run()


if __name__ == "__main__":
    main()