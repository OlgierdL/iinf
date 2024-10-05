import subprocess
import math
import Bio
import os
import glob
import platform
from Bio.PDB.PDBParser import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import shutil
import sys
import argparse
from itertools import permutations
import tempfile
from pathlib import Path
import csv
from itertools import product
import re


def clear(file):
    output_file = file[0:len(file) - 4] + "_tmp_m" + file[len(file) - 4:]
    with open(file, "r") as infile, open(output_file, "w") as outfile:
        no = None
        a = 0
        for line in infile:
            fields = line.strip().split()
            if fields[0] not in ["PFRMAT", "TARGET", "PARENT"] and fields[0] != "MODEL" and fields[0] != "END":
                outfile.write(line)
            elif fields[0] == "MODEL":
                if no != fields[1]:
                    outfile.write(line)
                    no = fields[1]
            elif fields[0] == "TER":
                a = 1
            elif fields[0] == "END":
                if a != 1:
                    outfile.write("TER\n")
                a = 0
        outfile.write("END\n")


def analyze(file):
    output = {}
    dna_dict = ["DT", "DA", "DC", "DG"]
    rna_dict = ["A", "C", "G", "U"]
    protein_dict = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                    "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    Select = Bio.PDB.Select
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure("str", file)

    distribution = {}
    molecule_chains = {}
    for model in structure:
        distribution[model.id] = {}
        molecule_chains[model.id] = {"RNA": [], "DNA": [], "Protein": []}
        for chain in model:
            distribution[model.id][chain.id] = {"RNA": 0, "DNA": 0, "Protein": 0}
            for residue in chain:
                res = residue.get_resname().strip()
                if res in protein_dict:
                    distribution[model.id][chain.id]["Protein"] = distribution[model.id][chain.id].pop("Protein") + 1
                elif res in dna_dict:
                    distribution[model.id][chain.id]["DNA"] = distribution[model.id][chain.id].pop("DNA") + 1
                elif res in rna_dict:
                    distribution[model.id][chain.id]["RNA"] = distribution[model.id][chain.id].pop("RNA") + 1
            molecule_type = max(distribution[model.id][chain.id], key=distribution[model.id][chain.id].get)
            molecule_chains[model.id][molecule_type].append(chain.id)
    for model in structure:
        output["model_no"] = format(model.serial_num)
        output["DNA"] = format(",".join(molecule_chains[model.id]["DNA"]))
        output["RNA"] = format(",".join(molecule_chains[model.id]["RNA"]))
        output["Protein"] = format(",".join(molecule_chains[model.id]["Protein"]))
    for model in structure:
        resnum_tot = 0
        for chain in model:
            residues = list(chain.get_residues())
            residues_no = len(residues)
            (field1, resseq1, icode1) = residues[0].id
            (field2, resseq2, icode2) = residues[residues_no - 1].id
            residues = [chain.id, resseq1, resseq2, residues_no]
            output[chain.id] = residues
            resnum_tot += residues_no
        output["residues_no"] = resnum_tot
    if (output["Protein"] == " " and output["RNA"] != "T"):
        output["Protein"] = "T"
    elif (output["Protein"] == " "):
        output["Protein"] = "S"
    if (len(output["Protein"].split(",")) == 1):
        output["single_protein"] = True
    else:
        output["single_protein"] = False
    return output


def separateChains(file, chain_ids, single_chain):
    with open("tmp/" + chain_ids[0] + chain_ids[1] + ".pdb", "w") as chains:
        file.seek(0)
        for line in file:
            line_tmp = " ".join(line.split())
            if (len(line_tmp.split()) >= 5):
                line_id = line_tmp.split()[4]
                if ((line_id in chain_ids or single_chain) and line_tmp.split()[0] == "ATOM"): chains.write(line)


def rmsd(native_coords, model_coords, rot, tran):
    model_coords_rotated = np.dot(model_coords, rot) + tran
    diff = native_coords - model_coords_rotated
    RMSD = np.sqrt(sum(sum(diff ** 2)) / native_coords.shape[0])
    return RMSD


def specific_align(aln_atoms, target, model):
    atom_types = ["CA", "N", "C", "O"]
    AA = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
          "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    p = PDBParser(QUIET=True)
    native = p.get_structure("target", target)
    model = p.get_structure("model", model)

    native_coords = [a.coord for a in native[0].get_atoms() if a.parent.resname in AA and a.name in atom_types]
    model_coords = [a.coord for a in model[0].get_atoms() if a.parent.resname in AA and a.name in atom_types]

    native_coords = np.array(native_coords)
    model_coords = np.array(model_coords)

    percentage_to_aln = int(aln_atoms * len(native_coords))

    si = SVDSuperimposer()
    si.set(native_coords[:percentage_to_aln], model_coords[:percentage_to_aln])
    si.run()

    if (native_coords.shape == model_coords.shape):
        RMSD = rmsd(native_coords, model_coords, si.rot, si.tran)
    else:
        RMSD = -1
    return [si, RMSD]


def get_mapping(targetData, modelData, target, model):
    for id in targetData["Protein"].split(","):
        separateChains(target, [targetData["RNA"], id], targetData["single_protein"])

    for id in modelData["Protein"].split(","):
        separateChains(model, [modelData["RNA"], id], modelData["single_protein"])

    rmsds = {}
    for id1 in targetData["Protein"].split(","):
        for id2 in modelData["Protein"].split(","):
            rmsd_value = \
                specific_align(1, "tmp/" + targetData["RNA"] + id1 + ".pdb", "tmp/" + modelData["RNA"] + id2 + ".pdb")[
                    1]
            rmsds[id1 + ":" + id2] = rmsd_value

    rmsds = dict(sorted(rmsds.items(), key=lambda x: x[1]))
    matched = []
    mapping = []
    for key in rmsds.keys():
        if (key[0] not in [match[0] for match in mapping] and key[2] not in [match[2] for match in mapping]):
            mapping.append(key)
    chains_mapping_target = ""
    for pair in mapping:
        if (len(chains_mapping_target) > 0): chains_mapping_target += ","
        chains_mapping_target += pair
    return chains_mapping_target


def run_hbplus(tmpdir, name1):
    command = f'bash -i -c "find {tmpdir} -name \'*.pdb\' -exec hbplus \\{{\\}} \\;"'
    subprocess.run(command, shell=True, cwd=tmpdir, stdout=subprocess.DEVNULL)
    name3 = name1[0:len(name1) - 3] + "hb2"
    target_HB2 = open(name3, "r")
    hb2_dict = {f.name: 0.0 for f in Path(tmpdir).glob("*.hb2") if f.stem != Path(name3).stem}
    return hb2_dict, target_HB2


def getid(don):
    return don[0] + " " + don[1:5] + " " + don[5]


def setDict(t):
    dict = {}
    for i in t:
        dict[i] = 1
    return dict


def remove_leading_0s(s: str) -> str:
    return s.lstrip("0") or "0"


def filter_pairs(file, rna_chains, protein_chains):
    a = []
    dna_str = ["DT", "DA", "DC", "DG"]
    rna_str = ["A", "C", "G", "U"]
    protein_str = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                   "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    dna_dict = setDict(dna_str)
    rna_dict = setDict(rna_str)
    protein_dict = setDict(protein_str)
    i = 0
    file.seek(0)
    for line in file:
        line.strip()
        if (len(line) != 76 or not line.split()[-1].isdigit()):
            pass
        else:
            donor = line[:9]
            acceptor = line[14:23]

            chain_d = donor[0].strip()

            chain_a = acceptor[0].strip()

            name_d = donor[6:9].strip()
            name_a = acceptor[6:9].strip()

            if (name_d[0] != "HOH" and name_a != "HOH") and \
                    (len(rna_chains) == 0 or (len(rna_chains) > 0 and (
                            rna_chains.find(chain_d) >= 0 or rna_chains.find(chain_a) >= 0))) and \
                    (name_d in rna_dict.keys() or name_a in rna_dict.keys()) and \
                    not ((chain_d in protein_chains.split(",") and chain_a in protein_chains.split(",")) or (
                            name_d in dna_dict.keys() and name_a in dna_dict.keys() > 0) or (
                                 chain_d in rna_chains.split(",") and chain_a in rna_chains.split(","))):
                don_id = getid(donor)
                acc_id = getid(acceptor)
                don_tmp1 = don_id[0:2]
                don_tmp2 = remove_leading_0s(don_id[2:9])
                don_id = don_tmp1 + don_tmp2
                acc_tmp1 = acc_id[0:2]
                acc_tmp2 = remove_leading_0s(acc_id[2:9])
                acc_id = acc_tmp1 + acc_tmp2
                don_acc_id = don_id + " " + acc_id
                if (don_acc_id not in a):
                    a.append(don_acc_id)
                    i += 1
    return a


def replace_chains(line, mapping):
    x = []
    y = []
    for swap in mapping:
        z = swap.split(":")
        x.append(z[0])
        y.append(z[1])

    fields = line.split()

    for k in range(len(x)):
        if fields[0] == x[k]:
            fields[0] = y[k]
            break

    for k in range(len(x)):
        if fields[3] == x[k]:
            fields[3] = y[k]
            break
    return (" ".join(fields))


def get_pairs(file, rna_chains, protein_chains, mapping):
    pairs = filter_pairs(file, rna_chains, protein_chains)
    sorted_pairs = sorted(pairs, key=lambda x: (x.split()[0], int(x.split()[1])))
    replaced_pairs = []
    for pair in sorted_pairs:
        if mapping != None:
            replaced_pairs.append(replace_chains(pair, mapping.split(',')))
        else:
            replaced_pairs.append(pair)
    for pair in replaced_pairs:
        parts = pair.split()
        a = parts[0] + parts[1]
        if parts[2] != "-":
            a += parts[2]
        b = parts[3] + parts[4]
        if parts[5] != "-":
            b += parts[5]
    return replaced_pairs


def get_inf(target, model):
    target_dict = {}
    i = 0
    for pair in target:
        i += 1
        target_dict[pair] = i

    model_dict = {}
    j = 0
    for pair in model:
        j += 1
        model_dict[pair] = j

    tp = 0
    fn = 0

    for k in range(i):
        if (target[k] in model_dict):
            tp += 1
        else:
            fn += 1
    fp = 0
    for k in range(j):
        if (not (model[k] in target_dict)):
            fp += 1

    ppv = tp / (tp + fp)
    tpr = tp / (tp + fn)
    inf = math.sqrt(ppv * tpr)

    return inf


def convert_to_wsl_path(path):
    if platform.system() == "Windows":
        return subprocess.check_output(["wsl", "wslpath", path.replace("\\", "/")]).decode().strip()
    else:
        return path


def rna_tools_renumerate(filename, output_file, edit_command):
    command = f'rna_pdb_tools.py --keep-hetatm --edit \'{edit_command}\' {filename} > {output_file}'
    subprocess.run(["bash", "-c", command])


def rna_tools_delete(filename, edit_command):
    for chain in edit_command.split(","):
        command = f'rna_pdb_tools.py --delete \'{chain}\' --inplace {filename}'
        subprocess.run(["bash", "-c", command])


def delete_residues(filename, edit_command):
    for chain_command in edit_command.split(","):
        if (len(chain_command.split(".")) == 1):
            chain = chain_command[0]
            start, end = chain_command.split('-')
            end = int(end)
            start = int(start[2:])
        else:
            chain = chain_command[0]
            start, end = chain_command.split('.')
            end = int(end)
            start = int(start[2:])
        with open(filename, 'r') as file:
            lines = file.readlines()

        with open(filename, 'w') as file:
            for line in lines:
                if len(line) >= 32:
                    if line[20:23].strip() == chain:
                        num_str = line[23:32].strip()
                        try:
                            num = int(num_str)
                        except ValueError:
                            continue
                        if start <= num <= end:
                            continue
                file.write(line)

def renumber_residues(filename, output_file, edit_command):
    for chain_command in edit_command.split(","):
        original, renumber = chain_command.split('>')
        if(len(original.split(".")) == 1 and len(original.split(".")) == 1):
            og_chain_id = original[0]
            ren_chain_id = renumber[0]
            og_start, og_end = original.split('-')
            og_end = int(og_end)
            og_start = int(og_start[2:])
            ren_start, ren_end = renumber.split('-')
            ren_start = int(ren_start[2:])
        else:
            og_chain_id = original[0]
            ren_chain_id = renumber[0]
            og_start, og_end = original.split('.')
            og_end = int(og_end)
            og_start = int(og_start[2:])
            ren_start, ren_end = renumber.split('.')
            ren_start = int(ren_start[2:])

        with open(filename, 'r') as file:
            lines = file.readlines()
        with open(filename, 'r') as file:
            file.seek(0)
            renum_index = ren_start
            last_num = og_start
            last_letter = ""
            surplus = 0
            for i, line in enumerate(lines):
                if len(line) >= 32:
                    if line[20:23].strip() == og_chain_id:
                        num_str = line[23:31].strip()
                        match = re.match(r'^(-?[0-9]+)([a-zA-Z]*)$', num_str)
                        if match:
                            num = match.group(1).strip()
                            letters = match.group(2).strip()
                        else: continue
                        try:
                            num = int(num)
                        except ValueError:
                            continue
                        if og_start <= num <= og_end + surplus:
                            if (last_num != num or last_letter != letters):
                                renum_index += 1
                                if(last_letter != letters):
                                    surplus += 1
                            spaces1 = str(((" ") * (4 - len(str(renum_index)))))
                            spaces2 = (" ") * (4)
                            new_line = line[:21] + ren_chain_id + spaces1 + str(renum_index) + spaces2 + line[30:]
                            last_num = num
                            last_letter = letters
                            lines[i] = new_line
            file.seek(0)
        with open(filename, 'w') as file:
            file.writelines(lines)



def find_next_identifier(identifiers):
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    identifiers_set = set(identifiers)
    for letter in alphabet:
        if letter not in identifiers_set:
            return letter
    for letters in product(alphabet, repeat=2):
        identifier = ''.join(letters)
        if identifier not in identifiers_set:
            return identifier


def add_protein_chain_identifier(input_file, output_file, rna_id):
    identifier = "A" if rna_id == "0" else find_next_identifier(rna_id)

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if len(line) > 21:
                id = line[21]
                if id == " ":
                    line = line[:21] + identifier + line[22:]
            outfile.write(line)


def auto_renumber(filename, chains, output):
    command = []
    for chain in chains:
        i = 1
        fragments, negatives= find_fragments(filename, chain)
        chain_mapping = []
        if negatives: sign = "."
        else: sign = "-"
        for fragment in fragments:
            fragment_mapping = chain + ":" + str(fragment[0]) + sign + str(fragment[1]) + ">" + chain + ":" + str(
                i) + sign + str(fragment[1] - fragment[0] + i)
            chain_mapping.append(fragment_mapping)
            i = fragment[1] - fragment[0] + i + 1
        command.append(",".join(chain_mapping))
    edit_command = ",".join(command)
    renumber_residues(filename, output, edit_command)


def find_fragments(filename, chain):
    numbers = []
    negatives_found = False
    with open(filename) as file:
        for line in file:
            if (len(line) > 40 and line[21] == chain):
                number = re.search(r'-?\d+', line[22:30].strip()).group()
                if(number.startswith("-")): negatives_found = True
                if (len(numbers) == 0):
                    i = number
                    numbers.append(int(number))
                if (number != i): numbers.append(int(number))
                i = number
        if (len(numbers) > 0): fragments = [[numbers[0]]]
        for i in range(1, len(numbers)):
            if (numbers[i] > numbers[i - 1] + 1):
                fragments[len(fragments) - 1].append(numbers[i - 1])
                fragments.append([numbers[i]])
        fragments[len(fragments) - 1].append(numbers[len(numbers) - 1])
        return fragments, negatives_found


def create_combinations(list1, list2):
    target_mappings = []
    for combination in permutations(list2, len(list1)):
        combined_str = ",".join(f"{a}:{b}" for a, b in zip(list1, combination))
        target_mappings.append(combined_str)
    return target_mappings


def copy_file_to_script_directory(source_path):
    script_directory = os.path.dirname(os.path.abspath(__file__))

    file_name = os.path.basename(source_path)

    destination_path = os.path.join(script_directory, file_name)

    if not os.path.exists(destination_path):
        try:
            shutil.copy2(source_path, destination_path)
        except FileNotFoundError:
            print(f"Source file not found: {source_path}")
        except PermissionError:
            print(f"Permission denied while copying the file.")
        except Exception as e:
            print(f"An error occurred: {e}")
    else:
        return False
    return True


def save_csv(result_file_path, data):
    with open(result_file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)


def copy_to_tmp(tmpdir, names):
    names_copied = []
    for name in names:
        name_copied = shutil.copy(name, tmpdir)
        names_copied.append(os.path.basename(name_copied))
    return names_copied


def single_chain_rename(name2, modelData, model):
    output_file = name2[0:len(name2) - 4] + "_id" + name2[len(name2) - 4:]
    add_protein_chain_identifier(name2, output_file, modelData["RNA"])
    model.close()
    os.remove(name2)
    new_name = name2[0:len(name2) - 4] + "_id" + name2[len(name2) - 4:]
    return new_name


def alpha_rename(name2, modelData, model):
    id = find_next_identifier(modelData["Protein"])
    rep = modelData["RNA"] + ">" + id
    output = name2[0:len(name2) - 4] + "_t" + name2[len(name2) - 4:]
    command = f'rna_pdb_tools.py --rename-chain \'{rep}\' {name2} > {output}'
    subprocess.run(["bash", "-c", command])
    model.close()
    os.remove(name2)
    return output


def back_rename(modelData, name2, model):
    rep = modelData["RNA"] + ">" + "0"
    output = name2[0:len(name2) - 4] + "_t" + name2[len(name2) - 4:]
    command = f'rna_pdb_tools.py --rename-chain \'{rep}\' {name2} > {output}'
    subprocess.run(["bash", "-c", command])
    model.close()
    os.remove(name2)
    return output


def custom_delete(done, name, delete):
    if (not done):
        delete_residues(name, delete)


def custom_renum(done, name, file, renum):
    if (not done):
        renumber_residues(name, name[0:len(name) - 4] + "_cus" + name[len(name) - 4:], renum)


def remove_om(file_path, rna_chains):
    modified_letters_3to1 = {'A23': 'A', 'A2L': 'A', 'A2M': 'A', 'A39': 'A',
                             'A3P': 'A', 'A44': 'A', 'A5O': 'A', 'A6A': 'A', 'A7E': 'A', 'A9Z': 'A',
                             'ADI': 'A', 'ADP': 'A', 'AET': 'A', 'AMD': 'A', 'AMO': 'A', 'AP7': 'A',
                             'AVC': 'A', 'MA6': 'A', 'MAD': 'A', 'MGQ': 'A', 'MIA': 'A', 'MTU': 'A',
                             'M7A': 'A', '26A': 'A', '2MA': 'A', '6IA': 'A', '6MA': 'A', '6MC': 'A',
                             '6MP': 'A', '6MT': 'A', '6MZ': 'A', '6NW': 'A', 'F3N': 'A', 'N79': 'A',
                             'RIA': 'A', 'V3L': 'A', 'ZAD': 'A', '31H': 'A', '31M': 'A', '7AT': 'A',
                             'O2Z': 'A', 'SRA': 'A', '00A': 'A', '45A': 'A', '8AN': 'A', 'LCA': 'A',
                             'P5P': 'A', 'PPU': 'A', 'PR5': 'A', 'PU': 'A', 'T6A': 'A', 'TBN': 'A',
                             'TXD': 'A', 'TXP': 'A', '12A': 'A', '1MA': 'A', '5FA': 'A', 'A6G': 'G',
                             'E6G': 'G', 'E7G': 'G', 'EQ4': 'G', 'IG': 'G', 'IMP': 'G', 'M2G': 'G',
                             'MGT': 'G', 'MGV': 'G', 'MHG': 'G', 'QUO': 'G', 'YG': 'G', 'YYG': 'G',
                             '23G': 'G', '2EG': 'G', '2MG': 'G', '2SG': 'G', 'B8K': 'G', 'B8W': 'G',
                             'B9B': 'G', 'BGH': 'G', 'N6G': 'G', 'RFJ': 'G', 'ZGU': 'G', '7MG': 'G',
                             'CG1': 'G', 'G1G': 'G', 'G25': 'G', 'G2L': 'G', 'G46': 'G', 'G48': 'G',
                             'G7M': 'G', 'GAO': 'G', 'GDO': 'G', 'GDP': 'G', 'GH3': 'G', 'GNG': 'G',
                             'GOM': 'G', 'GRB': 'G', 'GTP': 'G', 'KAG': 'G', 'KAK': 'G', 'O2G': 'G',
                             'OMG': 'G', '8AA': 'G', '8OS': 'G', 'LG': 'G', 'PGP': 'G', 'P7G': 'G',
                             'TPG': 'G', 'TG': 'G', 'XTS': 'G', '102': 'G', '18M': 'G', '1MG': 'G',
                             'A5M': 'C', 'A6C': 'C', 'E3C': 'C', 'IC': 'C', 'M4C': 'C', 'M5M': 'C',
                             '6OO': 'C', 'B8Q': 'C', 'B8T': 'C', 'B9H': 'C', 'JMH': 'C', 'N5M': 'C',
                             'RPC': 'C', 'RSP': 'C', 'RSQ': 'C', 'ZBC': 'C', 'ZCY': 'C', '73W': 'C',
                             'C25': 'C', 'C2L': 'C', 'C31': 'C', 'C43': 'C', 'C5L': 'C', 'CBV': 'C',
                             'CCC': 'C', 'CH': 'C', 'CSF': 'C', 'OMC': 'C', 'S4C': 'C', '4OC': 'C',
                             'LC': 'C', 'LHH': 'C', 'LV2': 'C', 'PMT': 'C', 'TC': 'C', '10C': 'C',
                             '1SC': 'C', '5HM': 'C', '5IC': 'C', '5MC': 'C', 'A6U': 'U', 'IU': 'U',
                             'I4U': 'U', 'MEP': 'U', 'MNU': 'U', 'U25': 'U', 'U2L': 'U', 'U2P': 'U',
                             'U31': 'U', 'U34': 'U', 'U36': 'U', 'U37': 'U', 'U8U': 'U', 'UAR': 'U',
                             'UBB': 'U', 'UBD': 'U', 'UD5': 'U', 'UPV': 'U', 'UR3': 'U', 'URD': 'U',
                             'US5': 'U', 'UZR': 'U', 'UMO': 'U', 'U23': 'U', '2AU': 'U', '2MU': 'U',
                             '2OM': 'U', 'B8H': 'U', 'FHU': 'U', 'FNU': 'U', 'F2T': 'U', 'RUS': 'U',
                             'ZBU': 'U', '3AU': 'U', '3ME': 'U', '3MU': 'U', '3TD': 'U', '70U': 'U',
                             '75B': 'U', 'CNU': 'U', 'OMU': 'U', 'ONE': 'U', 'S4U': 'U', 'SSU': 'U',
                             'SUR': 'U', '4SU': 'U', '85Y': 'U', 'DHU': 'U', 'H2U': 'U', 'LHU': 'U',
                             'PSU': 'U', 'PYO': 'U', 'P4U': 'U', 'T31': 'U', '125': 'U', '126': 'U',
                             '127': 'U', '1RN': 'U', '5BU': 'U', '5FU': 'U', '5MU': 'U', '9QV': 'U',
                             '5GP': 'G', 'ATP': 'A'}
    chains = rna_chains.split(",")
    invalid_lines = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for i in range(len(lines)):
        if len(lines[i]) >= 19:
            substring = lines[i][17:19]
            if lines[i][21] in chains and substring != "  ":
                if lines[i][17:20] in modified_letters_3to1.keys():
                    lines[i] = lines[i][:17] + "  " + modified_letters_3to1[lines[i][17:20]] + lines[i][20:]
                else:
                     print("Unrecognized nucleotide: " + substring + " in chain " + lines[i][21] + " in structure " + file_path.split("_")[0].split("/")[len(file_path.split("_")[0].split("/")) - 1] + " " + str(i))
                     invalid_lines.append(lines[i])
    for line in invalid_lines: lines.remove(line)
    with open(file_path, 'w') as file:
        file.writelines(lines)


def copy_to_script_dir(file_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_name = os.path.basename(file_path)
    dest_path = os.path.join(script_dir, file_name)
    shutil.copy(file_path, dest_path)

    print(f"File '{file_name}' copied to {script_dir}")


def check_for_negative(data):
    for id in data["RNA"].split(","):
        if (data[id][1] < 0): return True
    for id in data["Protein"].split(","):
        if (data[id][1] < 0): return True
    return False

def auto_renum(data, done, name, file):
    to_renum = []
    renumber = False

    for id in data["Protein"].split(","):
        if (data[id][2] - data[id][1] + 1 != data[id][3] or data[id][1] != 1 and not done):
            renumber = True
            to_renum.append(data[id][0])

    for id in data["RNA"].split(","):
        if (data[id][2] - data[id][1] + 1 != data[id][3] or data[id][1] != 1 and not done):
            renumber = True
            to_renum.append(data[id][0])

    if (renumber):
        auto_renumber(name, to_renum, name[0:len(name) - 4] + "_ren" + name[len(name) - 4:])


def get_chains(target_mappings, chains_mapping_model, is_target):
    result = ''
    for mapping in target_mappings:
        target_chains = mapping.split(':')
        if len(result) > 0:
            result = result + ','
        if (is_target):
            result = result + target_chains[0]
        else:
            result = result + target_chains[1]
    if is_target:
        model_chains = chains_mapping_model.split(':')[0].split(',')
    else:
        model_chains = chains_mapping_model.split(':')[1].split(',')
    for mapping in model_chains:
        if len(result) > 0:
            result = result + ','
        result = result + mapping
    return result


def get_mapping(is_target, model_chains_protein, target_chains_protein, model_chains_rna, target_chains_rna):
    result = ''
    for i in range(len(target_chains_protein)):
        if(len(model_chains_protein) > i):
            if (target_chains_protein[i] != model_chains_protein[i]):
                if len(result) > 0:
                    result = result + ','
                if (is_target):
                    result = result + target_chains_protein[i] + ':' + model_chains_protein[i]
                else:
                    result = result + model_chains_protein[i] + ':' + target_chains_protein[i]
    for i in range(len(target_chains_rna)):
        if(len(model_chains_rna) > i):
            if (target_chains_rna[i] != model_chains_rna[i]):
                if len(result) > 0:
                    result = result + ','
                if (is_target):
                    result = result + target_chains_rna[i] + ':' + model_chains_rna[i]
                else:
                    result = result + model_chains_rna[i] + ':' + target_chains_rna[i]
    if len(result) > 0:
        return result
    return None


def get_max_inf(target_mappings, chains_mapping_model, target_HB2, model_HB2, modelData, targetData):
    max_inf = 0
    model_chains = modelData["RNA"] + "," + modelData["Protein"]
    target_chains = targetData["RNA"] + "," + targetData["Protein"]
    model_mapping = get_mapping(False, modelData["Protein"], targetData["Protein"], modelData["RNA"], targetData["RNA"])
    target_pairs = get_pairs(target_HB2, targetData["RNA"], targetData["Protein"], None)
    model_pairs = get_pairs(model_HB2, modelData["RNA"], modelData["Protein"],  model_mapping)
    if (len(model_pairs) == 0):
        pass
    elif (len(target_pairs) == 0):
        pass
    else:
        inf = get_inf(target_pairs, model_pairs)
        if (inf > max_inf):
            max_inf = inf
    return max_inf


def shorten_for_output(input_string):
    tmp_index = input_string.find("_tmp")

    if tmp_index != -1:
        part_until_tmp = input_string[:tmp_index]
    else:
        part_until_tmp = input_string

    return part_until_tmp


def compare(name1, names2, custom_alignement, adj_inf, renumber, target_renum, model_renum, delete, target_delete,
            model_delete):
    infs = []
    target_done = False

    with tempfile.TemporaryDirectory() as tmpdir:
        name1_copied = shutil.copy(name1, tmpdir)
        name1_basename = os.path.basename(name1_copied)
        name1 = ""
        name2 = ""

        names_copied = copy_to_tmp(tmpdir, names2)
        for name in names_copied:
            name2 = os.path.join(tmpdir, name)
            if (not target_done):
                name1 = os.path.join(tmpdir, name1_basename)
            target = open(name1, "r")
            model = open(name2, "r")
            clear(name2)
            model.close()
            os.remove(name2)
            name2 = name2[0:len(name2) - 4] + "_tmp_m" + name2[len(name2) - 4:]
            model = open(name2, "r")

            if (not target_done): targetData = analyze(name1)
            modelData = analyze(name2)


            if (modelData["single_protein"] == True):
                name2 = single_chain_rename(name2, modelData, model)
                model = open(name2, "r")
                modelData = analyze(name2)

            if (not target_done and targetData["single_protein"] == True):
                name1 = single_chain_rename(name1, targetData, target)
                target = open(name1, "r")
                targetData = analyze(name1)

            remove_om(name1, targetData["RNA"])
            remove_om(name2, modelData["RNA"])

            if (delete):
                if (target_delete != ""): custom_delete(target_done, name1, target_delete); targetData = analyze(name1)
                if (name[0:len(name) - 4] in model_delete.keys()): custom_delete(False, name2, model_delete[
                    name[0:len(name) - 4]]); modelData = analyze(name2)

            if (renumber):
                if (custom_alignement):
                    if (target_renum != ""): custom_renum(target_done, name1, target,
                                                                          target_renum); targetData = analyze(name1)
                    if (name[0:len(name) - 4] in model_renum.keys()): custom_renum(False, name2, model,
                                                                                                  model_renum[name[
                                                                                                              0:len(
                                                                                                                  name) - 4]]); modelData = analyze(
                        name2)
                else:
                    auto_renum(targetData, target_done, name1, target)
                    auto_renum(modelData, False, name2, model)
                    targetData = analyze(name1)
                    modelData = analyze(name2)

            if(not target_done):
                if check_for_negative(targetData):
                    print("Negative residue numbers in target. Please renumber. You may use -r.")
                    exit()
            if check_for_negative(modelData):
                print("Negative residue numbers in " + name + " Please renumber. You may use -r.")
                exit()

            if (not target_done): target_done = True

        hb2_dict, target_HB2 = run_hbplus(tmpdir, name1)

        for model in hb2_dict.keys():
            model_HB2 = open(os.path.join(tmpdir, model), "r")

            os.makedirs("tmp_chains", exist_ok=True)

            chains_mapping_model = modelData["RNA"] + ":" + targetData["RNA"]

            target_mappings = create_combinations(targetData["Protein"].split(","), modelData["Protein"].split(","))

            inf = get_max_inf(target_mappings, chains_mapping_model, target_HB2, model_HB2, modelData, targetData)

            model_HB2.close()

            if (adj_inf):
                inf = inf * modelData["residues_no"] / targetData["residues_no"]
                hb2_dict[model] = inf
                model_name = shorten_for_output(model)
                infs.append([model_name, format(hb2_dict[model], ".3f")])
            else:
                hb2_dict[model] = inf
                model_name = shorten_for_output(model)
                infs.append([model_name, format(hb2_dict[model], ".3f")])

            tmp_files = glob.glob(os.path.join("tmp_chains", "*"))
            for tmp_file in tmp_files:
                pass
                try:
                    os.remove(tmp_file)
                except Exception as e:
                    print(f'Error deleting {tmp_file}: {e}')

            if not os.listdir("tmp_chains"):
                os.rmdir("tmp_chains")
    infs.sort()
    infs = [["model", "score"]] + infs
    return infs


def main(argv):
    parser = argparse.ArgumentParser(description="Compare model and target.")
    parser.add_argument("--target_path", type=str, help="target", required=True)
    parser.add_argument("--model_path", type=str, help="model", required=True)

    parser.add_argument("-a", "--adjust_inf", action="store_true", help="Adjust inf")
    parser.add_argument("-r", "--renumber_structures", action="store_true", help="Renumber chains")

    parser.add_argument("-c", "--custom_alignement", type=str, help="Custom renumbering", default=None)
    parser.add_argument("-d", "--custom_removal", type=str, help="Custom deleting", default=None)

    args = parser.parse_args()
    if not os.path.isfile(args.target_path):
        print(f"Error: Target path '{args.target_path}' does not exist or is not a file.")
        sys.exit(1)

    if not os.path.exists(args.model_path):
        print(f"Error: Model path '{args.model_path}' does not exist.")
        sys.exit(1)

    if os.path.isdir(args.model_path):
        files_to_compare = [os.path.join(args.model_path, f) for f in os.listdir(args.model_path) if f.endswith('.pdb')]
        if not files_to_compare:
            print(f"Error: No PDB files found in directory '{args.model_path}'.")
            sys.exit(1)
    else:
        files_to_compare = [args.model_path]

    custom_renumbering = False
    custom_target_renum = ""
    custom_model_renum = {}

    if (not args.custom_alignement is None):
        custom_renumbering = True
        custom_renums = args.custom_alignement.split(";")
        for renum in custom_renums:
            name, renum = renum.split("|")
            if any(os.path.splitext(os.path.basename(path))[0] == name for path in files_to_compare):
                custom_model_renum[name] = renum
            if (name == os.path.splitext(os.path.basename(args.target_path))[0]): custom_target_renum = renum

    custom_delete = False
    custom_target_delete = ""
    custom_model_delete = {}

    if (not args.custom_removal is None):
        custom_delete = True
        custom_deletes = args.custom_removal.split(";")
        for delete in custom_deletes:
            name, delete = delete.split("|")
            if any(os.path.splitext(os.path.basename(path))[0] == name for path in files_to_compare):
                custom_model_delete[name] = delete
            if (name == os.path.splitext(os.path.basename(args.target_path))[0]): custom_target_delete = delete

    infs = compare(args.target_path, files_to_compare, custom_renumbering, args.adjust_inf,
                   args.renumber_structures, custom_target_renum, custom_model_renum, custom_delete,
                   custom_target_delete, custom_model_delete)
    target_filename_without_ext = os.path.splitext(os.path.basename(args.target_path))[0]
    sorted_infs = [infs[0]] + sorted(infs[1:], key=lambda x: x[1], reverse=True)
    save_csv(os.path.join(os.path.dirname(args.target_path), '{}_ranking.csv'.format(target_filename_without_ext)),
             sorted_infs)


if __name__ == "__main__":
    main(sys.argv[1:])