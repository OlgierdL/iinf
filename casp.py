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

def clear(file):
    output_file = file[0:len(file)-4] + "m" + file[len(file)-4:]
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
    if(output["Protein"] == " " and output["RNA"] != "T"): output["Protein"] = "T"
    elif(output["Protein"] == " "): output["Protein"] = "S"
    if (len(output["Protein"].split(",")) == 1): output["single_protein"] = True
    else: output["single_protein"] = False
    return output


def separateChains(file, chain_ids, single_chain):
    with open("tmp/" + chain_ids[0] + chain_ids[1] + ".pdb", "w") as chains:
        file.seek(0)
        for line in file:
            line_tmp = " ".join(line.split())
            if(len(line_tmp.split()) >= 5):
                line_id = line_tmp.split()[4]
                if((line_id in chain_ids or single_chain) and line_tmp.split()[0] == "ATOM"): chains.write(line)


def rmsd(native_coords, model_coords, rot, tran):
    model_coords_rotated = np.dot(model_coords, rot) + tran
    diff = native_coords - model_coords_rotated
    RMSD = np.sqrt(sum(sum(diff**2))/native_coords.shape[0])
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

    if(native_coords.shape == model_coords.shape): RMSD = rmsd(native_coords, model_coords, si.rot, si.tran)
    else: RMSD = -1
    return [si, RMSD]


def get_mapping(targetData, modelData, target, model):
    for id in targetData["Protein"].split(","):
        separateChains(target, [targetData["RNA"], id], targetData["single_protein"])

    for id in modelData["Protein"].split(","):
        separateChains(model, [modelData["RNA"], id], modelData["single_protein"])

    rmsds = {}
    for id1 in targetData["Protein"].split(","):
        for id2 in modelData["Protein"].split(","):
            rmsd_value = specific_align(1, "tmp/" + targetData["RNA"] + id1 + ".pdb", "tmp/" + modelData["RNA"] + id2 + ".pdb")[1]
            rmsds[id1 + ":" + id2] = rmsd_value

    rmsds = dict(sorted(rmsds.items(), key=lambda x: x[1]))
    matched = []
    mapping = []
    for key in rmsds.keys():
        if(key[0] not in [match[0] for match in mapping] and key[2] not in [match[2] for match in mapping]):
            mapping.append(key)
    chains_mapping_target = ""
    for pair in mapping:
        if (len(chains_mapping_target) > 0): chains_mapping_target += ","
        chains_mapping_target += pair
    return chains_mapping_target


def run_hbplus(name1, name2, target_done):
    commandT = f'bash -i -c "find . -name \'{name1}\' -exec hbplus {{}} \\;"'
    commandM = f'bash -i -c "find . -name \'{name2}\' -exec hbplus {{}} \\;"'
    if(not target_done): subprocess.run(commandT, shell=True)
    subprocess.run(commandM, shell=True)


def trim(s):
    return s.strip()


def getid(don):
    return don[0] + " " + don[1:5] + " " + don[5]


def setDict(t):
    dict = {}
    for i in t:
        dict[i] = 1
    return dict


def remove_leading_0s(s: str) -> str:
    return s.lstrip("0") or "0"


def filter_pairs(file, rna_chains):
    a = []
    dna_str = ["DT", "DA", "DC", "DG"]
    rna_str = ["A", "C", "G", "U"]
    protein_str = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    dna_dict = setDict(dna_str)
    rna_dict = setDict(rna_str)
    protein_dict = setDict(protein_str)
    i=0
    file.seek(0)
    for line in file:
        line.strip()
        if (len(line) != 76 or not line.split()[-1].isdigit()):
            pass
        else:
            donor = line[:9]
            acceptor = line[14:23]

            chain_d = trim(donor[0])
            if (chain_d == "-"):
                chain_d = "A"

            chain_a = trim(acceptor[0])
            if (chain_a == "-"):
                chain_a = "A"

            name_d = trim(donor[6:9])
            name_a = trim(acceptor[6:9])

            if (name_d[0] != "HOH" and name_a != "HOH") and\
            (len(rna_chains) == 0 or (len(rna_chains) > 0 and (rna_chains.find(chain_d) >= 0 or rna_chains.find(chain_a) >= 0))) and\
            (name_d in rna_dict.keys() or name_a in rna_dict.keys()) and \
            not ((name_d in protein_dict.keys() and name_a in protein_dict.keys()) or (name_d in dna_dict.keys() and name_a in dna_dict.keys() > 0) or (name_d in rna_dict.keys() and name_a in rna_dict.keys())):
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


def get_pairs(file, rna_chains, mapping):
    pairs = filter_pairs(file, rna_chains)
    sorted_pairs = sorted(pairs, key=lambda x: (x.split()[0], int(x.split()[1])))
    replaced_pairs = []
    for pair in sorted_pairs:
        replaced_pairs.append(replace_chains(pair, mapping))
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
        if(target[k] in model_dict):
            tp += 1
        else:
            fn += 1
    fp = 0
    for k in range(j):
        if(not (model[k] in target_dict)):
            fp += 1

    ppv = tp / (tp+fp)
    tpr = tp / (tp+fn)
    inf = math.sqrt(ppv*tpr)

    return inf


def convert_to_wsl_path(path):
    if platform.system() == "Windows":
        return subprocess.check_output(["wsl", "wslpath", path.replace("\\", "/")]).decode().strip()
    else:
        return path


def rna_tools_renumerate(filename, output_file, edit_command):
    command = f'rna_pdb_tools.py --edit \'{edit_command}\' {filename} > {output_file}'
    subprocess.run(["bash", "-i", "-c", command])


def find_next_identifier(identifiers):
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    identifiers_set = set(identifiers)
    for letter in alphabet:
        if letter not in identifiers_set:
            return letter


def add_protein_chain_identifier(input_file, output_file, rna_id):
    identifier = "A" if rna_id == "0" else find_next_identifier(rna_id)

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM"):
                id = line[21]
                if id == " ":
                    line = line[:21] + identifier + line[22:]
            outfile.write(line)


def auto_renumber(filename, chains, output):
    command = []
    for chain in chains:
        i = 1
        fragments = find_fragments(filename, chain)
        chain_mapping = []
        for fragment in fragments:
            fragment_mapping = chain + ":" + str(fragment[0]) + "-" + str(fragment[1]) + ">" + chain + ":" + str(i) + "-" + str(fragment[1] - fragment[0] + i)
            chain_mapping.append(fragment_mapping)
            i = fragment[1] - fragment[0] + i + 1
        command.append(",".join(chain_mapping))
    edit_command = ",".join(command)
    rna_tools_renumerate(filename, output, edit_command)


def find_fragments(filename, chain):
    numbers = []
    with open(filename) as file:
        for line in file:
            line_tmp = " ".join(line.split())
            if (line_tmp.split()[0] == "ATOM" and line_tmp.split()[4] == chain):
                if(len(numbers) == 0):
                    i = line_tmp.split()[5]
                    numbers.append(int(line_tmp.split()[5]))
                if(line_tmp.split()[5] != i): numbers.append(int(line_tmp.split()[5]))
                i = line_tmp.split()[5]
        if(len(numbers) > 0): fragments = [[numbers[0]]]
        for i in range(1, len(numbers)):
            if(numbers[i] > numbers[i-1] + 1):
                fragments[len(fragments)-1].append(numbers[i-1])
                fragments.append([numbers[i]])
        fragments[len(fragments)-1].append(numbers[len(numbers)-1])
        return fragments


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


def compare (name1, names2, rmsd, clear_target, raw_inf):
    files_model = []
    files_target = []
    infs = []
    target_done = False

    for name in names2:
        name2 = name
        if(not target_done): name1_copied = copy_file_to_script_directory(name1)
        name2_copied = copy_file_to_script_directory(name2)
        if(name1_copied):
            name1 = os.path.basename(name1)
        if(name2_copied):
            name2 = os.path.basename(name2)

        target = open(name1, "r")
        model = open(name2, "r")
        if(name1_copied): files_target.append(target)
        if(name2_copied): files_model.append(model)

        if(clear_target and not target_done):
            clear(name1)
            name1 = name1[0:len(name1)-4] + "m" + name1[len(name1)-4:]
            target = open(name1, "r")
            files_target.append(target)

        clear(name2)
        name2 = name2[0:len(name2)-4] + "m" + name2[len(name2)-4:]
        model = open(name2, "r")
        files_model.append(model)

        print("Target:")
        targetData = analyze(name1)
        print("")
        print(targetData)
        print("")
        print("Model:")
        modelData = analyze(name2)
        print("")
        print(modelData)
        print("")

        if(modelData["single_protein"] == True):
            output_file = name2[0:len(name2)-4] + "_id" + name2[len(name2)-4:]
            add_protein_chain_identifier(name2, output_file, modelData["RNA"])
            name2 = name2[0:len(name2)-4] + "_id" + name2[len(name2)-4:]
            model = open(name2, "r")
            files_model.append(model)
            modelData = analyze(name2)
            print("")
            print(modelData)
            print("")

        if (not modelData["RNA"].isalpha()):
            id = find_next_identifier(modelData["Protein"])
            rep = modelData["RNA"] + ">" + id
            output = name2[0:len(name2)-4] + "_t" + name2[len(name2)-4:]
            command = f'rna_pdb_tools.py --rename-chain \'{rep}\' {name2} > {output}'
            subprocess.run(["bash", "-i", "-c", command])
            name2 = output
            model = open(name2, "r")
            files_model.append(model)
            modelData = analyze(name2)
            print("Model (renamed):")
            print("")
            print(modelData)
            print("")

        target_renum = []
        model_renum = []
        target_to_renum = []
        model_to_renum = []
        renumber_model = False
        renumber_target = False
        custom_renum_target = False
        custom_renum_model = False

        if ((targetData[targetData["RNA"]][2] - targetData[targetData["RNA"]][1] + 1 != targetData[targetData["RNA"]][3]) or
                targetData[targetData["RNA"]][1] != 1 and not target_done):
            renumber_target = True
            print("Mismatched RNA numbering (target). Auto-renumber? (Y/N)")
            ans = input()
            if (ans == "Y"):
                target_to_renum.append(targetData[targetData["RNA"]][0])
            else:
                print("Enter custom alignement:")
                custom_renum_target = True
                ans = input()
                renum = ans[0]
                model_renum.append(renum)

        if ((modelData[modelData["RNA"]][2] - modelData[modelData["RNA"]][1] + 1 != modelData[modelData["RNA"]][3] or
             modelData[modelData["RNA"]][1] != 1)):
            renumber_model = True
            print("Mismatched RNA numbering (model). Auto-renumber? (Y/N)")
            ans = input()
            if (ans == "Y"):
                model_to_renum.append(modelData[modelData["RNA"]][0])
            else:
                print("Enter custom alignement:")
                custom_renum_model = True
                ans = input()
                renum = ans[0]
                model_renum.append(renum)

        for id in modelData["Protein"].split(","):
            if (modelData[id][2] - modelData[id][1] + 1 != modelData[id][3] or modelData[id][1] != 1):
                renumber_model = True
                print("Mismatched Protein Chain " + id + " numbering (model). Auto-renumber? (Y/N)")
                ans = input()
                if (ans == "Y"):
                    model_to_renum.append(modelData[id][0])
                else:
                    print("Enter custom alignement:")
                    custom_renum_model = True
                    ans = input()
                    model_renum.append(ans)

        for id in targetData["Protein"].split(","):
            if (targetData[id][2] - targetData[id][1] + 1 != targetData[id][3] or targetData[id][1] != 1 and not target_done):
                renumber_target = True
                print("Mismatched Protein Chain " + id + " numbering (target) Auto-renumber? (Y/N)")
                ans = input()
                if (ans == "Y"):
                    target_to_renum.append(targetData[id][0])
                else:
                    print("Enter custom alignement:")
                    custom_renum_target = True
                    ans = input()
                    target_renum.append(ans)

        if (renumber_target):
            auto_renumber(name1, target_to_renum, name1[0:len(name1)-4] + "_ren" + name1[len(name1)-4:])
            name1 = name1[0:len(name1)-4] + "_ren" + name1[len(name1)-4:]
            target = open(name1, "r")
            files_target.append(target)

        if (renumber_model):
            auto_renumber(name2, model_to_renum, name2[0:len(name2)-4] + "_ren" + name2[len(name2)-4:])
            name2 = name2[0:len(name2)-4] + "_ren" + name2[len(name2)-4:]
            model = open(name2, "r")
            files_model.append(model)

        if(custom_renum_target):
            rna_tools_renumerate(name2, name1[0:len(name1)-4] + "_cus" + name1[len(name1)-4:], ",".join(target_renum))
            name1 = name1[0:len(name1)-4] + "_cus" + name1[len(name1)-4:]
            target = open(name1, "r")
            files_target.append(target)

        if (custom_renum_model):
            rna_tools_renumerate(name2, name2[0:len(name2)-4] + "_cus" + name2[len(name2)-4:], ",".join(model_renum))
            name2 = name2[0:len(name2)-4] + "_cus" + name2[len(name2)-4:]
            model = open(name2, "r")
            files_model.append(model)

        rep = modelData["RNA"] + ">" + "0"
        output = name2[0:len(name2)-4] + "_t" + name2[len(name2)-4:]
        command = f'rna_pdb_tools.py --rename-chain \'{rep}\' {name2} > {output}'
        subprocess.run(["bash", "-i", "-c", command])
        name2 = output
        model = open(name2, "r")
        files_model.append(model)
        modelData = analyze(name2)

        run_hbplus(name1, name2, target_done)
        name3 = name1[0:len(name1) - 3] + "hb2"
        name4 = name2[0:len(name2) - 3] + "hb2"
        target_HB2 = open(name3, "r")
        if(not target_done): target_done = True
        model_HB2 = open(name4, "r")
        files_target.append(target_HB2)
        files_model.append(model_HB2)

        os.makedirs("tmp", exist_ok=True)

        chains_mapping_model = modelData["RNA"] + ":" + targetData["RNA"]

        if(not rmsd):
            target_mappings = create_combinations(targetData["Protein"].split(","), modelData["Protein"].split(","))
            max_inf = 0
            max_mapping = ""
            for mapping in target_mappings:
                target_pairs = get_pairs(target_HB2, "A", mapping.split(","))
                model_pairs = get_pairs(model_HB2, "0", chains_mapping_model.split(","))
                if (len(model_pairs) == 0):
                    pass
                elif (len(target_pairs) == 0):
                    pass
                else:
                    inf = get_inf(target_pairs, model_pairs)
                    if(inf > max_inf):
                        max_inf = inf
                        max_mapping = mapping

            inf = max_inf
        else:
            chains_mapping_target = get_mapping(targetData, modelData, target, model)
            target_pairs = get_pairs(target_HB2, "A", chains_mapping_target.split(","))
            model_pairs = get_pairs(model_HB2, "0", chains_mapping_model.split(","))
            inf = get_inf(target_pairs, model_pairs)
        print("")
        if(raw_inf):
            print(format(inf, ".3f"))
            infs.append([name, format(inf, ".3f")])
        else:
            inf = inf * modelData["residues_no"] / targetData["residues_no"]
            print(format(inf, ".3f"))
            infs.append([name, format(inf, ".3f")])

        tmp_files = glob.glob(os.path.join("tmp", "*"))
        for tmp_file in tmp_files:
            pass
            try:
                os.remove(tmp_file)
            except Exception as e:
                print(f'Error deleting {tmp_file}: {e}')

        if not os.listdir("tmp"):
            os.rmdir("tmp")

        for file in files_model:
            file.close()
            file_name = file.name
            try:
                os.remove(file_name)
                files_model.remove(file)
            except Exception as e:
                print(f"Error deleting {file_name}: {e}")

        try:
            os.remove("hbdebug.dat")
        except Exception as e:
            print(f"Error deleting {'hbdebug.dat'}: {e}")

    for file in files_target:
        file.close()
        file_name = file.name
        try:
            os.remove(file_name)
        except Exception as e:
            print(f"Error deleting {file_name}: {e}")

    return infs


name1 = input()
name2 = input()

adj = False
print("In case of mismatched residues number, adjust INF (INF = INF * res_model/res_target)?: (Y/N)")
ans = input()
if(ans == "Y"): adj = True

if(os.path.isdir(name2)):
    files_to_compare = [os.path.join(name2, f) for f in os.listdir(name2) if f.endswith('.pdb')]
    print(files_to_compare)
    infs = compare(name1, files_to_compare, False, False, adj)
else:
    infs = compare(name1, [name2], False, False, False)
print("")
print(infs)