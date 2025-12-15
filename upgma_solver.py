import itertools
from io import StringIO
from Bio import Phylo
import matplotlib.pyplot as plt
import re
from prettytable import PrettyTable
import numpy as np
import json
import sys
import os

def format_number(num):
    if num == 0:
        return "0"
    s = f"{num:.2f}".rstrip('0').rstrip('.')
    return s

def print_distance_matrix(matrix, step_num):
    print(f"\n{'='*60}")
    print(f"ÉTAPE {step_num}")
    print('='*60)
    print("\nMatrice distance\n")
    labels = sorted(matrix.keys())
    
    table = PrettyTable()
    
    field_names = ["D"] + [str(label) for label in labels]
    table.field_names = field_names
    
    for i, label in enumerate(labels):
        row = [str(label)]
        for j, other in enumerate(labels):
            if i > j:
                row.append("")
            elif label == other:
                row.append("0")
            elif other in matrix[label]:
                row.append(format_number(matrix[label][other]))
            else:
                row.append("N/A")
        table.add_row(row)
    
    print(table)

def print_clusters(clusters):
    print(f"Clusters actuels ({len(clusters)}):")
    table = PrettyTable()
    table.field_names = ["#", "Cluster"]
    table.align["Cluster"] = "l"
    
    for i, cluster in enumerate(sorted(clusters), 1):
        table.add_row([i, cluster])
    
    print(table)
    print()

def print_clusters_with_trees(cluster_list):
    table = PrettyTable()
    table.field_names = ["Nom", "Contenu"]
    table.align["Nom"] = "l"
    table.align["Contenu"] = "l"
    
    for cluster_name, tree_str in cluster_list:
        table.add_row([cluster_name, tree_str + ";"])
    
    print(table)

def newick_without_distances(newick_str):
    return re.sub(r':[0-9.]+', '', newick_str)

def calculate_tree_height(clade):
    if clade.is_terminal():
        return clade.branch_length if clade.branch_length is not None else 0.0
    else:
        branch_len = clade.branch_length if clade.branch_length is not None else 0.0
        max_child_height = max([calculate_tree_height(child) for child in clade.clades], default=0.0)
        return branch_len + max_child_height

def print_tree_with_lengths(tree_obj, cluster_name, newick_str):
    print(f"\n{cluster_name}:")
    
    root = tree_obj.root
    tree_height = calculate_tree_height(root)
    
    if tree_height == 0.0:
        lengths = re.findall(r':([0-9.]+)', newick_str)
        if lengths:
            length_values = [float(l) for l in lengths]
            if len(length_values) == 2:
                tree_height = max(length_values)
            else:
                main_branch_lengths = re.findall(r'\):([0-9.]+)', newick_str)
                if main_branch_lengths:
                    tree_height = sum([float(l) for l in main_branch_lengths])
                else:
                    tree_height = max(length_values)
    
    if tree_height > 0:
        length_str = format_number(tree_height)
        
        from io import StringIO
        import sys
        old_stdout = sys.stdout
        sys.stdout = tree_output = StringIO()
        Phylo.draw_ascii(tree_obj)
        sys.stdout = old_stdout
        tree_lines = tree_output.getvalue().split('\n')
        
        max_width = 0
        for line in tree_lines:
            if line.strip():
                line_width = len(line.rstrip())
                max_width = max(max_width, line_width)
        
        if max_width == 0:
            max_width = 60
        
        total_dashes = max_width - len(length_str) - 15
        if total_dashes < 0:
            total_dashes = 0
        
        dashes_left = "-" * (total_dashes // 2)
        dashes_right = "-" * (total_dashes - len(dashes_left))
        
        print(f"  <{dashes_left}    {length_str}    {dashes_right}>")
    
    Phylo.draw_ascii(tree_obj)

def print_tree(newick_str, step_num):
    print(f"Arbre actuel (étape {step_num}):")
    print(f"  Newick avec distances: {newick_str}")
    newick_no_dist = newick_without_distances(newick_str)
    print(f"  Newick sans distances: {newick_no_dist}")
    
    try:
        tree_obj = Phylo.read(StringIO(newick_str), "newick")
        print("\n  Vue ASCII:")
        Phylo.draw_ascii(tree_obj)
    except:
        pass
    print()

def get_all_upgma_trees(distance_matrix):
    labels = list(distance_matrix.keys())
    initial_sizes = {label: 1 for label in labels}
    initial_heights = {label: 0.0 for label in labels}
    found_trees = set()
    step_counter = [0]
    cluster_counter = [0]

    def solve(curr_matrix, curr_sizes, curr_heights, tree_map=None, display_names=None):
        step_counter[0] += 1
        step = step_counter[0]
        
        if tree_map is None:
            tree_map = {label: label for label in curr_matrix.keys()}
        if display_names is None:
            display_names = {label: label for label in curr_matrix.keys()}
        reverse_display_names = {v: k for k, v in display_names.items()}
        
        display_matrix = {}
        for internal_name in curr_matrix.keys():
            display_name = display_names.get(internal_name, internal_name)
            display_matrix[display_name] = {}
            for other_internal in curr_matrix[internal_name]:
                other_display = display_names.get(other_internal, other_internal)
                display_matrix[display_name][other_display] = curr_matrix[internal_name][other_internal]
        
        print_distance_matrix(display_matrix, step)
        
        cluster_list = [(display_names.get(name, name), tree_map.get(name, name)) 
                        for name in sorted(curr_matrix.keys())]
        print("\nGroupements (clusters)\n")
        print_clusters_with_trees(cluster_list)
        
        if len(curr_matrix) > 1:
            print("Arbres courants:")
            for cluster_name, tree_str in cluster_list:
                if "(" not in tree_str and ")" not in tree_str:
                    print(f"{cluster_name}: {tree_str}")
                else:
                    try:
                        tree_obj = Phylo.read(StringIO(tree_str + ";"), "newick")
                        print_tree_with_lengths(tree_obj, cluster_name, tree_str)
                    except:
                        pass
            print()
        elif len(curr_matrix) == 1:
            print("Arbres courants:")
            cluster_name, tree_str = cluster_list[0]
            try:
                tree_obj = Phylo.read(StringIO(tree_str + ";"), "newick")
                print_tree_with_lengths(tree_obj, cluster_name, tree_str)
            except:
                pass
            print()
        
        if len(curr_matrix) == 1:
            root = list(curr_matrix.keys())[0]
            final_tree = tree_map[root] + ";"
            found_trees.add(final_tree)
            
            root_display_name = display_names.get(root, root)
            print(f"\nARBRE FINAL (Cluster: {root_display_name}):")
            print(f"  Newick avec distances: {final_tree}")
            print(f"  Newick sans distances: {newick_without_distances(final_tree)}")
            
            return True

        labels = list(curr_matrix.keys())
        min_dist = float('inf')
        
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                d = curr_matrix[labels[i]][labels[j]]
                if d < min_dist: min_dist = d
        
        ties = []
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if curr_matrix[labels[i]][labels[j]] == min_dist:
                    ties.append((labels[i], labels[j]))
        
        ties_display = [(display_names.get(u, u), display_names.get(v, v)) 
                       for (u, v) in ties]
        print(f"Distance minimale: {min_dist:.2f} | Fusions possibles: {ties_display}")

        for (u, v) in ties:
            new_height = min_dist / 2.0
            len_u = new_height - curr_heights[u]
            len_v = new_height - curr_heights[v]
            
            tree_u = tree_map[u]
            tree_v = tree_map[v]
            
            child_strs = sorted([f"{tree_u}:{len_u:g}", f"{tree_v}:{len_v:g}"])
            new_name = f"({child_strs[0]},{child_strs[1]})"
            new_size = curr_sizes[u] + curr_sizes[v]
            
            u_display = display_names.get(u, u)
            v_display = display_names.get(v, v)
            
            if u_display.startswith('C') and v_display.startswith('C'):
                cluster_counter[0] += 1
                new_display_name = f"C{cluster_counter[0]}"
            elif u_display.startswith('C'):
                new_display_name = u_display
            elif v_display.startswith('C'):
                new_display_name = v_display
            else:
                cluster_counter[0] += 1
                new_display_name = f"C{cluster_counter[0]}"
            
            print(f"Fusion choisie: {u_display} + {v_display} -> {new_display_name} (Newick: {newick_without_distances(new_name)})")
            
            new_mat = {}
            new_sizes = curr_sizes.copy()
            new_heights = curr_heights.copy()
            new_tree_map = tree_map.copy()
            new_display_names = display_names.copy()
            
            new_sizes[new_name] = new_size
            new_heights[new_name] = new_height
            new_tree_map[new_name] = new_name
            new_display_names[new_name] = new_display_name
            
            del new_sizes[u], new_sizes[v]
            del new_heights[u], new_heights[v]
            del new_tree_map[u], new_tree_map[v]
            del new_display_names[u], new_display_names[v]
            
            others = [k for k in labels if k != u and k != v]
            new_mat[new_name] = {}
            new_mat[new_name][new_name] = 0.0
            
            if others:
                print(f"  Calcul des distances pour {new_display_name}:")
            
            for other in others:
                if other not in new_mat: 
                    new_mat[other] = {}
                    new_mat[other][other] = 0.0
                dist_u = curr_matrix[u][other]
                dist_v = curr_matrix[v][other]
                dist_new = (dist_u * curr_sizes[u] + dist_v * curr_sizes[v]) / new_size
                new_mat[new_name][other] = dist_new
                new_mat[other][new_name] = dist_new
                
                other_display = display_names.get(other, other)
                u_display = display_names.get(u, u)
                v_display = display_names.get(v, v)
                print(f"    d({new_display_name}, {other_display}) = (d({u_display}, {other_display}) × |{u_display}| + d({v_display}, {other_display}) × |{v_display}|) / |{new_display_name}| = ({format_number(dist_u)} × {curr_sizes[u]} + {format_number(dist_v)} × {curr_sizes[v]}) / {new_size} = {format_number(dist_new)}")
                
                for other2 in others:
                    if other2 in curr_matrix[other]:
                        new_mat[other][other2] = curr_matrix[other][other2]
            
            if others:
                print()

            if solve(new_mat, new_sizes, new_heights, new_tree_map, new_display_names):
                return True
        
        return False

    solve(distance_matrix, initial_sizes, initial_heights)
    return list(found_trees)

def load_distance_matrix(filename):
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        
        labels = sorted(data.keys())
        
        for i, label1 in enumerate(labels):
            if label1 not in data[label1] or data[label1][label1] != 0:
                if label1 not in data:
                    data[label1] = {}
                data[label1][label1] = 0
            
            for j, label2 in enumerate(labels):
                if i != j:
                    if label2 not in data[label1]:
                        if label1 in data.get(label2, {}):
                            if label1 not in data:
                                data[label1] = {}
                            data[label1][label2] = data[label2][label1]
                        else:
                            print(f"Erreur: Distance manquante entre '{label1}' et '{label2}'")
                            sys.exit(1)
                    if label1 not in data.get(label2, {}):
                        if label2 not in data:
                            data[label2] = {}
                        data[label2][label1] = data[label1][label2]
        
        return data
    except FileNotFoundError:
        print(f"Erreur: Le fichier '{filename}' n'existe pas.")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Erreur: Le fichier '{filename}' n'est pas un JSON valide.")
        sys.exit(1)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier: {e}")
        sys.exit(1)

if len(sys.argv) > 1:
    matrix_file = sys.argv[1]
else:
    matrix_file = "distance_matrix.json"

data = load_distance_matrix(matrix_file)

print("DÉBUT DU CALCUL UPGMA\n")
trees = get_all_upgma_trees(data)

print(f"\nRÉSUMÉ: {len(trees)} arbre(s) unique(s) trouvé(s)\n")

for i, newick_tree in enumerate(trees):
    print(f"ARBRE FINAL #{i+1}:")
    print(f"  Newick avec distances: {newick_tree}")
    newick_no_dist = newick_without_distances(newick_tree)
    print(f"  Newick sans distances: {newick_no_dist}")
    
    tree_obj = Phylo.read(StringIO(newick_tree), "newick")
    
    print("\nVue ASCII:")
    newick_for_lengths = newick_tree.rstrip(';')
    print_tree_with_lengths(tree_obj, f"Arbre #{i+1}", newick_for_lengths)
    
    base_name = os.path.splitext(os.path.basename(matrix_file))[0]
    if len(trees) > 1:
        image_filename = f"{base_name}_tree_{i+1}.png"
    else:
        image_filename = f"{base_name}_tree.png"
    
    print(f"\nSauvegarde: {image_filename}")
    fig, ax = plt.subplots(figsize=(6, 4))
    Phylo.draw(tree_obj, axes=ax, do_show=False)
    plt.title(f"UPGMA Tree #{i+1}")
    plt.savefig(image_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
