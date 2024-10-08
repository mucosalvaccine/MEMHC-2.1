# root.mainloop()
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText
import pandas as pd
from mhcflurry import Class1AffinityPredictor
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import threading
import random
import math

# Global variables for cumulative data
cumulative_x_data = []
cumulative_y_data = []
data_store = {}


# Utility function to log messages in the UI
def log_message(message):
    log_text.insert(tk.END, message + "\n")
    log_text.see(tk.END)


# Function to generate peptides from the protein sequence
def generate_peptides(protein_sequence, min_length, max_length):
    peptides = []
    start_positions = []
    for i in range(len(protein_sequence)):
        for j in range(i + min_length, min(i + max_length + 1, len(protein_sequence) + 1)):
            peptides.append(protein_sequence[i:j])
            start_positions.append(i)
    return peptides, start_positions


# Function to update the cumulative coverage chart
def update_chart(x_data, y_data):
    global cumulative_x_data, cumulative_y_data

    ax.clear()
    cumulative_x_data.extend(x_data)
    cumulative_y_data.extend(y_data)

    ax.plot(cumulative_x_data, cumulative_y_data, color='blue', marker='o', linestyle='-', linewidth=2)
    ax.fill_between(cumulative_x_data, cumulative_y_data, color='blue', alpha=0.7)

    ax.set_xlabel('Percentage of Protein Area Covered')
    ax.set_ylabel('Percentage of HLA Coverage')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_title("Protein vs HLA Coverage")

    canvas.draw()


# Function to display the coverage of peptides on the protein sequence
def display_protein_matrix(protein_sequence, covered_positions):
    coverage_array = [1 if i in covered_positions else 0 for i in range(len(protein_sequence))]

    matrix_ax.clear()
    matrix_ax.imshow([coverage_array], cmap='Blues', aspect='auto', interpolation='nearest')

    matrix_ax.set_xticks(range(len(protein_sequence)))
    matrix_ax.set_xticklabels([aa for aa in protein_sequence])
    matrix_ax.set_yticks([])

    matrix_ax.set_title("Peptide Coverage in Protein Sequence")
    matrix_canvas.draw()


# Function to update the results table
def update_result_table(pep, hits, absolute_coverage, cumulative_coverage, protein_coverage):
    tree.insert("", "end", values=(
        pep, ", ".join(hits), f"{absolute_coverage:.2f}", f"{cumulative_coverage:.2f}", f"{protein_coverage:.2f}"))


# File dialog to browse for HLA file
def browse_file():
    filename = filedialog.askopenfilename()
    hla_file_entry.delete(0, tk.END)
    hla_file_entry.insert(0, filename)


# Save results to a file (CSV or Excel)
def save_file():
    rows = tree.get_children()
    data = []
    for row in rows:
        data.append(tree.item(row)["values"])

    if not data:
        messagebox.showerror("Error", "No data available to save.")
        return

    df_result = pd.DataFrame(data, columns=['Peptide', 'HLA Hits', 'Absolute Coverage (%)', 'Cumulative Coverage (%)',
                                            'Protein Coverage (%)'])
    filetypes = [('CSV Files', '*.csv'), ('Excel Files', '*.xlsx')]
    filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=filetypes)
    if not filepath:
        return

    if filepath.endswith('.csv'):
        df_result.to_csv(filepath, index=False)
    else:
        df_result.to_excel(filepath, index=False)

    messagebox.showinfo("Save Successful", f"File saved as {filepath}")


# Function to run Simulated Annealing for peptide selection
def simulated_annealing(peptides, hla_hits, hla_alleles, max_iterations=1000, start_temp=1000, end_temp=1, alpha=0.99):
    def calculate_coverage(selected_peptides):
        covered_hla = set()
        for pep in selected_peptides:
            covered_hla.update(hla_hits.get(pep, []))
        return covered_hla

    def swap_peptides(selected_peptides, peptides_pool):
        new_selection = selected_peptides.copy()
        if len(selected_peptides) > 0:
            peptide_to_remove = random.choice(list(selected_peptides))
            new_selection.remove(peptide_to_remove)
        peptide_to_add = random.choice(peptides_pool)
        new_selection.add(peptide_to_add)
        return new_selection

    current_peptides = set(random.sample(peptides, len(hla_alleles)))
    current_coverage = calculate_coverage(current_peptides)
    best_peptides = current_peptides
    best_coverage = current_coverage

    current_temp = start_temp
    iteration = 0

    while current_temp > end_temp and iteration < max_iterations:
        iteration += 1
        new_peptides = swap_peptides(current_peptides, peptides)
        new_coverage = calculate_coverage(new_peptides)

        delta_coverage = len(new_coverage) - len(current_coverage)
        if delta_coverage > 0:
            acceptance_probability = 1.0
        else:
            acceptance_probability = math.exp(delta_coverage / current_temp)

        if delta_coverage > 0 or random.random() < acceptance_probability:
            current_peptides = new_peptides
            current_coverage = new_coverage

            if len(current_coverage) > len(best_coverage):
                best_peptides = current_peptides
                best_coverage = current_coverage

        current_temp *= alpha

    return best_peptides, best_coverage


# Function to run the main Simulated Annealing optimization
def run_simulated_annealing():
    global cumulative_x_data, cumulative_y_data
    cumulative_x_data.clear()
    cumulative_y_data.clear()

    log_text.delete(1.0, tk.END)
    log_message("Script started with Simulated Annealing.")

    for row in tree.get_children():
        tree.delete(row)

    protein_sequence = protein_entry.get().upper()
    min_length = int(min_length_entry.get())
    max_length = int(max_length_entry.get())
    hla_file = hla_file_entry.get()
    kd_threshold = int(kd_threshold_entry.get())
    hla_coverage_limit = int(hla_coverage_limit_entry.get())

    if not protein_sequence or not hla_file:
        messagebox.showerror("Input Error", "Please provide all required inputs.")
        return

    accepted_aas = set("ACDEFGHIKLMNPQRSTVWY")
    if set(protein_sequence) - accepted_aas:
        messagebox.showerror("Input Error", "Protein sequence contains invalid amino acids.")
        return

    log_message("Generating peptides...")
    peptides, start_positions = generate_peptides(protein_sequence, min_length, max_length)

    log_message("Predicting binding affinities...")
    predictor = Class1AffinityPredictor.load()
    HLA_type_I = pd.read_csv(hla_file)
    HLA_type_I_list = list(HLA_type_I.columns)
    total_hla_alleles = len(HLA_type_I_list)

    binding_predictions = {}
    for mhc_allele in HLA_type_I_list:
        try:
            predicted_affinities = predictor.predict(peptides, [mhc_allele] * len(peptides))
            for pep, affinity in zip(peptides, predicted_affinities):
                if pep not in binding_predictions:
                    binding_predictions[pep] = {}
                binding_predictions[pep][mhc_allele] = affinity
        except ValueError:
            pass

    hla_hits = {}
    for pep in peptides:
        hla_hits[pep] = [mhc_allele for mhc_allele in HLA_type_I_list if
                         binding_predictions[pep].get(mhc_allele, float('inf')) <= kd_threshold]

    best_peptides, best_coverage = simulated_annealing(peptides, hla_hits, HLA_type_I_list)

    covered_positions = set()
    total_covered_alleles = best_coverage
    cumulative_protein_coverage = 0

    for pep in best_peptides:
        start_pos = start_positions[peptides.index(pep)]
        new_positions = set(range(start_pos, start_pos + len(pep)))

        unique_new_positions = new_positions - covered_positions
        covered_positions.update(unique_new_positions)

        protein_coverage = (len(unique_new_positions) / len(protein_sequence)) * 100
        cumulative_protein_coverage = (len(covered_positions) / len(protein_sequence)) * 100

        cumulative_coverage = len(total_covered_alleles) / total_hla_alleles * 100

        update_result_table(pep, hla_hits[pep], len(hla_hits[pep]) / total_hla_alleles * 100, cumulative_coverage, cumulative_protein_coverage)
        update_chart([cumulative_protein_coverage], [cumulative_coverage])
        display_protein_matrix(protein_sequence, covered_positions)

    log_message("Simulated Annealing finished successfully.")
    messagebox.showinfo("Success", "The Simulated Annealing optimization has completed.")


# GUI setup with Tkinter
root = tk.Tk()
root.title("Peptide Generator and HLA Affinity Predictor")
root.geometry("1024x800")

# Input fields
tk.Label(root, text="Protein Sequence:").grid(row=0, column=0, sticky="e")
protein_entry = tk.Entry(root, width=50)
protein_entry.grid(row=0, column=1, padx=10, pady=10)

tk.Label(root, text="Min Peptide Length:").grid(row=1, column=0, sticky="e")
min_length_entry = tk.Entry(root)
min_length_entry.grid(row=1, column=1, padx=10, pady=10)
min_length_entry.insert(0, "8")

tk.Label(root, text="Max Peptide Length:").grid(row=2, column=0, sticky="e")
max_length_entry = tk.Entry(root)
max_length_entry.grid(row=2, column=1, padx=10, pady=10)
max_length_entry.insert(0, "9")

tk.Label(root, text="HLA File:").grid(row=3, column=0, sticky="e")
hla_file_entry = tk.Entry(root, width=50)
hla_file_entry.grid(row=3, column=1, padx=10, pady=10)
tk.Button(root, text="Browse", command=browse_file).grid(row=3, column=2, padx=10, pady=10)

tk.Label(root, text="Kd Threshold (nM):").grid(row=4, column=0, sticky="e")
kd_threshold_entry = tk.Entry(root)
kd_threshold_entry.grid(row=4, column=1, padx=10, pady=10)
kd_threshold_entry.insert(0, "500")

tk.Label(root, text="HLA Coverage Limit:").grid(row=5, column=0, sticky="e")
hla_coverage_limit_entry = tk.Entry(root)
hla_coverage_limit_entry.grid(row=5, column=1, padx=10, pady=10)
hla_coverage_limit_entry.insert(0, "1")

# Buttons
run_button = tk.Button(root, text="Run Simulated Annealing", command=lambda: threading.Thread(target=run_simulated_annealing).start())
run_button.grid(row=7, column=0, sticky="w", padx=10, pady=10)

abort_button = tk.Button(root, text="Abort", command=root.quit)
abort_button.grid(row=7, column=1, sticky="e", padx=10, pady=10)

# Treeview to display the results
tree = ttk.Treeview(root, height=10)
tree["columns"] = ("Peptide", "HLA Hits", "Absolute Coverage (%)", "Cumulative Coverage (%)", "Protein Coverage (%)")
tree.column("#0", width=0, stretch=tk.NO)
tree.column("Peptide", anchor=tk.W, width=150)
tree.column("HLA Hits", anchor=tk.W, width=300)
tree.column("Absolute Coverage (%)", anchor=tk.CENTER, width=150)
tree.column("Cumulative Coverage (%)", anchor=tk.CENTER, width=150)
tree.column("Protein Coverage (%)", anchor=tk.CENTER, width=150)

tree.heading("Peptide", text="Peptide", anchor=tk.W)
tree.heading("HLA Hits", text="HLA Hits", anchor=tk.W)
tree.heading("Absolute Coverage (%)", text="Absolute Coverage (%)", anchor=tk.CENTER)
tree.heading("Cumulative Coverage (%)", text="Cumulative Coverage (%)", anchor=tk.CENTER)
tree.heading("Protein Coverage (%)", text="Protein Coverage (%)", anchor=tk.CENTER)
tree.grid(row=9, column=0, columnspan=2, padx=10, pady=10, sticky="nsew")

# Save button
save_button = tk.Button(root, text="Download", command=save_file)
save_button.grid(row=10, column=1, padx=10, pady=10)

# Log panel
log_frame = tk.Frame(root)
log_text = ScrolledText(log_frame, height=10, state='normal')
log_text.pack(fill=tk.BOTH, expand=True)
log_frame.grid(row=11, column=0, columnspan=2, sticky="nsew")

# Live updating chart for HLA coverage vs Protein Area
fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=2, column=3, rowspan=6, padx=10, pady=10, sticky="nsew")

# Add another figure for the protein position matrix
matrix_fig, matrix_ax = plt.subplots(figsize=(10, 1), dpi=100)
matrix_canvas = FigureCanvasTkAgg(matrix_fig, master=root)
matrix_canvas.get_tk_widget().grid(row=8, column=3, rowspan=2, padx=10, pady=10, sticky="nsew")

# Start Tkinter event loop
root.mainloop()
