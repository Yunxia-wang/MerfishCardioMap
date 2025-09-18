
# 📘 MerfishCardioMap

### Overview  
**MerfishCardioMap** is a research toolkit for integrating **scRNA-seq** data with **MERFISH spatial transcriptomics** data to explore spatial gene expression patterns in cardiac development. The repository provides data analysis scripts and utilities for processing, visualization, and downstream biological interpretation.  

Key features include:  
- Integration with scRNA-seq reference datasets  
- Preprocessing of MERFISH spatial data    
- Visualization tools for single-cell and spatial domains  

---

### 🛠 Installation Guide  

#### Clone the repository  
```bash
git clone https://github.com/Yunxia-wang/MerfishCardioMap.git
cd MerfishCardioMap
```

#### Create the conda environment  
```bash
conda env create -f environment_merfish.yml
conda activate merfish
```

Dependencies include:  
- Python (≥3.9)  
- numpy, pandas, scanpy, anndata  
- seaborn, matplotlib  
- jupyter, notebook  

---

### ⚙️ Setting up the development environment  
1. Install [VS Code](https://code.visualstudio.com/) or use JupyterLab.  
2. Open the repo folder in VS Code.  
3. Select the `merfish` conda environment as your Python interpreter.  
4. Run Jupyter notebooks such as:  
   - `Main_get_single_cell_region.ipynb` → main entry point for analysis  
   - Additional notebooks inside `/scRNA_seq/` and `/MERFISH/` folders  

---

### 📂 Repository Structure  

```
MerfishCardioMap/
│
├── MERFISH/                     # MERFISH data and analysis scripts
├── scRNA_seq/                   # scRNA-seq data and integration
├── images/                      # Figures and plots
├── environment_merfish.yml      # Conda environment file
├── Main_get_single_cell_region.ipynb  # Main notebook
├── README.md                    # Project documentation
├── LICENSE                      # License (MIT, GPL, etc.)
│
├── cell_by_gene_3D_embryo4.csv        # Example MERFISH cell-by-gene matrix
├── cell_metadata_3D_embryo4.csv       # Cell-level metadata
├── cellpose_micron_space.parquet      # Processed cell spatial positions
├── multipolygon_plot.pdf              # Example visualization
```

---

### 📊 Data position  

- **`cell_by_gene_3D_embryo4.csv`** → Raw MERFISH counts (cells × genes)  
- **`cell_metadata_3D_embryo4.csv`** → Cell annotations and metadata  
- **`cellpose_micron_space.parquet`** → Spatial coordinates of cells  
- **`images/`** → Example output figures  

> ⚠️ Data included here are processed subsets for reproducibility. Single-cell RNA sequencing data have been deposited in the GEO under accession number GSE294086. MERFISH data have been deposited in the GEO under accession number GSE294469. 

---

### 🚀 Usage Example  

Run the main notebook to perform spatial integration and clustering:  
```bash
jupyter notebook Main_get_single_cell_region.ipynb
```

Inside the notebook, you can:  
- Load MERFISH and scRNA-seq data  
- Visualize UMAPs and spatial plots  

---

### 📜 License  
This project is licensed under the [Apache 2.0 LICENSE](LICENSE).  

---

