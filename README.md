# 🧬 AlphaFold Structure Downloader & Viewer

This is a lightweight script to fetch **AlphaFold-predicted protein structures** using UniProt IDs and visualize them in 3D — right from your terminal.

No deep learning, no setup hassle — just instant access to cutting-edge structural biology.


## ⚙️ How to Use

### 1. Clone the Repository
```bash
git clone https://github.com/dhanisthsarawagi/Alpha-fold-protien.git
cd Alpha-fold-protien
2. Install Dependencies
pip install py3Dmol
3. Run the Script
python alpha.py

💡change the uniprot_id directly in the script if you want to:

# Examples:
P38398  → BRCA1 (Breast cancer protein)
P04637  → p53 (Tumor suppressor)
P01308  → Insulin
P69905  → Hemoglobin subunit alpha
4. Open the 3D Viewer
After the script runs, it generates:
✅ PDB file: P38398.pdb
✅ 3D HTML viewer: P38398_structure.html

👉 Open the .html file in your browser to explore the interactive protein structure.
🔬 Example Output (P38398 – BRCA1)
Atoms: 14,550
Residues: 1,863
Chains: A
Domains: 5
Confidence (PLDDT):
Very High (>90): 200 residues
High (70–90): 125
Low/Very Low: 1,500+

✅ Why This Exists
This tool helps you:
Explore how AlphaFold models real biological proteins
Understand the structure–function relationship hands-on
Build your own dataset of folded proteins for downstream analysis or ML tasks
Whether you're a student, biologist, or AI researcher — this is your entry point into computational biology.
