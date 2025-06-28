import requests
import py3Dmol  # pip install py3Dmol
import math
from collections import defaultdict

def fetch_pdb(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Structure not found: {uniprot_id}")

def save_structure_html(pdb_data, uniprot_id, filename=None):
    if filename is None:
        filename = f"{uniprot_id}_structure.html"
    
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(pdb_data, "pdb")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})
    viewer.zoomTo()
    html_content = viewer._make_html()
    
    with open(filename, 'w') as f:
        f.write(html_content)
    
    print(f"3D structure saved {filename}")
    return filename

def save_pdb_file(pdb_data, uniprot_id, filename=None):
    if filename is None:
        filename = f"{uniprot_id}.pdb"
    
    with open(filename, 'w') as f:
        f.write(pdb_data)    
    print(f"PDB file saved {filename}")
    return filename

def analyze_structure_info(pdb_data, uniprot_id):
    lines = pdb_data.split('\n')
    atom_lines = [line for line in lines if line.startswith('ATOM')]
    num_atoms = len(atom_lines)
    header_lines = [line for line in lines if line.startswith('HEADER')]
    title_lines = [line for line in lines if line.startswith('TITLE')]
    
    print(f"\n STRUCTURE ANALYSIS for {uniprot_id}:")
    print("="*50)
    print(f"Number of atoms {num_atoms:,}")
    
    if header_lines:
        print(f"Header: {header_lines[0][10:].strip()}")
    
    if title_lines:
        for title_line in title_lines:
            print(f"Title: {title_line[10:].strip()}")
    
    chains = set()
    for line in atom_lines:
        if len(line) > 21:
            chain_id = line[21]
            chains.add(chain_id)    
    print(f"Chains found: {', '.join(sorted(chains))}")
    
    residues = set()
    for line in atom_lines:
        if len(line) > 26:
            res_num = line[22:26].strip()
            if res_num:
                residues.add(res_num)
    print(f"Approximate residues: {len(residues)}")

class ProteinAnalyzer:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        self.pdb_data = None
        self.atoms = []
        
    def fetch_structure(self):
        """Download structure from AlphaFold"""
        url = f"https://alphafold.ebi.ac.uk/files/AF-{self.uniprot_id}-F1-model_v4.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            self.pdb_data = response.text
            self._parse_atoms()
            return True
        return False
    
    def _parse_atoms(self):
        """Parse PDB data into structured format"""
        self.atoms = []
        for line in self.pdb_data.split('\n'):
            if line.startswith('ATOM'):
                atom = {
                    'atom_num': int(line[6:11]),
                    'atom_name': line[12:16].strip(),
                    'residue': line[17:20].strip(),
                    'chain': line[21],
                    'res_num': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'confidence': float(line[60:66])
                }
                self.atoms.append(atom)
    
    def get_confidence_regions(self):
        """Analyze AlphaFold confidence scores"""
        confidence_ranges = {
            'Very High (>90)': [],
            'High (70-90)': [],
            'Low (50-70)': [],
            'Very Low (<50)': []
        }
        
        for atom in self.atoms:
            conf = atom['confidence']
            if conf > 90:
                confidence_ranges['Very High (>90)'].append(atom['res_num'])
            elif conf > 70:
                confidence_ranges['High (70-90)'].append(atom['res_num'])
            elif conf > 50:
                confidence_ranges['Low (50-70)'].append(atom['res_num'])
            else:
                confidence_ranges['Very Low (<50)'].append(atom['res_num'])
        
        return confidence_ranges
    
    def find_domains(self):
        high_conf_residues = set()
        for atom in self.atoms:
            if atom['confidence'] > 80 and atom['atom_name'] == 'CA':
                high_conf_residues.add(atom['res_num'])
        
        sorted_residues = sorted(high_conf_residues)
        domains = []
        current_domain = [sorted_residues[0]] if sorted_residues else []
        
        for i in range(1, len(sorted_residues)):
            if sorted_residues[i] - sorted_residues[i-1] <= 3: 
                current_domain.append(sorted_residues[i])
            else:
                if len(current_domain) >= 20: 
                    domains.append((min(current_domain), max(current_domain)))
                current_domain = [sorted_residues[i]]
        
        if len(current_domain) >= 20:
            domains.append((min(current_domain), max(current_domain)))
        
        return domains

def main():
    # Example UniProt IDs to try
    examples = {
        "P38398": "BRCA1 (Breast cancer 1 protein)",
        "P04637": "p53 (Tumor protein p53)", 
        "P01308": "Insulin",
        "P69905": "Hemoglobin subunit alpha",
        "P02768": "Human serum albumin"
    }
    
    print("ALPHAFOLD STRUCTURE DOWNLOADER")
    print("="*40)
    print("\nAvailable examples:")
    for uid, name in examples.items():
        print(f"  {uid}: {name}")
    uniprot_id = "P38398" # You can change this to any UniProt ID
    protein_name = examples.get(uniprot_id, "Unknown protein")
    print(f"\nFetching {uniprot_id} ({protein_name})...")
    
    try:
        analyzer = ProteinAnalyzer(uniprot_id)
        if analyzer.fetch_structure():
            print(f"Successfully downloaded structure for {uniprot_id}")
            analyze_structure_info(analyzer.pdb_data, uniprot_id)
            html_file = save_structure_html(analyzer.pdb_data, uniprot_id)
            pdb_file = save_pdb_file(analyzer.pdb_data, uniprot_id)
            print(f"• {html_file} (3D visualization)")
            print(f"• {pdb_file} (PDB structure file)")
            
            conf_regions = analyzer.get_confidence_regions()
            domains = analyzer.find_domains()
            
            print("Confidence Analysis:")
            for region, residues in conf_regions.items():
                print(f"{region}: {len(set(residues))} residues")
            
            print("Predicted Domains:")
            for i, (start, end) in enumerate(domains, 1):
                print(f"Domain {i}: residues {start}-{end} ({end-start+1} residues)")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
