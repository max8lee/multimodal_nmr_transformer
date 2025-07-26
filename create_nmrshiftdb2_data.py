import xml.etree.ElementTree as ET
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional, Set
import re
from sklearn.model_selection import train_test_split
from pathlib import Path

# RDKit import with fallback
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Install with: conda install -c conda-forge rdkit")

class InvalidTokenError(Exception):
    """Raised when invalid tokens are found during tokenization."""
    def __init__(self, invalid_tokens: List[str], smiles: str):
        self.invalid_tokens = invalid_tokens
        self.smiles = smiles
        super().__init__(f"Invalid tokens found in SMILES '{smiles}': {invalid_tokens}")

def load_vocab(vocab_file: str) -> Set[str]:
    """Load vocabulary from a vocab file and return set of valid tokens."""
    vocab = set()
    with open(vocab_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                # Split on whitespace and take the first part as the token
                token = line.split()[0]
                vocab.add(token)
    return vocab

def tokenize_smiles_vocab_aligned(smiles: str, vocab_set: Optional[Set[str]] = None):
    pattern = r"""
        \[[^\[\]]+\]         |  # bracketed expressions
        Cl|Br                |  # two-letter atoms
        B|C|N|O|S|P|F|I|H    |  # one-letter atoms
        b|c|n|o|s|p          |  # aromatic atoms
        \(|\)|=|\#|\+|-|\\|/ |  # bonds and branches
        [0-9]                   # ring closure digits
    """
    regex = re.compile(pattern, re.X)
    tokens = regex.findall(smiles)
    invalid_tokens = [token for token in tokens if token not in vocab_set]
    
    if invalid_tokens:
        raise InvalidTokenError(invalid_tokens, smiles)
    return tokens

def cml_to_rdkit_molecule_3d(molecule, ns):
    """Convert CML molecule (atoms + bonds) to RDKit molecule object with stereochemistry"""
    if not RDKIT_AVAILABLE:
        return None
    
    mol = Chem.RWMol()
    atom_id_to_idx = {}
    
    # Add atoms
    atom_array = molecule.find('.//cml:atomArray', ns)
    if atom_array is None:
        return None
        
    atoms = atom_array.findall('.//cml:atom', ns)
    if not atoms:
        return None
    
    for atom in atoms:
        atom_id = atom.get('id')
        element_type = atom.get('elementType')
        formal_charge = int(atom.get('formalCharge', '0'))
        
        if not atom_id or not element_type:
            continue
            
        # Create RDKit atom
        rdkit_atom = Chem.Atom(element_type)
        rdkit_atom.SetFormalCharge(formal_charge)
        
        atom_idx = mol.AddAtom(rdkit_atom)
        atom_id_to_idx[atom_id] = atom_idx
    
    # Add bonds with stereochemistry information
    bond_array = molecule.find('.//cml:bondArray', ns)
    if bond_array is not None:
        bonds = bond_array.findall('.//cml:bond', ns)
        
        for bond in bonds:
            atom_refs = bond.get('atomRefs2')
            bond_order = bond.get('order', 'S')
            
            if not atom_refs:
                continue
                
            atom_ids = atom_refs.strip().split()
            if len(atom_ids) != 2:
                continue
                
            atom1_id, atom2_id = atom_ids
            
            if atom1_id not in atom_id_to_idx or atom2_id not in atom_id_to_idx:
                continue
                
            atom1_idx = atom_id_to_idx[atom1_id]
            atom2_idx = atom_id_to_idx[atom2_id]
            
            # Convert bond order
            bond_type_map = {
                'S': Chem.BondType.SINGLE, '1': Chem.BondType.SINGLE,
                'D': Chem.BondType.DOUBLE, '2': Chem.BondType.DOUBLE,
                'T': Chem.BondType.TRIPLE, '3': Chem.BondType.TRIPLE,
                'A': Chem.BondType.AROMATIC
            }
            rdkit_bond_type = bond_type_map.get(bond_order, Chem.BondType.SINGLE)
            
            mol.AddBond(atom1_idx, atom2_idx, rdkit_bond_type)
            bond_idx = mol.GetNumBonds() - 1

            bond_stereo_elem = bond.find('.//cml:bondStereo', ns)
            if bond_stereo_elem is not None:
                bond_stereo = bond_stereo_elem.text.strip() if bond_stereo_elem.text else None
                
                if bond_stereo:
                    rdkit_bond = mol.GetBondWithIdx(bond_idx)
                    if bond_stereo == 'W':  # Wedge (up)
                        rdkit_bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                        mol.GetAtomWithIdx(atom1_idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                    elif bond_stereo == 'H':  # Hatch/Dash (down)
                        rdkit_bond.SetBondDir(Chem.BondDir.BEGINDASH)
                        mol.GetAtomWithIdx(atom1_idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

    mol = mol.GetMol()
    if mol is None:
        return None
    
    try:
        # Assign stereochemistry from 3D coordinates if available
        before_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        #print(f"Chiral centers before: {before_centers}")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)                
        # Alternative approach for better stereo assignment
        Chem.AssignAtomChiralTagsFromStructure(mol)
        after_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        #print(f"Chiral centers after: {after_centers}")                
    except:
        print("Error assigning stereochemistry")
        pass
    return mol

def extract_smiles_from_structure(molecule, ns, vocab_set) -> Optional[str]:
    """Extract SMILES from molecular structure"""
    if not RDKIT_AVAILABLE:
        return None
        
    rdkit_mol = cml_to_rdkit_molecule_3d(molecule, ns)
    if rdkit_mol is None:
        return None

    try:
        rdkit_mol = Chem.RemoveHs(rdkit_mol)
    except:
        print("Skip: H removal error")
        return None
    try:
        smiles = Chem.MolToSmiles(rdkit_mol, isomericSmiles=True, canonical=True)
        tokens = tokenize_smiles_vocab_aligned(smiles, vocab_set)
        return ' '.join(tokens)
    except:
        print("Skip: " + Chem.MolToSmiles(rdkit_mol) + "contains excluded element")
        return None

def extract_molecular_formula(molecule, ns) -> Optional[str]:
    """Extract molecular formula from molecule"""
    atom_counts = defaultdict(int)
    
    atom_array = molecule.find('.//cml:atomArray', ns)
    if atom_array is None:
        return None
        
    atoms = atom_array.findall('.//cml:atom', ns)
    for atom in atoms:
        element = atom.get('elementType')
        if element:
            atom_counts[element] += 1
    
    if not atom_counts:
        return None
    
    # Format as molecular formula (C, H, then alphabetical)
    formula_parts = []
    
    # Carbon first
    if 'C' in atom_counts:
        count = atom_counts['C']
        formula_parts.append(f"C {count}" if count > 1 else "C")
        del atom_counts['C']
    
    # Hydrogen second
    if 'H' in atom_counts:
        count = atom_counts['H']
        formula_parts.append(f"H {count}" if count > 1 else "H")
        del atom_counts['H']
    
    # Others alphabetically
    for element in sorted(atom_counts.keys()):
        count = atom_counts[element]
        formula_parts.append(f"{element} {count}" if count > 1 else element)
    
    return " ".join(formula_parts)

def extract_nmr_data(spectrum, ns, vocab_tokens) -> Dict:
    """Extract NMR peak data from spectrum"""
    nmr_data = {'nucleus': None, 'peaks': []}
    
    # Get observed nucleus
    metadata_list = spectrum.find('.//cml:metadataList', ns)
    if metadata_list is not None:
        for metadata in metadata_list.findall('.//cml:metadata', ns):
            name = metadata.get('name', '')
            if 'OBSERVENUCLEUS' in name:
                nmr_data['nucleus'] = metadata.get('content', '')
                break
    
    # Get peaks
    peak_list = spectrum.find('.//cml:peakList', ns)
    if peak_list is not None:
        peaks = peak_list.findall('.//cml:peak', ns)
        
        for peak in peaks:
            x_value = peak.get('xValue')
            if not x_value:
                continue 
            try:
                x_val_float = float(x_value)
            except:
                continue

            atom_refs = peak.get('atomRefs', '')
            atom_count = len(atom_refs.split()) if atom_refs else 0
            atom_list = atom_refs.strip().split() if atom_refs else []

            multiplicity = peak.get('peakMultiplicity').lower()
            if not multiplicity or multiplicity not in vocab_tokens:
                multiplicity = 's' if atom_count <= 1 else 'm' # Create fake multiplicity when necessary

            peak_data = {
                'xValue': x_val_float,
                'multiplicity': multiplicity,
                'atomRefs': atom_list,
                'height': peak.get('peakHeight', '1'),
                'coupling': []
            }

            seen_jvalues = set()
            for struct in peak.findall('.//cml:peakStructure', ns):
                if struct.get('type') == 'coupling':
                    val = struct.get('value')
                    if val is not None:
                        try:
                            j = float(val)
                            seen_jvalues.add(j)
                        except ValueError:
                            continue
            peak_data["coupling"] = sorted(seen_jvalues)
            nmr_data['peaks'].append(peak_data)
    
    return nmr_data if nmr_data['peaks'] else None

def format_1h_nmr(peaks: List[Dict]) -> str:
    """Format 1H NMR peaks in the required format"""
    if not peaks:
        return ""
    nmr_parts = []
    for peak in peaks:
        ppm = peak['xValue']
        multiplicity = peak.get('multiplicity', 's')
        integration = len(peak.get('atomRefs', [])) or 1
        coupling = peak.get('coupling', [])
        ppm_min = ppm - 0.02 # Creating fake peak start
        ppm_max = ppm + 0.02 # Creating fake peak end

        peak_str = f"{ppm_max:.2f} {ppm_min:.2f} {multiplicity} {integration}H"
        if coupling:
            j_str = " ".join(f"{j:.2f}" for j in coupling)
            peak_str += f" J {j_str}"
        nmr_parts.append(peak_str)
    return "1HNMR " + " | ".join(nmr_parts)

def format_13c_nmr(peaks: List[Dict]) -> str:
    """Format 13C NMR peaks in the required format"""
    if not peaks:
        return ""
    # Sort peaks by chemical shift (descending)
    sorted_peaks = sorted(peaks, key=lambda p: p['xValue'], reverse=True)
    shifts = []
    for peak in sorted_peaks:
        shift = round(peak['xValue'], 1)
        shifts.append(str(shift))
    return f"13CNMR {' '.join(shifts)}"

def format_source_line(mol_formula: str, nmr_data_list: List[Dict]) -> str:
    """Format the source line combining multiple NMR types"""
    parts = [mol_formula]
    
    # Separate 1H and 13C data
    h_data = []
    c_data = []

    for nmr_data in nmr_data_list:
        nucleus = nmr_data.get('nucleus', '')
        if nucleus == '1H' and not h_data:  # Only take the first 1H spectrum. Can be changed to take average spectra
            h_data = nmr_data['peaks']
        elif nucleus == '13C' and not c_data:  # Only take the first 13C spectrum. Can be changed to take average spectra
            c_data = nmr_data['peaks']
        else:
            print("Skip: " + nucleus)
    
    # Add 1H NMR if available
    if h_data:
        h_nmr = format_1h_nmr(h_data)
    else:
        h_nmr = "1HNMR"
    parts.append(h_nmr)
    
    # Add 13C NMR if available
    if c_data:
        c_nmr = format_13c_nmr(c_data)
    else:
        c_nmr = "13CNMR"
    parts.append(c_nmr)
    
    return " ".join(parts)

def extract_nmr_dataset(file_path: str, seed: str, out_path: str, 
                       max_records: Optional[int] = None, verbose: bool = True, vocab_tokens: set = None):
    """
    Main function to extract NMR dataset from NMRShiftDB2 XML
    """
    
    if not RDKIT_AVAILABLE:
        print("ERROR: RDKit is required for SMILES generation!")
        print("Install with: conda install -c conda-forge rdkit")
        return
    
    tree = ET.parse(file_path)
    root = tree.getroot()
    ns = {'cml': 'http://www.xml-cml.org/schema'}
    
    if verbose:
        print(f"Processing {len(root)} elements...")
    
    # Create molecule lookup
    molecules = {}
    for mol in root.findall('.//cml:molecule', ns):
        mol_id = mol.get('id')
        if mol_id:
            molecules[mol_id] = mol
    
    if verbose:
        print(f"Found {len(molecules)} unique molecules")
    
    # Group spectra by molecule
    molecule_spectra = defaultdict(list)
    for spectrum in root.findall('.//cml:spectrum', ns):
        molecule_id = spectrum.get('moleculeRef')
        if molecule_id and molecule_id in molecules:
            molecule_spectra[molecule_id].append(spectrum)
    
    if verbose:
        print(f"Found {len(molecule_spectra)} molecules with spectra")
    
    # Process each molecule
    src_lines = []
    tgt_lines = []
    processed_count = 0
    skipped_count = 0
    
    for mol_id, spectra in molecule_spectra.items():
        if max_records and processed_count >= max_records:
            break
            
        molecule = molecules[mol_id]
        
        try:
            # Extract molecular formula
            mol_formula = extract_molecular_formula(molecule, ns)
            if not mol_formula:
                skipped_count += 1
                continue
            
            # Extract SMILES
            vocab_set = load_vocab("runs/runs_f_groups/all_modalities_v2/data/vocab/vocab.tgt")
            smiles = extract_smiles_from_structure(molecule, ns, vocab_set)
            if not smiles:
                skipped_count += 1
                continue
            
            # Extract all NMR data for this molecule
            nmr_data_list = []
            for spectrum in spectra:
                nmr_data = extract_nmr_data(spectrum, ns, vocab_tokens)
                if nmr_data:
                    nmr_data_list.append(nmr_data)
            if not nmr_data_list:
                skipped_count += 1
                continue
            has_valid_nucleus = any(
                nmr_data.get('nucleus', '') in ['1H', '13C'] for nmr_data in nmr_data_list
            )
            if not has_valid_nucleus:
                skipped_count += 1
                continue
            
            # Format source line
            src_line = format_source_line(mol_formula, nmr_data_list)
            
            # Add to outputs
            src_lines.append(src_line)
            tgt_lines.append(smiles)
            processed_count += 1
            
            if verbose and processed_count % 100 == 0:
                print(f"Processed {processed_count} molecules...")
                
        except Exception as e:
            if verbose:
                print(f"Error processing molecule {mol_id}: {e}")
            skipped_count += 1
            continue
    
    # Split
    src_train, src_test, tgt_train, tgt_test = train_test_split(
        src_lines, tgt_lines, test_size=0.1, random_state=seed, shuffle=True
    )

    src_train, src_val, tgt_train, tgt_val = train_test_split(
        src_train, tgt_train, test_size=0.05, random_state=seed, shuffle=True
    )

    # Save helper
    def write_lines(lines, path):
        with open(path, "w", encoding="utf-8") as f:
            f.write('\n'.join(lines))

    # Save all sets
    out_data_path = Path(out_path) / "data"
    out_data_path.mkdir(parents=True, exist_ok=True)

    write_lines(src_train, out_data_path / "src-train.txt")
    write_lines(tgt_train, out_data_path / "tgt-train.txt")
    write_lines(src_val, out_data_path / "src-val.txt")
    write_lines(tgt_val, out_data_path / "tgt-val.txt")
    write_lines(src_test, out_data_path / "src-test.txt")
    write_lines(tgt_test, out_data_path / "tgt-test.txt")
    
    if verbose:
        print(f"\n=== EXTRACTION COMPLETE ===")
        print(f"Successfully processed: {processed_count} molecules")
        print(f"Skipped: {skipped_count} molecules")
        print(f"Success rate: {processed_count/(processed_count+skipped_count)*100:.1f}%")
    
    return processed_count, skipped_count

# Usage
if __name__ == "__main__":
    file_path = "nmrshiftdb2.xml"
    
    try:
        with open("runs/h_and_c_nmr_v1/data/vocab/vocab.src", "r", encoding="utf-8") as f:
            vocab_tokens = set(line.strip().split('\t')[0] for line in f if line.strip())
        
        print("\nStarting extraction...")
        extract_nmr_dataset(
            file_path=file_path,
            seed=3245,
            out_path="runs/finetune_v1",
            max_records=10,
            verbose=True,
            vocab_tokens=vocab_tokens
        )
        
    except FileNotFoundError:
        print(f"File not found! Please check if you have nmrshiftdb2.xml and runs/h_and_c_nmr_v1/data/vocab/vocab.src")
    except Exception as e:
        print(f"Error: {e}")