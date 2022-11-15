#!/usr/bin/env python3

# space separate tokens found in DeepSMILES

import argparse, deepsmiles, random, rdkit, re, sys, time

from rdkit import Chem

# DeepSMILES: no rings neither branches opening/closing
to_deep_smiles = deepsmiles.Converter(rings=True, branches=True)

# space-separate all DeepSMILES tokens corresponding to given SMILES
def tokenize_one(smi):
    assert(smi.find('.') == -1) # enforce standardization/salt removal
    mol = Chem.MolFromSmiles(smi)
    # don't canonicalize: the input SMILES might have been randomized on purpose
    protected_smi = Chem.MolToSmiles(mol, allHsExplicit=True, canonical=False)
    protected_dsmi = to_deep_smiles.encode(protected_smi)
    # print("pdsmi: '%s'" % protected_dsmi)
    # space before [ space after ]
    pdsmi = re.sub(r"(\[[^\]]+\])", r" \1 ", protected_dsmi)
    # space before %
    pdsmi = pdsmi.replace('%', ' %')
    # protect branch closings (no branch openings in DeepSMILES)
    pdsmi = pdsmi.replace(')', ' ) ')
    # protect bonds
    pdsmi = pdsmi.replace('-', ' - ')
    # protect - when it is a formal charge
    pdsmi = re.sub(' - (\d)\]', '-\1]', pdsmi)
    pdsmi = re.sub(' - \]', '-]', pdsmi)
    pdsmi = pdsmi.replace('=', ' = ')
    pdsmi = pdsmi.replace('#', ' # ')
    pdsmi = pdsmi.replace('$', ' $ ')
    pdsmi = pdsmi.replace(':', ' : ')
    # protect long numbers (prefixed by % in SMILES)
    pdsmi = re.sub(r"%(\d)(\d)", r" %\1\2 ", pdsmi)
    # single digit numbers are separate words
    pdsmi = re.sub(r" (\d)(\d)", r" \1 \2", pdsmi)
    # protect stereo bonds
    pdsmi = pdsmi.replace("/", " / ")
    pdsmi = pdsmi.replace("\\", " \\ ")
    # several spaces to one
    pdsmi = re.sub('[ ]+', ' ', pdsmi)
    # rm leading/trailing whitespaces
    pdsmi = pdsmi.strip()
    # print("pdsmi: '%s'" % pdsmi)
    return pdsmi

def random_reorder_atoms(mol):
    rand_order = list(range(mol.GetNumAtoms()))
    random.shuffle(rand_order)
    rand_mol = Chem.RenumberAtoms(mol, newOrder=rand_order)
    return rand_mol

# return n random versions of smi
def smi_randomize(smi, n):
    res = []
    mol = Chem.MolFromSmiles(smi)
    for i in range(n):
        rand_mol = random_reorder_atoms(mol)
        rand_smi = Chem.MolToSmiles(rand_mol, canonical=False)
        res.append(rand_smi)
    return res

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "DeepSMILES tokenizer")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output_std.smi", dest = "output_fn",
                        help = "molecules output file")
    parser.add_argument("--seed", dest = "rng_seed", default = 12345,
                        type = int, help = "RNG seed")
    parser.add_argument("--rand", dest = "randomize", default = False,
                        action = "store_true",
                        help = "randomize SMILES")
    parser.add_argument("-n", dest = "data_augment", default = 1,
                        type = int, help = "number of randomized SMILES output \
                        per input one")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    randomize = args.randomize
    augment_n = args.data_augment
    rng_seed = args.rng_seed
    # parse CLI end -----------------------------------------------------------
    random.seed(rng_seed)
    with open(output_fn, 'w') as output:
        for line in open(input_fn, 'r').readlines():
            smi_name = line.strip().split()
            # print("stripped: '%s'" % stripped)
            smi = smi_name[0]
            # print("smi: '%s'" % smi)
            name = smi_name[1]
            # print("name: '%s'" % name)
            if randomize or augment_n > 1:
                rand_smiles = smi_randomize(smi, augment_n)
                for rand_smi in rand_smiles:
                    tokenized = tokenize_one(rand_smi)
                    print(tokenized, file=output)
            else:
                tokenized = tokenize_one(smi)
                print(tokenized, file=output)
