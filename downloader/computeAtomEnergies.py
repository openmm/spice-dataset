import psi4
from openmm.app import element

# This script uses Psi4 to compute the per-atom reference energies that appear in the downloader script.
# It is included here in case we want to add more elements in the future.

charges = {'Br': (-1, 0), 'C': (-1, 0, 1), 'Ca': (2,), 'Cl': (-1, 0), 'F': (-1, 0), 'H': (-1, 0, 1), 'I': (-1, 0),
           'K': (1,), 'Li': (1,), 'Mg': (2,), 'N': (-1, 0, 1), 'Na': (1,), 'O': (-1, 0, 1), 'P': (0, 1), 'S': (-1, 0, 1)}

psi4.set_options({'reference': 'uhf'})
energies = {}
for symbol in charges:
    energies[symbol] = {}
    for charge in charges[symbol]:
        # Try all multiplicities up to 4 to find the one with lowest energy.

        electrons = element.get_by_symbol(symbol).atomic_number+charge
        multiplicity = 1 if electrons%2 == 0 else 2
        best_energy = None
        while multiplicity <= 4:
            try:
                mol = psi4.geometry(f"""
                {charge} {multiplicity}
                {symbol}
                """)
                energy = psi4.energy('wB97M-D3BJ/def2-TZVPPD', molecule=mol)
                if best_energy is None or energy < best_energy:
                    best_energy = energy
                    best_multiplicity = multiplicity
            except:
                pass
            multiplicity += 2
        energies[symbol][charge] = best_energy
print(energies)
