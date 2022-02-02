Hierarchy iterators are ignored when converting to other packages, only the info in transferred

Molecule conversion spec

Example hierarchy conversion table
| At idx                     | 0   | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10  | 11  | 12  | 13  | 14  | 15  | 16 | 17  | 18    | 19    | 20    | 21    |
|----------------------------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|----|-----|-------|-------|-------|-------|
| OFFMol                     | 0   | 0   | 1   | 1   | 1   | 1   | 2   | 3   | 4   | 4   | 4   | 4   | 5   | 6   | 7   | 7   | 8  | 9   | 10    | 11    | 11    | 11    |
| resname                    | AAA | AAA | AAA | BBB | BBB | BBB | CCC | CCC | CCC | DDD | EEE | EEE | FFF | GGG | GGG | HHH | YZ | AAA | unset | unset | AAA   | unset |
| resnum                     | 1   | 1   | 1   | 1   | 2   | 2   | 2   | 3   | 3   | 4   | 4   | 5   | 6   | 6   | 7   | 8   | 9  | 1   | unset | unset | unset | unset |
| chain                      | A   | A   | A   | A   | A   | B   | B   | B   | C   | C   | D   | E   | E   | F   | G   | H   | I  | A   | unset | unset | unset | A     |
| Omm chain                  |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
| Omm res                    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
| OEHierChain                |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
| OEHierFragment             |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
| OEHierResidue              |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
| Chem.SplitMolByPDBChainId  |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
| Chem.SplitMolByPDBResidues |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
|                            |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
|                            |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |
|                            |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    |     |       |       |       |       |