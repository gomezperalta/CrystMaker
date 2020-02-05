# CrystMaker
This is a python program to generate xyz-files of crystal compounds. As input data you only need a cif file. 

The code was done in Linux and uses the Python module 'pymatgen'. Please, check in the pymatgen documentation if this module is available for MacOs or Windows.

The program has three options:

<ul>
  <li>Positions: The atoms in the Unit Cell are generated without those atoms repeated in faces, edges or verteces.</li>
  <li>UnitCell: The atoms in the Unit Cell are generated, including all the repeated atoms in faces, edges or verteces.</li>
  <li>CrystMaker: You create a Crystal nxnxn times the Unit Cell.</li>
</ul>

Drawbacks:
<ol>
  <li>In case there were vacancies in the compound, the program will draw generate all the atoms ignoring the existing vacancies.</li>
  <li>In case there were two atom types in the same Wyckoff site, the program crashes. This normally happens with all programs when they deal with solid solutions.</li>
</ol>
