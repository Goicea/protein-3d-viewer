# Protein 3D Viewer

> Aplicație interactivă pentru vizualizarea structurilor 3D ale proteinelor și acizilor nucleici.

![Python](https://img.shields.io/badge/Python-3.9-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-1.0-red)
![Biopython](https://img.shields.io/badge/Biopython-3D-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

## Descriere
**Vizualizatorul 3D de proteine** este o aplicație web care permite descărcarea, analiza și vzualizarea structurilor conținute în fișierele PDB, adică proteine și acizi nucleici.
Structurile sunt preluate din baza de date **RCSB Protein Data Bank**, parse-uite folosind **Biopython** și afișate folosind **Py3Dmol**.
Interfața aplicației este construită cu **Streamlit**.

## Interfață
![](C:\Users\Matei\Desktop\Proiect info vara - proteine\Documentatie\1.png)
![](C:\Users\Matei\Desktop\Proiect info vara - proteine\Documentatie\2.png)
![](C:\Users\Matei\Desktop\Proiect info vara - proteine\Documentatie\3.png)

## Funcționalități

- Descărcarea structurilor folosind PDB ID
- Parsarea fișierelor structurale
- Numărarea atomilor, reziduurilor și lanțurilor
- Vizualizare 3D interactivă
- Identificarea legăturilor de hidrogen
- Identificarea punților disulfidice

## Structură și biblioteci necesare
Codul aplicației e scris în mai multe fișiere pentru că organizarea asta permite să fie scris mai curat și o înțelegere mai simplă.

-pdb_utils.py - funcții pentru descărcat, parsat și extras informația din fișiere .pdb
-statistics.py - gestionează toate calculele pentru statistici, de la numărul atomilor la condițiile pe care trebuie să le îndeplinească o legătură de hidrogen validă
-vizualizator.py - creează fereastra de vizualizare, randează structura atomică și o afișează; se ocupă și cu alegerea unui mod de vizualizare și a colorării structurii
-app.py - integrează funcțiile din celelalte fisiere și construieste interfața web cu *streamlit* 

Pentru că întreg proiectul e scris în python folosirea unor biblioteci diverse și eficiente e foarte confortabilă. E posibilă și dezvoltarea aplicației în viitor folosind alte biblioteci sau funcții noi care apar în cele deja folosite, pentru că la python comunitatea este foarte activă și deschisă să ajute alți membri cu muncă din trecut.

| Bibliotecă suplimentară | Utilizare |
| Biopython | Parsarea și referirea la orice obiect din structuri |
| Streamlit | Interfața web |
| Py3Dmol | Generarea structurii 3D și funcții de zoom, rotire |
| stmol | Afișarea și integrarea Py3Dmol cu Streamlit |
| NumPy | Calcule matematice, în special vectoriale |

## Instalare
Clonează repository-ul:
```bash
git clone ADRESA_REPOSITORY
```
Intră în folder:
```bash
cd protein_viewer
```
Instalează dependențele:
```bash
pip install -r requirements.txt
```

## Rulare
Pornește aplicația folosind:
```bash
streamlit run app.py
```
Sau deschide launcherul batch file.

## Limitări

- Structurile foarte mari pot necesita mai mult timp pentru procesare.
- Randarea 3D depinde de performanța browserului și a sistemului.
- Calcularea interacțiunilor moleculare poate deveni costisitoare
  pentru structuri cu un număr foarte mare de atomi.

  ## Îmbunătățiri viitoare

- Optimizarea detectării legăturilor de hidrogen pentru structuri foarte mari
- Suport complet pentru formatul PDBx/mmCIF
- Selectarea individuală a reziduurilor în vizualizator
- Evidențierea legăturilor de hidrogen în reprezentarea 3D
