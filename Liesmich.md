# SpatialData Zelltyp-Annotation - Kurze Anleitung

Hier ist eine schnelle Anleitung, wie du das Zelltyp-Annotations-Tool verwendest und mit Git arbeitest.

## Das Tool: celltype_annotate_cli_v2.py

### Was macht es?
Das Tool nimmt deine Xenium-Daten und annotiert die Zelltypen mit Hilfe von Single-Cell-Referenzdaten. Am Ende bekommst du:
- Annotierte Spatial-Daten als `.zarr` Datei
- Schöne Visualisierungen (räumliche Zelltyp-Karten, UMAP-Plots, etc.)
- Einen HTML-Report mit allen Ergebnissen
- Ordentlich organisierte Ausgabe-Ordner

### Schnellstart
```bash
# Einfachste Version
uv run python celltype_annotate_cli_v2.py meine_xenium_daten.zarr meine_referenz.h5ad

# Mit mehr Optionen
uv run python celltype_annotate_cli_v2.py meine_xenium_daten.zarr meine_referenz.h5ad \
    --results-dir /pfad/zu/meinen/ergebnissen \
    --min-clusters 5 \
    --max-clusters 15 \
    --show-unknown-cells
```

### Wichtige Optionen

**Ausgabe-Verwaltung:**
- `--results-dir /pfad/zur/ausgabe` - Wo die Ergebnisse gespeichert werden sollen
- `--consolidate-data` - Mehr Infos über Datenportabilität (für später nützlich)

**Visualisierung:**
- `--show-unknown-cells` - Zeigt auch unbekannte Zelltypen in den räumlichen Plots
- Standardmäßig werden unbekannte Zellen ausgeblendet für sauberere Bilder

**Clustering:**
- `--min-clusters 5` - Mindestanzahl Cluster
- `--max-clusters 15` - Höchstanzahl Cluster

**Andere nützliche Optionen:**
- `--overwrite` - Überschreibt existierende Ausgaben
- `--log-level debug` - Mehr Details im Log
- `--help` - Zeigt alle verfügbaren Optionen

### Ausgabe-Struktur
Das Tool erstellt Ordner mit dem Format:
```
ergebnisse/xenium_dateiname___referenz_dateiname___20250724_120628/
├── data/                    # Annotierte .zarr Datei
├── plots/                   # Alle Visualisierungen
├── logs/                    # Log-Dateien
└── annotation_report.html  # Schöner HTML-Report
```

## Git-Workflow für das Projekt

### Projekt auschecken
```bash
# Repository klonen
git clone https://github.com/chrismuke/spatialdata-notebooks.git
cd spatialdata-notebooks

# Auf den richtigen Branch wechseln
git checkout feature/celltype-annotation-cli

# Dependencies installieren
uv sync
```

### Branch-Management

**Neuen Branch erstellen:**
```bash
# Neuen Branch von aktueller Position erstellen
git checkout -b mein-neuer-feature-branch

# Oder von einem bestimmten Branch
git checkout -b mein-feature feature/celltype-annotation-cli
```

**Branch wechseln:**
```bash
# Zu existierendem Branch wechseln
git checkout feature/celltype-annotation-cli
git checkout main

# Alle verfügbaren Branches anzeigen
git branch -a
```

**Branches mergen:**
```bash
# Wechsel zum Ziel-Branch (meist main)
git checkout main

# Branch mergen
git merge mein-feature-branch

# Oder interaktiver merge
git merge --no-ff mein-feature-branch
```

### Änderungen committen und pushen

**Basis-Workflow:**
```bash
# Status checken
git status

# Dateien zur Staging Area hinzufügen
git add .                    # Alle Änderungen
git add datei1.py datei2.md  # Spezifische Dateien

# Commit erstellen
git commit -m "Kurze Beschreibung der Änderungen"

# Zu GitHub pushen
git push origin mein-branch-name
```

**Für längere Commit-Messages:**
```bash
git commit -m "Kurze Zusammenfassung

Längere Beschreibung was gemacht wurde:
- Feature X hinzugefügt
- Bug Y behoben
- Dokumentation erweitert"
```

### GitHub Desktop Alternative

Falls du nicht so gerne mit der Kommandozeile arbeitest, ist **GitHub Desktop** super praktisch:

1. **Download:** https://desktop.github.com/
2. **Repository klonen:** File → Clone Repository → URL eingeben
3. **Branches:** Current Branch Dropdown → New Branch / Switch Branch
4. **Commits:** Änderungen auswählen → Summary eingeben → Commit
5. **Push:** "Push origin" Button klicken
6. **Pull Requests:** Branch → Create Pull Request

GitHub Desktop zeigt auch schön die Diffs an und macht Merging einfacher.

### Typischer Arbeitsablauf

1. **Feature-Branch erstellen:**
   ```bash
   git checkout -b neue-funktion
   ```

2. **Code schreiben und testen:**
   ```bash
   # Änderungen machen...
   uv run python celltype_annotate_cli_v2.py --help  # testen
   ```

3. **Committen:**
   ```bash
   git add .
   git commit -m "Neue Funktion XYZ hinzugefügt"
   ```

4. **Pushen:**
   ```bash
   git push origin neue-funktion
   ```

5. **Pull Request auf GitHub erstellen** (über die Web-UI oder GitHub Desktop)

6. **Nach Review mergen** und lokalen Branch aufräumen:
   ```bash
   git checkout main
   git pull origin main
   git branch -d neue-funktion  # lokalen Branch löschen
   ```

## Häufige Probleme

**"Not self-contained" Warnungen:**
Das ist normal! Das Tool spart Speicherplatz, indem es große Bilddateien nicht kopiert. Deine Analyse funktioniert trotzdem perfekt.

**Spatial Celltype Map ist leer:**
Das ist jetzt behoben! Das Tool filtert standardmäßig unbekannte Zellen raus. Mit `--show-unknown-cells` kannst du alle sehen.

**Git merge conflicts:**
Am besten mit GitHub Desktop oder einem Merge-Tool wie VSCode lösen. Oder frag einfach! 😊

## Hilfe

- **Tool-Hilfe:** `uv run python celltype_annotate_cli_v2.py --help`
- **Git-Hilfe:** `git help <command>` oder Google
- **GitHub Desktop:** Hat eine gute integrierte Hilfe
- **Bei Problemen:** Einfach fragen oder Issue auf GitHub erstellen
