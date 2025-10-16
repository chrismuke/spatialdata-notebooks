# Cell Type Explanations

This document provides biological explanations for all cell types found in the published annotations and predicted cell types from the Xenium Human Breast Cancer dataset.

---

## Published Cell Type Annotations (19 types)

### Tumor/Cancer Cells

#### **Invasive Tumor**
Malignant breast cancer cells that have invaded beyond the basement membrane into surrounding tissue. These are aggressive cancer cells capable of metastasis. Most abundant cell type in the dataset (4,568 cells).

#### **DCIS 1** (Ductal Carcinoma In Situ, Type 1)
Pre-invasive breast cancer cells confined within milk ducts. These are early-stage cancer cells that haven't broken through the duct walls. Subtype 1 classification based on molecular/spatial characteristics.

#### **DCIS 2** (Ductal Carcinoma In Situ, Type 2)
Similar to DCIS 1 but represents a different molecular or spatial subtype of ductal carcinoma in situ. May have different gene expression patterns or location within the tissue.

#### **Prolif Invasive Tumor** (Proliferative Invasive Tumor)
Actively dividing/proliferating invasive cancer cells. These tumor cells are in high mitotic state, indicating rapid growth and aggressive behavior.

### Epithelial Cells

#### **Myoepi ACTA2+** (Myoepithelial cells, ACTA2-positive)
Contractile epithelial cells that surround mammary ducts and alveoli. Express ACTA2 (alpha smooth muscle actin), indicating contractile properties. Help squeeze milk during lactation and provide structural support.

#### **Myoepi KRT15+** (Myoepithelial cells, KRT15-positive)
Myoepithelial cells marked by keratin 15 expression. Represent a stem/progenitor-like subpopulation of myoepithelial cells with potential regenerative capacity.

### Immune Cells - Myeloid Lineage

#### **Macrophages 1**
Tissue-resident immune cells that phagocytose (engulf) debris, pathogens, and dead cells. In tumors, can have pro-tumor or anti-tumor roles. Type 1 likely represents a specific activation state or location (2,964 cells).

#### **Macrophages 2**
Second subtype of macrophages with distinct molecular signature or spatial distribution. May represent differently activated macrophages (M1 vs M2 polarization) or location-specific populations.

#### **IRF7+ DCs** (IRF7-positive Dendritic Cells)
Specialized dendritic cells expressing Interferon Regulatory Factor 7, involved in antiviral responses and type I interferon production. Play key role in initiating adaptive immune responses.

#### **LAMP3+ DCs** (LAMP3-positive Dendritic Cells)
Mature, activated dendritic cells expressing LAMP3 (lysosomal-associated membrane protein 3). Often called "mature DCs" or "migratory DCs" - actively migrating to lymph nodes to present antigens to T cells.

#### **Mast Cells**
Immune cells containing granules filled with histamine and heparin. Involved in allergic responses and inflammation. In tumors, can influence angiogenesis and immune responses.

### Immune Cells - Lymphoid Lineage

#### **CD4+ T Cells**
Helper T cells that coordinate immune responses by activating other immune cells. Express CD4 surface marker. Critical for anti-tumor immunity and orchestrating adaptive immune responses (2,899 cells).

#### **CD8+ T Cells**
Cytotoxic T cells that directly kill infected or cancerous cells. Express CD8 surface marker. Major effector cells in anti-tumor immunity (1,843 cells).

#### **B Cells**
Lymphocytes that produce antibodies. Can differentiate into plasma cells. In tumors, form tertiary lymphoid structures and contribute to anti-tumor immunity.

### Stromal/Support Cells

#### **Stromal**
Connective tissue cells, primarily fibroblasts, that provide structural framework. Produce extracellular matrix proteins. In tumors, can become cancer-associated fibroblasts (CAFs) that support tumor growth (2,611 cells).

#### **Endothelial**
Cells lining blood vessels. Form the inner layer of blood vessel walls, controlling exchange of materials between blood and tissue. Essential for tumor angiogenesis.

#### **Perivascular-Like**
Cells surrounding blood vessels, likely including pericytes and smooth muscle cells. Provide structural support to blood vessels and regulate blood flow.

### Hybrid/Doublet Cells

#### **T Cell & Tumor Hybrid**
Likely technical artifact representing doublets (two cells captured together) - one T cell and one tumor cell. Could also represent very close spatial proximity or cell-cell interactions (1,003 cells).

#### **Stromal & T Cell Hybrid**
Similar doublet artifact with one stromal cell and one T cell captured together. Alternatively, represents extremely close spatial association between these cell types.

---

## Predicted Cell Type Annotations (20 types)

### Epithelial Cells

#### **keratinocyte** (66,586 cells - most abundant)
Epithelial cells that produce keratin proteins. In breast tissue context, this annotation likely represents epithelial cells including both normal ductal/lobular epithelium and transformed tumor cells. Vastly overrepresented compared to published annotations.

#### **melanocyte** (4,951 cells)
Pigment-producing cells normally found in skin. Unexpected in breast tissue - may represent misclassification or cross-tissue annotation from skin-trained reference.

#### **Langerhans cell** (8 cells)
Specialized dendritic cells typically found in skin. Act as antigen-presenting cells. Very rare in this dataset, possibly misclassified or contaminant cells.

### Fibroblast/Stromal Cells

#### **skin fibroblast** (38,960 cells - second most abundant)
Connective tissue cells from skin. In breast tissue, these likely represent general fibroblasts/stromal cells. The "skin" prefix suggests the reference dataset was skin-derived.

### Immune Cells - Myeloid Lineage

#### **macrophage** (20,499 cells)
Large phagocytic immune cells that engulf debris and pathogens. Can be pro-inflammatory (M1) or anti-inflammatory/pro-tumor (M2). Key components of tumor microenvironment.

#### **dendritic cell** (5,018 cells)
Professional antigen-presenting cells that bridge innate and adaptive immunity. Capture antigens and present them to T cells to initiate immune responses.

#### **monocyte-derived dendritic cell** (603 cells)
Dendritic cells differentiated from monocytes, typically during inflammation. Represent activated state in response to tissue damage or infection.

#### **conventional dendritic cell** (307 cells)
Classical dendritic cells that develop directly from bone marrow precursors (as opposed to monocyte-derived). Specialized in antigen presentation.

#### **inflammatory macrophage** (329 cells)
Activated macrophages in pro-inflammatory state (likely M1 polarization). Produce inflammatory cytokines and help fight infections/tumors.

#### **mast cell** (142 cells)
Granule-containing immune cells involved in allergic reactions and inflammation. Release histamine and other mediators.

### Immune Cells - Lymphoid Lineage

#### **helper T cell** (10,061 cells)
CD4+ T cells that coordinate immune responses. Activate other immune cells through cytokine secretion. Essential for adaptive immunity.

#### **cytotoxic T cell** (4,365 cells)
CD8+ T cells that directly kill target cells. Major effector cells against cancer and viral infections.

#### **regulatory T cell** (1,380 cells)
Specialized T cells (Tregs) that suppress immune responses to maintain homeostasis and prevent autoimmunity. In tumors, can suppress anti-tumor immunity.

#### **plasma cell** (856 cells)
Antibody-secreting cells differentiated from B cells. Produce large amounts of specific antibodies.

#### **natural killer cell** (679 cells)
Innate immune cells that kill virus-infected or cancerous cells without prior sensitization. Provide rapid immune response.

#### **innate lymphoid cell** (592 cells)
Recently discovered family of immune cells that resemble lymphocytes but lack antigen-specific receptors. Involved in tissue homeostasis and early immune responses.

### Vascular Cells

#### **endothelial cell of vascular tree** (7,433 cells)
Cells lining blood vessels. Form the inner layer of arteries, veins, and capillaries. Essential for angiogenesis and nutrient delivery.

#### **endothelial cell of lymphatic vessel** (311 cells)
Cells lining lymphatic vessels. Specialized for fluid drainage and immune cell trafficking. Distinct from blood vascular endothelium.

#### **pericyte** (4,267 cells)
Contractile cells that wrap around blood vessels. Provide structural support, regulate blood flow, and maintain blood-brain barrier. Important for vessel stability.

### Neural Cells

#### **Schwann cell** (433 cells)
Glial cells that produce myelin sheath around peripheral nerves. Unexpected in breast tissue - may represent nerve-associated cells or misclassification.

---

## Notes on Mapping Differences

The predicted cell types appear to be based on a **skin-derived reference atlas**, explaining annotations like "skin fibroblast," "melanocyte," "Langerhans cell," and "Schwann cell" that are atypical for breast tissue.

The published annotations are **breast cancer-specific** with detailed tumor cell classifications (DCIS subtypes, invasive tumor, proliferative states) that are clinically relevant.

**Key discrepancies:**
- Predicted model groups all epithelial/tumor cells as "keratinocyte" (66,586 cells)
- Published data distinguishes 6 different epithelial/tumor subtypes (11,758 total cells)
- This suggests the predicted model lacks tumor-specific training data
- Many predicted cell types (melanocyte, Schwann cell, skin fibroblast) reflect skin tissue origin rather than breast biology
