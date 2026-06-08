use std::collections::HashMap;

#[derive(Debug)]
pub struct Biomarker {
    pub gene: String,
    pub cancer_type: String,
    pub rna_name: String,
    pub direction: String, // "up" or "down" in cancer
    pub clinical_significance: String,
}

pub fn get_cancer_panels() -> HashMap<String, Vec<Biomarker>> {
    let mut panels: HashMap<String, Vec<Biomarker>> = HashMap::new();

    // Breast cancer panel
    panels.insert("breast".to_string(), vec![
        Biomarker {
            gene: "ESR1".to_string(),
            cancer_type: "breast".to_string(),
            rna_name: "ENST00000342189".to_string(),
            direction: "up".to_string(),
            clinical_significance: "ER+ subtype, tamoxifen response".to_string(),
        },
        Biomarker {
            gene: "ERBB2".to_string(),
            cancer_type: "breast".to_string(),
            rna_name: "ENST00000357659".to_string(),
            direction: "up".to_string(),
            clinical_significance: "HER2+ subtype, trastuzumab target".to_string(),
        },
        Biomarker {
            gene: "MKI67".to_string(),
            cancer_type: "breast".to_string(),
            rna_name: "ENST00000379384".to_string(),
            direction: "up".to_string(),
            clinical_significance: "Proliferation marker, Ki-67 index".to_string(),
        },
        Biomarker {
            gene: "GSTM1".to_string(),
            cancer_type: "breast".to_string(),
            rna_name: "ENST00000338643".to_string(),
            direction: "down".to_string(),
            clinical_significance: "Luminal A marker".to_string(),
        },
        Biomarker {
            gene: "FGFR1".to_string(),
            cancer_type: "breast".to_string(),
            rna_name: "ENST00000262430".to_string(),
            direction: "up".to_string(),
            clinical_significance: "Amplification, poor prognosis".to_string(),
        },
    ]);

    // Lung cancer panel
    panels.insert("lung".to_string(), vec![
        Biomarker {
            gene: "EGFR".to_string(),
            cancer_type: "lung".to_string(),
            rna_name: "ENST00000257586".to_string(),
            direction: "up".to_string(),
            clinical_significance: "TKI response predictor".to_string(),
        },
        Biomarker {
            gene: "ALK".to_string(),
            cancer_type: "lung".to_string(),
            rna_name: "ENST00000379321".to_string(),
            direction: "up".to_string(),
            clinical_significance: "ALK inhibitor target".to_string(),
        },
        Biomarker {
            gene: "KRAS".to_string(),
            cancer_type: "lung".to_string(),
            rna_name: "ENST00000274472".to_string(),
            direction: "up".to_string(),
            clinical_significance: "Mutation, poor prognosis".to_string(),
        },
        Biomarker {
            gene: "TP53".to_string(),
            cancer_type: "lung".to_string(),
            rna_name: "ENST00000269305".to_string(),
            direction: "down".to_string(),
            clinical_significance: "Tumor suppressor, mutation".to_string(),
        },
    ]);

    // Prostate cancer panel
    panels.insert("prostate".to_string(), vec![
        Biomarker {
            gene: "KLK3".to_string(),
            cancer_type: "prostate".to_string(),
            rna_name: "ENST00000261782".to_string(),
            direction: "up".to_string(),
            clinical_significance: "PSA, diagnostic marker".to_string(),
        },
        Biomarker {
            gene: "TMPRSS2".to_string(),
            cancer_type: "prostate".to_string(),
            rna_name: "ENST00000389727".to_string(),
            direction: "up".to_string(),
            clinical_significance: "ERG fusion partner".to_string(),
        },
        Biomarker {
            gene: "ERGF".to_string(),
            cancer_type: "prostate".to_string(),
            rna_name: "ENST00000373736".to_string(),
            direction: "up".to_string(),
            clinical_significance: "Fusion oncogene".to_string(),
        },
    ]);

    // Colorectal cancer panel
    panels.insert("colorectal".to_string(), vec![
        Biomarker {
            gene: "APC".to_string(),
            cancer_type: "colorectal".to_string(),
            rna_name: "ENST00000246161".to_string(),
            direction: "down".to_string(),
            clinical_significance: "Wnt pathway, mutation".to_string(),
        },
        Biomarker {
            gene: "KRAS".to_string(),
            cancer_type: "colorectal".to_string(),
            rna_name: "ENST00000274472".to_string(),
            direction: "up".to_string(),
            clinical_significance: "Anti-EGFR resistance".to_string(),
        },
        Biomarker {
            gene: "BRAF".to_string(),
            cancer_type: "colorectal".to_string(),
            rna_name: "ENST00000259934".to_string(),
            direction: "up".to_string(),
            clinical_significance: "V600E, poor prognosis".to_string(),
        },
        Biomarker {
            gene: "MSH2".to_string(),
            cancer_type: "colorectal".to_string(),
            rna_name: "ENST00000347413".to_string(),
            direction: "down".to_string(),
            clinical_significance: "MMR deficiency, MSI".to_string(),
        },
    ]);

    // Pan-cancer通用 markers
    panels.insert("pancancer".to_string(), vec![
        Biomarker {
            gene: "TP53".to_string(),
            cancer_type: "pancancer".to_string(),
            rna_name: "ENST00000269305".to_string(),
            direction: "down".to_string(),
            clinical_significance: "Most common cancer mutation".to_string(),
        },
        Biomarker {
            gene: "MYC".to_string(),
            cancer_type: "pancancer".to_string(),
            rna_name: "ENST00000394617".to_string(),
            direction: "up".to_string(),
            clinical_significance: "Oncogene, amplification".to_string(),
        },
        Biomarker {
            gene: "PTEN".to_string(),
            cancer_type: "pancancer".to_string(),
            rna_name: "ENST00000262620".to_string(),
            direction: "down".to_string(),
            clinical_significance: "Tumor suppressor, PI3K pathway".to_string(),
        },
        Biomarker {
            gene: "VEGFA".to_string(),
            cancer_type: "pancancer".to_string(),
            rna_name: "ENST00000253637".to_string(),
            direction: "up".to_string(),
            clinical_significance: "Angiogenesis, anti-VEGF target".to_string(),
        },
        Biomarker {
            gene: "MDM2".to_string(),
            cancer_type: "pancancer".to_string(),
            rna_name: "ENST00000263496".to_string(),
            direction: "up".to_string(),
            clinical_significance: "p53 inhibitor, amplification".to_string(),
        },
    ]);

    panels
}

pub fn evaluate_biomarkers(
    abundances: &HashMap<String, f64>,
    panel: &str,
) -> Vec<(String, f64, String, String)> {
    let panels = get_cancer_panels();
    let mut results = Vec::new();

    if let Some(biomarkers) = panels.get(panel) {
        for bm in biomarkers {
            if let Some(&abundance) = abundances.get(&bm.rna_name) {
                results.push((
                    bm.gene.clone(),
                    abundance,
                    bm.direction.clone(),
                    bm.clinical_significance.clone(),
                ));
            }
        }
    }

    results
}
